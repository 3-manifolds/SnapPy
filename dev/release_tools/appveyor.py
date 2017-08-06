#! /usr/bin/env python

import os, sys, json, requests
from future.builtins import input
from configparser import ConfigParser

class AppveyorREST(object):
    """
    Interacts with Appveyor via their REST API.
    """
    host = 'ci.appveyor.com'
    
    def __init__(self, project, username, api_token):
        self.project, self.username, self.api_token = project, username, api_token
        self.json_decoder = json.JSONDecoder()
        
    def jobs(self):
        """
        Return the list of job ids from the latest build.
        """
        template = 'https://ci.appveyor.com/api/projects/{0.username}/{0.project}'
        url = template.format(self)
        headers = {'Authorization': 'Bearer {0.api_token}'.format(self)}
        response = requests.get(url, headers=headers)
        info = self.json_decoder.decode(response.content.decode('utf-8'))
        return [job['jobId'] for job in info['build']['jobs']]

    def artifacts(self, jobId):
        """
        Return a list of the artifact filenames from the latest build.
        """
        template = 'https://ci.appveyor.com/api/buildjobs/{job}/artifacts'
        url = template.format(job=jobId)
        headers = {'Authorization': 'Bearer {0.api_token}'.format(self)}
        response = requests.get(url, headers=headers)
        arts = self.json_decoder.decode(response.content.decode('utf-8'))
        return [art['fileName'] for art in arts if art['type'] == 'File']
        
    def download(self, dir):
        """
        Download all of the artifacts from the latest build.
        """
        if not os.path.exists(dir):
            print('Directory {dir} does not exist.'.format(dir=dir))
            return
        template = 'https://ci.appveyor.com/api/buildjobs/{jobId}/artifacts/{path}'
        jobIds = self.jobs()
        for jobId in jobIds:
            paths = self.artifacts(jobId)
            for path in paths:
                url = template.format(jobId=jobId, path=path)
                headers = {'Authorization': 'Bearer {0.api_token}'.format(self)}
                response = requests.get(url, headers=headers)
                filename = os.path.split(path)[-1] 
                with open(os.path.join(dir, filename), 'wb') as output:
                    output.write(response.content)
                print(filename)

    def clear_cache(self):
        """
        Clear the build cache.
        """
        template = 'https://ci.appveyor.com/api/projects/{0.username}/{0.project}/buildcache'
        url = template.format(self)
        headers = {'Authorization': 'Bearer {0.api_token}'.format(self)}
        response = requests.delete(url, headers=headers)
        if response.status_code == requests.code.ok:
            print('Cache deleted.')
        else:
            response.raise_for_status()

    def run_build(self):
        """
        Start a new build
        """
        url = 'https://ci.appveyor.com/api/builds'
        headers = {'Authorization': 'Bearer {0.api_token}'.format(self),
                   'Content-type': 'application/json'}
        data = {
            'accountName': self.username,
            'projectSlug': self.project,
            'branch': 'default',
            'environmentVariables': {}
        }
        response = requests.post(url, data=json.dumps(data), headers=headers)
        response.raise_for_status()
        info = self.json_decoder.decode(response.content.decode('utf-8'))
        print('Started build#{build_id}.'.format(build_id=info['buildId']))

def main():
    if len(sys.argv) == 1 or set(['help', '-help', '--help', '-h']).intersection(sys.argv):
        print('usage: python appveryor.py [help|build|clear|download <to_dir>]')
        sys.exit()
        
    config = ConfigParser()
    if os.path.exists('appveyor.cfg'):
        config.read('appveyor.cfg')
    else:
        print('appveyor.cfg not found - creating ...')
        config.add_section('appveyor')
        config.set('appveyor', 'project', input('Project name: '))
        config.set('appveyor', 'username', input('Username: '))
        config.set('appveyor', 'token', input('Api token: '))
        with open('appveyor.cfg', 'wb') as output:
            config.write(output)
            
    A = AppveyorREST(config.get('appveyor', 'project'),
                     config.get('appveyor', 'username'),
                     config.get('appveyor', 'token'))

    skip = False
    for n, command in enumerate(sys.argv[1:]):
        if skip:
            skip = False
            continue
        if command == 'build':
            A.run_build()
        elif command == 'clear':
            A.clear_cache()
        elif command == 'download':
            try:
                dir = sys.argv[n+2]
            except IndexError:
                print('Please supply a download directory.')
                sys.exit()
            A.download(dir)
            skip = True
        else:
            print('unknown command {cmd}'.format(cmd=command))
            sys.exit()
            
if __name__ == '__main__':
    main()
