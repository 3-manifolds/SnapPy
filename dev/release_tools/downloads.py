#! /usr/bin/env python3

import os, sys, json, requests

class GithubREST(object):
    """
    Interacts with Github via their REST API.
    """
    
    def __init__(self, project):
        self.project = project
        self.json_decoder = json.JSONDecoder()
        
    def json_info(self):
        """
        Return the download counts for all assets in the release.
        """
        template = 'https://api.github.com/repos/3-manifolds/{0.project}/releases'
        url = template.format(self)
        response = requests.get(url)
        return self.json_decoder.decode(response.content.decode('utf-8'))

def main():
    if (len(sys.argv) == 1 or
        set(['help', '-help', '--help', '-h']).intersection(sys.argv)):
        print('usage: python downloads.py project')
        sys.exit()

    print('Download Counts for %s:'%sys.argv[1])
    for info_dict in GithubREST(sys.argv[1]).json_info():
        print('    Release: %s'%info_dict['tag_name'])
        assets = info_dict['assets']
        for asset in assets:
            print('    %5d    %s'%(asset['download_count'], asset['name']))
            
if __name__ == '__main__':
    main()
