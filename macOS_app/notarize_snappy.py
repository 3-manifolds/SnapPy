import release
from notabot import Notarizer

class SnapPyNotarizer(Notarizer):

    def build_dmg(self):
        app_name = self.config['app']['app_name']
        dmg_file = self.config['app']['dmg_path']
        print('building dmg %s for %s'%(dmg_file, app_name))
        release.package_app(app_name)

if __name__ == '__main__':
    notarizer = SnapPyNotarizer('notabot.cfg')
    notarizer.run()
        
