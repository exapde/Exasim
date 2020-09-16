from numpy import *

def checkapp(app1,app2):

    diff = ones(13)*100;
    diff[0] = max(abs(app1['nsize'].flatten('F')-app2['nsize'].flatten('F')));
    diff[1] = max(abs(app1['ndims'].flatten('F')-app2['ndims'].flatten('F')));
    diff[2] = max(abs(app1['flag'].flatten('F')-app2['flag'].flatten('F')));
    diff[3] = max(abs(app1['problem'].flatten('F')-app2['problem'].flatten('F')));
    diff[4] = max(abs(app1['uinf'].flatten('F')-app2['uinf'].flatten('F')));
    diff[5] = max(abs(app1['dt'].flatten('F')-app2['dt'].flatten('F')));
    diff[6] = max(abs(app1['factor'].flatten('F')-app2['factor'].flatten('F')));
    diff[7] = max(abs(app1['physicsparam'].flatten('F')-app2['physicsparam'].flatten('F')));
    diff[8] = max(abs(app1['solversparam'].flatten('F')-app2['solversparam'].flatten('F')));
    diff[9] = max(abs(app1['tau'].flatten('F')-app2['tau'].flatten('F')));
    diff[10] = max(abs(app1['stgdata'].flatten('F')-app2['stgdata'].flatten('F')));
    diff[11] = max(abs(app1['stgparam'].flatten('F')-app2['stgparam'].flatten('F')));
    diff[12] = max(abs(app1['stgib'].flatten('F')-app2['stgib'].flatten('F')));

    return diff
