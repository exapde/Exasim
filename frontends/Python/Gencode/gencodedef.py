import pdeapp
import numpy as np

#def gencode(app):

param = np.array([1.0]);
xdg = np.array([0.0, 0.0]);
udg = np.array([2.0, 2.0, 2.0]);
odg = np.array([0.0]);
wdg = np.array([0.0]);
uinf = np.array([0.0]);
time = np.array([0.0]);
print(pdeapp.flux(xdg, udg, odg, wdg, uinf, param, time))
    #Flux = getattr(pdeapp, app['Flux'])
    #print(Flux(xdg, udg, odg, wdg, uinf, param, time))

#    return 0;
