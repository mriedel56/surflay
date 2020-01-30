
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

tmp = np.ones([256,4])
tmp[:,1]  = tmp[:,1]*np.linspace(0,1,256)
magenta_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[128:256,0] = tmp[128:256,0]*np.linspace(1,0.5,128)
tmp[:,1]  =tmp[:,1]*0
tmp[128:256,2] = tmp[128:256,2]*np.linspace(1,0.5,128)
magentablack_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[128:256,0]  = tmp[128:256,0]*np.linspace(1,0,128)
tmp[128:256,1]  = tmp[128:256,1]*np.linspace(1,0,128)
whiteblue_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[128:256,1]  = tmp[128:256,1]*np.linspace(1,0,128)
whitemagenta_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[128:256,2]  = tmp[128:256,1]*np.linspace(1,0,128)
whiteyellow_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[128:256,0]  = tmp[128:256,1]*np.linspace(1,0,128)
whiteaqua_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[128:256,1]  = tmp[128:256,1]*np.linspace(1,0,128)
tmp[128:256,2]  = tmp[128:256,2]*np.linspace(1,0,128)
whitered_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[128:256,0]  = tmp[128:256,0]*np.linspace(1,0,128)
tmp[128:256,2]  = tmp[128:256,2]*np.linspace(1,0,128)
whitegreen_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[0:128,0] = tmp[0:128,0]*np.linspace(0,1,128)
tmp[128:256,1] = tmp[128:256,1]*np.linspace(1,0,128)
tmp[:,2] = tmp[:,2]*0
greenyellowred_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[0:128,0] = tmp[0:128,0]*np.linspace(0,1,128)
tmp[128:256,1] = tmp[128:256,1]*np.linspace(1,0,128)
tmp[0:128,2] = tmp[0:128,2]*np.linspace(0,1,128)
tmp[128:256,2] = tmp[128:256,2]*np.linspace(1,0,128)
greenwhitered_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[0:128,0] = tmp[0:128,0]*np.linspace(0,1,128)
tmp[128:256,0] = tmp[128:256,0]*np.linspace(1,0,128)
tmp[128:256,1] = tmp[128:256,1]*np.linspace(1,0,128)
tmp[0:128,2] = tmp[0:128,2]*np.linspace(0,1,128)
greenwhiteblue_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[:,0] = tmp[:,0]*0
tmp[128:256,1] = tmp[128:256,1]*np.linspace(1,0,128)
tmp[0:128,2] = tmp[0:128,2]*np.linspace(0,1,128)
greenaquablue_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[0:128,0] = tmp[0:128,0]*np.linspace(0,1,128)
tmp[0:128,1] = tmp[0:128,1]*np.linspace(0,1,128)
tmp[128:256,1] = tmp[128:256,1]*np.linspace(1,0,128)
tmp[128:256,2] = tmp[128:256,2]*np.linspace(1,0,128)
bluewhitered_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[128:256,0] = tmp[128:256,0]*np.linspace(1,0,128)
tmp[0:128,1] = tmp[0:128,1]*np.linspace(0,1,128)
tmp[128:256,1] = tmp[128:256,1]*np.linspace(1,0,128)
tmp[0:128,2] = tmp[0:128,2]*np.linspace(0,1,128)
redwhiteblue_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[128:256,0] = tmp[128:256,0]*np.linspace(1,0,128)
tmp[0:128,1] = tmp[0:128,1]*np.linspace(0,1,128)
tmp[0:128,2] = tmp[0:128,2]*np.linspace(0,1,128)
tmp[128:256,2] = tmp[128:256,2]*np.linspace(1,0,128)
redwhitegreen_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[0:128,0] = tmp[0:128,0]*np.linspace(0,1,128)
tmp[:,1] = tmp[:,1]*0
tmp[128:256,2] = tmp[128:256,2]*np.linspace(1,0,128)
bluemagentared_cm = ListedColormap(tmp)

tmp = np.ones([256,4])
tmp[0:86,0] = tmp[0:86,0]*0
tmp[:,1] = tmp[:,1]*0
tmp[172:256,2] = tmp[172:256,2]*0
bluemagentarednospect_cm = ListedColormap(tmp)

tmp = cm.get_cmap('jet', 256)
newred = tmp._segmentdata['red']
newred = list(newred)
for i in range(len(newred)):
    tmpcol = list(newred[i])
    tmpcol[0] = tmpcol[0]/2+0.5
    newred[i] = tuple(tmpcol)
newred = tuple([tuple([0,0,0])] + newred)

newgreen = tmp._segmentdata['green']
newgreen = list(newgreen)
for i in range(len(newgreen)):
    tmpcol = list(newgreen[i])
    tmpcol[0] = tmpcol[0]/2+0.5
    newgreen[i] = tuple(tmpcol)
newgreen = tuple([tuple([0,0,0])] + newgreen)

newblue = tmp._segmentdata['blue']
newblue = list(newblue)
for i in range(len(newblue)):
    tmpcol = list(newblue[i])
    tmpcol[0] = tmpcol[0]/2+0.5
    newblue[i] = tuple(tmpcol)
newblue = tuple([tuple([0,0,0])] + newblue)

cdict = {'red': newred, 'green': newgreen, 'blue': newblue}
spectrum_cm = LinearSegmentedColormap('spectrum', segmentdata=cdict, N=256)

color_dict = {'magenta': magenta_cm, 'magentablack': magentablack_cm, 'whitered': whitered_cm, 'whitemagenta': whitemagenta_cm, 'whiteyellow': whiteyellow_cm, 'whiteaqua': whiteaqua_cm, 'whiteblue': whiteblue_cm, 'whitegreen': whitegreen_cm, 'greenwhitered': greenwhitered_cm, 'greenyellowred': greenyellowred_cm, 'greenwhiteblue': greenwhiteblue_cm, 'greenaquablue': greenaquablue_cm, 'redwhiteblue': redwhiteblue_cm, 'redwhitegreen': redwhitegreen_cm, 'bluewhitered': bluewhitered_cm, 'bluemagentared': bluemagentared_cm, 'bluemagentarednospect': bluemagentarednospect_cm, 'spectrum': spectrum_cm}

