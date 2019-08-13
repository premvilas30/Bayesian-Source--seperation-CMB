import subprocess
import os,json
from os import listdir
from os.path import isfile, join
import glob
import json
import subprocess
import os
import matplotlib.pyplot as plt
import numpy as np
import pyfits as pf
import healpy as hp
from matplotlib import cm
import matplotlib
matplotlib.use('Agg')

mylist = [f for f in glob.glob("../WMAP_result/*.fits")]

#print(mylist)
input = input("Press 1 for old output, and 2 for new execution")

cmd = "./test.out "
with open('config.json') as config_file:
    data = json.load(config_file)
    for i in data['input']['Files']:
        cmd = cmd+' "'+data["input"]["source_location"]+i+'"'
#print(cmd)
mypath= '../WMAP_plot/'
#print([f for f in listdir(mypath) if isfile(join(mypath, f))])

if str(input) == "1":
     os.system(cmd)
elif str(input) == "2":
     if os.system('make wmap_analysis_9yr') == 0:
        os.system(cmd)
     else:
        print("Comilation error")
else:
    print("wrong input")

inp_plots = input(" do you want to create maps of results (Y or N)")
cmap = cm.bwr
cmap = cm.binary
cmap.set_under('green')
cmap.set_bad('gray')
for i in mylist:
    print("Saving plot for :",i)
    path = mypath+str(i.split('/')[-1])[:-5]+'.pdf'
    #print(path)
    cmap = cm.bwr
    cmap = cm.binary
    cmap.set_under('green')
    cmap.set_bad('gray')
    mapin = hp.read_map(i)
    hp.mollview(mapin, cmap=cmap)
    plt.savefig(path)
    plt.close("all")

#if str(input) == "1":
#     os.system(cmd)
#elif str(input) == "2":
#     if os.system('make wmap_analysis_9yr') == 0:
#        os.system('./test.out "../WMAP_data/wmap_band_imap_r9_9yr_K_v5.fits" "../WMAP_data/wmap_band_imap_r9_9yr_Ka_v5.fits" "../WMAP_data/wmap_band_imap_r9_9yr_Q_v5.fits" "../WMAP_data/wmap_band_imap_r9_9yr_V_v5.fits" "../WMAP_data/wmap_band_imap_r9_9yr_W_v5.fits"')
#     else:
#        print("Comilation error")
#else:
#    print("wrong input")

