# import packages
import re,csv
import numpy as np
import h5py
import imageio
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib import gridspec
from time import time
import numpy as np
from scipy.spatial import cKDTree
import glob, os
# can200k  fb10200k  IMBH200k  kick200k  Z200k
path = '/fred/oz038/mm_nbody6/200k/fb10200k/rdv_analysis/'
tstart = 500
tend = 19000
gap = 500
## n around is how many more data points either side
n_around = 20
# ngap is how many Myr between steps
n_gap = 9
n_space = n_around*2+1

f = h5py.File('/fred/oz038/mm_nbody6/200k/fb10200k/rdv_analysis/snapdata.hdf5', 'r')
print('loaded the hdf5 file')


def get_time_index(path,f):
    """Calculate the number and mass of objects around the centre of mass fro a fixed bubble size.

    Parameters
    ----------
    path
        Location of hdf5 file from NBODY.
    f
        hdf5 file


    Returns
    -------
    hdf_index

        The index in hdf5 file corresponding to desired time array.


    """


    ## load the data files



    snap_list = list(f)
    snap_len = len(snap_list)

    step_num = int((tend-tstart)/gap)+1
    time_list_escape = np.linspace(np.float(tstart+1), tend+1, step_num)

    ### get the timesteps of interest
    times_tmp = []


    ### now do a for loop for every single time in the current list.
    for item in time_list_escape:
        around = np.linspace(item-n_around*n_gap,item+n_around*n_gap,n_space)
        for value in around:
            times_tmp.append(value)

    orig_times = time_list_escape
    time_list_escape = times_tmp
    ### create your special time list, then loop over that (only ~10 elements anyway)
    time_len = len(time_list_escape)
    hdf_index = []
    j = 0
    i = 0
    while j < time_len:
        ### create your special time list, then loop over that (only ~10 elements anyway)
        ### set a distance for the closest time
        dist = 10
        i = 0
        while i < snap_len:
            data = f.get(snap_list[i])
            if abs(float(data['t'].value)-float(time_list_escape[j]))<dist:

                ### if statement to check that our m1 and m2 will even be in the output file (or ejected)
                if float(data['t'].value)-float(time_list_escape[j])<0:
                    dist = abs(float(data['t'].value)-float(time_list_escape[j]))
                    snap_index = i
                else:
                    foundit = 0
                    while foundit == 0:
                        dist = abs(float(data['t'].value)-float(time_list_escape[j]))
                        if dist > 50:
                            snap_index = i+1
                            data = f.get(snap_list[snap_index])
                        else:
                            foundit = 1
        ### find closest value to special time and access properties of this
            i +=1
        ### here you want to check if snap_index
        hdf_index.append(snap_index)
        data = f.get(snap_list[hdf_index[j]])
        ### here is where you can check
        print("printing time")
        print(data['t'].value)
        j+=1
    # have hdf_index from using time_list_escape
    out_path = '/home/lmcneill/nbody/batch_output/200k/fb10200k/times_for_200k_fb10200k.txt'
    out_path2 = '/home/lmcneill/nbody/batch_output/200k/fb10200k/orig_times_for_200k_fb10200k.txt'
    ### save hdf_index for troubleshooting.
    np.savetxt(out_path,hdf_index)
    np.savetxt(out_path2,orig_times)
    return hdf_index


hdf_indices = get_time_index(path,f)
