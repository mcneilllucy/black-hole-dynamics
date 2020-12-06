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

parent_dir = 'data/200k'

plt.rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 40})

# import data

# these are summary files of important events.
time_list = []
m1_list = []
m2_list = []
id1_ej = []
id2_ej = []
id1_list = []
id2_list = []
P_list = []
time_list = []
id3_list = []
id4_list = []
time_list_escape = []
id1_list_escape = []
id2_list_escape = []

event_list = []

### get time arrays. Also .csv with a summary.

## outputs: time_list_escape (from this summary file)

filenames = []

for file in glob.glob(os.path.join(parent_dir, '*_output.dat')):
    filenames.append(file)
filenames = sorted(filenames)
print(filenames)

for name in filenames:
    file = open(name,'r',errors = 'replace')
    lines_list = file.readlines()

    for line in lines_list:
        if line[10:14] == 'TIME':
            time = line[16:24]
        if line[1+7:7+7] == 'ESCAPE':

            if line[45+7:47+7] == '14':
                if line[48+7:50+7] == '14':
                    m1_list.append(line[62+7:66+7])
                    m2_list.append(line[67+7:72+7])
                    id1_list.append(line[26+7:32+7])
                    id2_list.append(line[32+7:39+7])
                    id1_ej.append(line[26+7:32+7])
                    id2_ej.append(line[32+7:39+7])
                    id3_list.append('N/A')
                    id4_list.append('N/A')
                    P_list.append(line[143+7:150+7])
                    time_list.append(time)
                    event_list.append('BH ESCAPE')
                    time_list_escape.append(time)
                    id1_list_escape.append(line[26+7:32+7])
                    id2_list_escape.append(line[32+7:39+7])


        if line[1:9] == 'EXCHANGE':
            m1_list.append('?')
            m2_list.append('?')
            time_list.append(time)
            id1_list.append(line[45:51])
            id2_list.append(line[51:57])
            id3_list.append(line[57:63])
            id4_list.append(line[63:69])
            P_list.append('?')
            event_list.append('EXCHANGE')

file.close()

### now go into large hdf5 file to locate the index of the important times which we will use.

f = h5py.File('../snapdata.hdf5', 'r')


snap_list = list(f)

snap_len = len(snap_list)


### create your special time list, then loop over that (only ~10 elements anyway)
####
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
                dist = abs(float(data['t'].value)-float(time_list_escape[j]))
                snap_index = i-1

    ### find closest value to special time and access properties of this
        i +=1
    hdf_index.append(snap_index)
    data = f.get(snap_list[hdf_index[j]])
    j+=1

# have hdf_index from using time_list_escape

## write .csv
## use this .csv to determine the properties of three body interactions and ejections
## see e.g. nbody_finding_m3
i = 0
with open('data/200k/BH_eject_exchange_properties.csv', 'w') as f:
    w = csv.writer(f, delimiter=',')
    w.writerow(['event type', 'mass 1', 'mass 2', 'id 1', 'id 2','id 3', 'id 4', 'period (days)','time (Myr)'])

    while i<len(m1_list):
            w.writerow([event_list[i],m1_list[i], m2_list[i], id1_list[i], id2_list[i], id3_list[i], id4_list[i],P_list[i],time_list[i]])
            i+=1

f.close()

## don't change anything about BH_time without changing snap_list

### now into hdf5 file
f = h5py.File('../snapdata.hdf5', 'r')
BH_time = []
for times in hdf_index:
    data = f.get(snap_list[times])
    BH_time.append(data['t'].value)

    ### want both M(r) for the radial density profile, and also N(m)?

### for each time step, we want N as function of m. N as a function of r. and  M as a function of r.

def calculate_N_M_radius(xyz, mass, bubble_size,COM):
    """Calculate the number and mass of objects around the centre of mass fro a fixed bubble size.

    Parameters
    ----------
    xyz
        Position of particles.
    mass
        Mass of particles.

    bubble_size
        The size of bubbles around each black hole.

    Returns
    -------
    mass_tot
        The stellar mass inside a bubble with radius bubble_size.

    num_tot
        The number of stars inside a bubble with radius bubble_size.
    """

    # Size of bubbles around particles
    bubble_volume = 4 * np.pi / 3 * bubble_size ** 3


    # Build kd-tree for all particles
    tree = cKDTree(xyz)

    # Get stellar neighbours of every black hole
    neighbours = tree.query_ball_point((COM), bubble_size)
    mass_tot = 0
    n_tot = 0

    for neigh in neighbours:
        mass_tot += (mass[neigh])
        n_tot+=1
    return [mass_tot,n_tot]

bubble_arr = np.logspace(-2,1,50)
dens_time_array = np.zeros((len(hdf_index),len(bubble_arr)))
density_arr = np.zeros((len(hdf_index),len(bubble_arr),2))
time_idx = 0
for times in hdf_index:
    data = f.get(snap_list[times])
    print(int(data['t'].value),str('Myr'))
## identify black holes
    jj=0
    data = f.get(snap_list[times])
#### locate each black hole, contruct a bubble around it..
    N_start = len(np.array(data['m'].value[:]))
    N_start2 = N_start
    dens_rad = []
    mass_rad = []
    num_rad = []
    #### add array which identifies the BHs to later overwrite
    while jj<len(bubble_arr):
        #densities.append(data['t'].value)
        ###### make an array with all x y z positions
        bubble_size = bubble_arr[jj]
        all_pos = np.zeros((N_start, 3))
        all_pos[:,0] = data['x'].value[:]
        all_pos[:,1] = data['y'].value[:]
        all_pos[:,2] = data['z'].value[:]
        mass = data['m'].value[:]
        ### centre of mass
        CM = np.average(all_pos, axis=0, weights=mass)
        mass_tot,n_tot = calculate_N_M_radius(all_pos, mass,bubble_size,CM)
        #dens_rad.append(calculate_density_radius(all_pos, mass,bubble_size,CM))
        mass_rad.append(mass_tot)
        num_rad.append(n_tot)
        jj = jj+1
    #### now append density and mass arrays to final one

    ### compute derivative here
    ## rho = dM/dr /(4 pi r^2)
    mass_dens = np.gradient(mass_rad,bubble_arr)/(4*np.pi*(bubble_arr)**2)
    num_dens = np.gradient(num_rad,bubble_arr)/(4*np.pi*(bubble_arr)**2)
    density_arr[time_idx,:,0]= mass_dens
    density_arr[time_idx,:,1]= num_dens
    time_idx+=1


### plot mass vs local density.
plt_idx = 0
while plt_idx<len(BH_time):
    time_BH = BH_time[plt_idx]
    print(time_BH)
    densities = density_arr[plt_idx,:,1]
    fig = plt.figure()
    axes = plt.gca()
    line1 = plt.plot(3.10*bubble_arr,np.array(densities)/3.10**3, color='black', linestyle='-', linewidth=5)
    fig.set_size_inches(15, 10)
    axes.set_xlabel('radius (pc)')
    axes.set_ylabel('number density (/pc$^3$)')
    axes.set_xscale('log')
    #axes.set_yscale('log')
    plt.grid(True)
    plt.rcParams.update({'font.size': 30})
    axes.set_xlim([0.01*3, 13])
    axes.set_ylim([0, 7000])
    plt.legend(('t='+str(int(time_BH))+' Myr','Entropy'), loc='upper right')
    fig.savefig('num_dens_200k_can_radius'+str(int(time_BH))+'.png')
    plt_idx+=1

plt_idx=0
images = []
filenames = []
while plt_idx<len(BH_time):
    time_BH = BH_time[plt_idx]
    filenames.append('num_dens_200k_can_radius'+str(int(time_BH))+'.png')
    plt_idx+=1
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave('num-densities-radius-200k-can.gif', images,fps=2)


scaling = 127918.2*(3.10)**(-3)
### plot mass vs local density.
plt_idx = 0
while plt_idx<len(BH_time):
    time_BH = BH_time[plt_idx]
    print(time_BH)
    densities = density_arr[plt_idx,:,0]



    fig = plt.figure()
    axes = plt.gca()
    line1 = plt.plot(3.10*bubble_arr,scaling*np.array(densities), color='black', linestyle='-', linewidth=5)
    fig.set_size_inches(15, 10)
    axes.set_xlabel('radius (pc)')
    axes.set_ylabel('Density ($M_\odot/$pc$^3$)')
    axes.set_xscale('log')
    #axes.set_yscale('log')
    plt.grid(True)
    plt.rcParams.update({'font.size': 30})
    axes.set_xlim([0.01*3, 13])
    axes.set_ylim([0, 7000])
    plt.legend(('t='+str(int(time_BH))+' Myr','Entropy'), loc='upper right')
    fig.savefig('mass_dens_200k_can_radius'+str(int(time_BH))+'.png')
    plt_idx+=1

plt_idx=0
images = []
filenames = []
while plt_idx<len(BH_time):
    time_BH = BH_time[plt_idx]
    filenames.append('mass_dens_200k_can_radius'+str(int(time_BH))+'.png')
    plt_idx+=1
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave('mass-densities-radius-200k-can.gif', images,fps=2)
