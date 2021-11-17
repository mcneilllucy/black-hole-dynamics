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

path = '/fred/oz038/mm_nbody6/200k/Z200k/rdv_analysis/'

f = h5py.File('/fred/oz038/mm_nbody6/200k/Z200k/rdv_analysis/snapdata.hdf5', 'r')
print('loaded the hdf5 file')

testing_times = np.loadtxt('/home/lmcneill/nbody/batch_output/200k/Z200/times_for_200k_Z200.txt')



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





def time_avg(array,xwindow):
    """gets the spatial average at the current timestep.

    Parameters
    ----------
    array
        data which needs spatial averaging.
    xwindow
        size of window.


    Returns
    -------
    new_array_t
        spatially smoothed curve
    """
    twindow = len(array[1,:])
    new_array=np.zeros((len(array),twindow))
    new_array_t = []

    t_i = 0
    while t_i<twindow:
        i = 0
        while i<len(array):
            if i<int(xwindow/2):
                averaged = np.mean(array[0:i+int(xwindow/2),t_i])
                new_array[i,t_i] = averaged
                i+=1
            elif i<len(array)-int(xwindow/2):
                averaged = np.mean(array[i-int(xwindow/2):i+int(xwindow/2),t_i])
                new_array[i,t_i] = averaged
                i+=1
            else:
                averaged = np.mean(array[i-int(xwindow/2):len(array),t_i])
                new_array[i,t_i] = averaged
                i+=1
        t_i+=1
    #### now take average of all your arrays
    i = 0
    while i<len(array):
        new_array_t.append(np.mean(new_array[i,:]))
        i+=1

    return(new_array_t)


def get_densities_all(path,f, hdf_index):
    """Calculate the number and mass of objects around the centre of mass fro a fixed bubble size.

    Parameters
    ----------
    path
        Location of hdf5 file from NBODY.
    f
        hdf5 file


    Returns
    -------
    arrays of densities

        The mass and number density


    """
    ## Calculate the local mass and number densities. Just of the BH population.

## do we want to save these as output? To compute power laws etc? Sure, if it's just a 2d array.
## save local density array
    hdf_index = testing_times
    #bubble_arr = np.logspace(-1.5,0.5,50)
    bubble_arr = np.logspace(-1.5,-0.5,10)
    dens_time_array = np.zeros((len(hdf_index),len(bubble_arr)))
    density_arr = np.zeros((len(hdf_index),len(bubble_arr),2))
    time_idx = 0
    snap_list = list(f)
    snap_len = len(snap_list)
    for times in hdf_index:
        data = f.get(snap_list[int(times)])
        print(int(data['t'].value),str('Myr'))
    ## identify black holes
        jj=0
        data = f.get(snap_list[int(times)])
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

    np.savetxt('/home/lmcneill/nbody/batch_output/200k/Z200k/density_M_all.txt',density_arr[:,:,0])
    np.savetxt('/home/lmcneill/nbody/batch_output/200k/Z200k/density_N_all.txt',density_arr[:,:,1])

    ## load the data files


def get_densities_BH_only(path,f,hdf_index):
    """Calculate the number and mass of objects around the centre of mass fro a fixed bubble size.

    Parameters
    ----------
    path
        Location of hdf5 file from NBODY.
    f
        hdf5 file


    Returns
    -------
    arrays of densities

        The mass and number density


    """

    bubble_arr_BH = np.logspace(-1.5,-0.5,50)
    density_arr_BH = np.zeros((len(hdf_index),len(bubble_arr_BH),2))
    time_idx = 0
    snap_list = list(f)
    snap_len = len(snap_list)
    for times in hdf_index:
        data = f.get(snap_list[int(times)])
        print(int(data['t'].value),str('Myr'))
        ## identify black holes
        jj=0
        data = f.get(snap_list[int(times)])
        #### locate each black hole, contruct a bubble around it..
        mask = ((data['kw'].value[:] == 14))
        selected = data['kw'].value[mask]
        N_start_BH = np.array(data['m'].value[mask])
        N_start = len(np.array(data['m'].value[mask]))
        N_start2 = N_start
        N_start_all = len(np.array(data['m'].value[:]))
        dens_rad_BH = []
        mass_rad_BH = []
        num_rad_BH = []
        #### add array which identifies the BHs to later overwrite
        while jj<len(bubble_arr_BH):
            #densities.append(data['t'].value)
            ###### make an array with all x y z positions
            bubble_size = bubble_arr_BH[jj]
            all_pos = np.zeros((N_start_all, 3))
            all_pos[:,0] = data['x'].value[:]
            all_pos[:,1] = data['y'].value[:]
            all_pos[:,2] = data['z'].value[:]
            mass = data['m'].value[:]
            ### centre of mass
            CM = np.average(all_pos, axis=0, weights=mass)
            ### now just BH things
            all_pos_BH = np.zeros((N_start, 3))
            all_pos_BH[:,0] = data['x'].value[mask]
            all_pos_BH[:,1] = data['y'].value[mask]
            all_pos_BH[:,2] = data['z'].value[mask]
            mass_BH = data['m'].value[mask]
            mass_tot,n_tot = calculate_N_M_radius(all_pos_BH, mass_BH,bubble_size,CM)
            #dens_rad_BH.append(calculate_density_radius(all_pos, mass,bubble_size,CM))
            mass_rad_BH.append(mass_tot)
            num_rad_BH.append(n_tot)
            jj = jj+1
        #### now append density and mass arrays to final one

        ### compute derivative here
        ## rho = dM/dr /(4 pi r^2)
        mass_dens_BH = np.gradient(mass_rad_BH,bubble_arr_BH)/(4*np.pi*(bubble_arr_BH)**2)
        num_dens_BH = np.gradient(num_rad_BH,bubble_arr_BH)/(4*np.pi*(bubble_arr_BH)**2)
        density_arr_BH[time_idx,:,0]= mass_dens_BH
        density_arr_BH[time_idx,:,1]= num_dens_BH
        time_idx+=1

    np.savetxt('/home/lmcneill/nbody/batch_output/200k/Z200k/density_M_BH.txt',density_arr_BH[:,:,0])
    np.savetxt('/home/lmcneill/nbody/batch_output/200k/Z200k/density_N_BH.txt',density_arr_BH[:,:,1])

def get_densities_non_BH_only(path,f,hdf_index):
    """Calculate the number and mass of objects around the centre of mass fro a fixed bubble size.

    Parameters
    ----------
    path
        Location of hdf5 file from NBODY.
    f
        hdf5 file


    Returns
    -------
    arrays of densities

        The mass and number density


    """
    #bubble_arr_nonBH = np.logspace(-1.5,0.5,50)
    bubble_arr_nonBH = np.logspace(-1.5,-0.5,10)
    density_arr_nonBH = np.zeros((len(hdf_index),len(bubble_arr_nonBH),2))
    time_idx = 0
    snap_list = list(f)
    snap_len = len(snap_list)
    for times in hdf_index:
        data = f.get(snap_list[int(times)])
        print(int(data['t'].value),str('Myr'))
    ## identify black holes
        jj=0
        data = f.get(snap_list[int(times)])
    #### locate each black hole, contruct a bubble around it..
        mask = ((data['kw'].value[:] != 14))
        selected = data['kw'].value[mask]
        N_start_nonBH = np.array(data['m'].value[mask])
        N_start = len(np.array(data['m'].value[mask]))
        N_start2 = N_start
        dens_rad_nonBH = []
        mass_rad_nonBH = []
        num_rad_nonBH = []
        #### add array which identifies the BHs to later overwrite
        while jj<len(bubble_arr_nonBH):
            #densities.append(data['t'].value)
            ###### make an array with all x y z positions
            bubble_size = bubble_arr_nonBH[jj]
            all_pos = np.zeros((N_start, 3))
            all_pos[:,0] = data['x'].value[mask]
            all_pos[:,1] = data['y'].value[mask]
            all_pos[:,2] = data['z'].value[mask]
            mass = data['m'].value[mask]
            ### centre of mass
            CM = np.average(all_pos, axis=0, weights=mass)
            mass_tot,n_tot = calculate_N_M_radius(all_pos, mass,bubble_size,CM)
            #dens_rad_nonBH.append(calculate_density_radius(all_pos, mass,bubble_size,CM))
            mass_rad_nonBH.append(mass_tot)
            num_rad_nonBH.append(n_tot)
            jj = jj+1
        #### now append density and mass arrays to final one

        ### compute derivative here
        ## rho = dM/dr /(4 pi r^2)
        mass_dens_nonBH = np.gradient(mass_rad_nonBH,bubble_arr_nonBH)/(4*np.pi*(bubble_arr_nonBH)**2)
        num_dens_nonBH = np.gradient(num_rad_nonBH,bubble_arr_nonBH)/(4*np.pi*(bubble_arr_nonBH)**2)
        density_arr_nonBH[time_idx,:,0]= mass_dens_nonBH
        density_arr_nonBH[time_idx,:,1]= num_dens_nonBH
        time_idx+=1
    np.savetxt('/home/lmcneill/nbody/batch_output/200k/Z200k/density_M_non_BH.txt',density_arr_nonBH[:,:,0])
    np.savetxt('/home/lmcneill/nbody/batch_output/200k/Z200k/density_N_non_BH.txt',density_arr_nonBH[:,:,1])


get_densities_all(path,f,testing_times)
get_densities_BH_only(path,f,testing_times)
get_densities_non_BH_only(path,f,testing_times)
