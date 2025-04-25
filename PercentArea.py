import numpy as np
#import healpy as hp #this does skymap plotting (based on matplotlib)
#from healpy import nside2pixarea
import pandas as pd

from astropy.io import fits
from astropy.table import Column, Table, join

from ligo.skymap.moc import rasterize #this processes skymaps 
from ligo.skymap.io.fits import read_sky_map #this reads skymaps

#for i in range(len(O4a_superevents)):
def get_area(file_location, event_array, area_percent=90, onepixelarea=0.052455852825697924):
    #file_location = '/home/sumer/Downloads/Problem2_folder_copy/O4a_fits_fornotes/'
    superevent_fits=file_location+event_array[0]+"_fits,"+str(event_array[1])
    M1_data=read_sky_map(superevent_fits, nest=False, distances=True, moc=True)
    M1_nested_data=rasterize(M1_data,order=8)
        #when rasterize: nside = 2^order 
        #number of pixels = 12*(n_side^2) 
        #order shouldn't go more than 9 since we want it fast. 
    #onepixelarea=nside2pixarea(256, degrees = True)  
    M1_probdensity = M1_nested_data ['PROBDENSITY']
    M1_mu=M1_nested_data['DISTMU']
    M1_sigma=M1_nested_data['DISTSIGMA']
    M1_sort_prob=np.sort(M1_probdensity/sum(M1_probdensity))[::-1]
    
    #get number of pixels for the 90& probability region
    #M1_count_90=0
    #M1_count_50=0
    M1_count=0
    M1_total=0
    area_value=area_percent/100
    for i in range(0,len(M1_sort_prob)): 
        M1_total+=M1_sort_prob[i]
        if M1_total<=area_value:
            M1_count+=1
        #if M1_total <=0.5:
        #    M1_count_50+=1
        #    M1_count_90+=1
        #elif M1_total<=0.9 and M1_total>0.5: 
        #    M1_count_90+=1
        else: 
            break
    print("superevent_fits: ", superevent_fits)
    #print("appending M1_90prob")
    #M1_90prob.append(M1_count_90*onepixelarea)
    #M1_50prob.append(M1_count_50*onepixelarea)
    M1_prob.append(M1_count*onepixelarea)
    ind.append(event_array[0])
