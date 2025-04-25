
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from hpmoc import PartialUniqSkymap
from hpmoc.plot import get_wcs, plot, gridplot
import sys
from DetectorMinimaFolder.DetectorMinima import *
from astropy.coordinates import angular_separation
from astropy.coordinates import SkyCoord
#from PercentAreaFolder.PercentArea import *
from O4a_fits.O4aSuperevents import *
#from CoordsProbabilityFolder.CoordsProbabilityModule import *
def single_pixel_to_probability(skymap_coord):
    #print("printing skymap coord")
    #print(skymap_coord)
    #print("Done printing skymap coord")
    if type(skymap_coord) == str:
        skymap_coord_str=skymap_coord
    else:
        skymap_coord_str=str(skymap_coord)
    #print(skymap_coord_str)
    number_str=""
    exponent_str=""
    isExponent=False
    numbers=["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
    str_portion=0
    len_string=len(skymap_coord_str)
    for i in range(len_string):
        #print("i:", i, "str_portion:", str_portion, "str_val", skymap_coord_str[i])
        if str_portion==4:
            exponent_str=exponent_str+skymap_coord_str[i]
        if str_portion==3:
            if skymap_coord_str[i]=="e":
                str_portion=4
                isExponent=True
            else:
                #print("adding")
                number_str=number_str+skymap_coord_str[i]
                #print(number_str)
        if str_portion==2:
            if skymap_coord_str[i]==" ":
                str_portion=3
        if str_portion==1:
            if skymap_coord_str[i] in numbers:
                str_portion=2
                #print("not adding")
                #number_str=number_str+skymap_coord_str[i]
                #print(number_str)
        if str_portion==0:
            if skymap_coord_str[i]=="-":
                str_portion=1
    #print(number_str)
    number_float=float(number_str)
    if isExponent==True:
        exponent_float=float(exponent_str)
        result=number_float*10**exponent_float
    else:
        result=number_float
    return result

def rank_points(skymap):
    ra, dec = skymap.coords()
    n_points = len(ra)
    points_listed=np.zeros((n_points,4))
    for i in range(n_points):
        points_listed[i][1]=i
        skymap_coord=skymap[i]
        points_listed[i][0]=single_pixel_to_probability(skymap_coord)
        points_listed[i][2]=dec[i].value
        points_listed[i][3]=ra[i].value
        #points_listed[i][4]=str(skymap_coord)
    indices_sorted=np.lexsort((points_listed[:,1], points_listed[:,0]))
    points_sorted = points_listed[indices_sorted]
    return points_sorted
    #points_sorted_0=np.array(points_listed, dtype=dtype_test)
    #return np.sort(points_sorted_0, order=['index'])

def max_point(skymap):
    points_sorted=rank_points(skymap)
    return points_sorted[-1]

def max_portion(skymap,portion):
    points_sorted=rank_points(skymap)
    len_ps=len(points_sorted)
    n=int(portion*len_ps)
    return points_sorted[-n:]

def angular_separation_deg(point_ranked, point_test):
    dec_rad=np.deg2rad(point_test[0])
    ra_rad=np.deg2rad(point_test[1])
    point_dec_rad=np.deg2rad(point_ranked[2])
    point_ra_rad=np.deg2rad(point_ranked[3])
    sep_rad=angular_separation(dec_rad, ra_rad, point_dec_rad, point_ra_rad)
    return np.rad2deg(sep_rad)

def angular_separation_deg_area(points_sorted_l, point_test):
    n_l = len(points_sorted_l)
    angular_separations_list=np.zeroes((n_l, 2))
    for i in range(n_l):
        angular_separations_list[i][0]=angular_separation_deg(points_sorted_l[i], point_test)
        angular_separations_list[i][1]=points_sorted_l[i][1]
    return angular_separation_deg_area

def angular_separation_deg_sorted(points_sorted_l, point_test):
    n_l = len(points_sorted_l)
    angular_separations_list=np.zeroes((n_l, 2))
    for i in range(n_l):
        angular_separations_list[i][1]=angular_separation_deg(points_sorted_l[i], point_test)
        angular_separations_list[i][0]=points_sorted_l[i][1] 
    indices_sorted=np.lexsort((angular_separations_list[:,1], angular_separations_list[:,0]))
    angular_separations_sorted = angular_separations_list[indices_sorted]
    return angular_separations_sorted

def angular_separation_min_in_list(points_sorted_l, point_test):
    n_l = len(points_sorted_l)
    angular_separations_list=np.zeroes((n_l, 2))
    for i in range(n_l):
        angular_separations_list[i][1]=angular_separation_deg(points_sorted_l[i], point_test)
        angular_separations_list[i][0]=points_sorted_l[i][1] 
    indices_sorted=np.lexsort((angular_separations_list[:,1], angular_separations_list[:,0]))
    angular_separations_sorted = angular_separations_list[indices_sorted]
    return angular_separations_sorted[0][0]
    
    

file_location="./testing/"
superevent_fits=file_location+'S230627c'+"_fits,"+str(1)+".fits"
superevent_skymap=PartialUniqSkymap.read(superevent_fits, strategy="ligo")
print(rank_points(superevent_skymap))
maxpoint=max_point(superevent_skymap)
print(maxpoint)
print(angular_separation_deg(maxpoint, [0,0]))

