import numpy as np
import math
import sys
from hpmoc import PartialUniqSkymap
from hpmoc.plot import get_wcs, plot, gridplot
from pathlib import Path
from matplotlib import pyplot as plt
from astropy.table import Table
from hpmoc.points import PointsTuple, Rgba
#import astropy.units
#from astropy.coordinates import AltAz, EarthLocation, SkyCoord
#from astropy.time import Time

def CartesianToLatLong(coords):
    '''
    Takes cartesian coordinates and converts them to coordinates in \ndegrees in [latitude, longitude] format.
    '''
    theta=np.arctan2(coords[1], coords[0])
    phi=np.arccos(coords[2]/np.linalg.norm(coords))
    longitude=theta
    latitude=(np.pi/2)-phi
    direction=[np.rad2deg(latitude), np.rad2deg(longitude)]
    return direction
    

def getMinimaRotatingFrame (latitude, longitude, arm_orient_1, arm_orient_2, debug_vectors=False, debug_results=False, debug_dots=False): #arm_orient measured in clockwise angle from the North Pole
    '''
    Finds the location of detector minima in an Earth-centred rotating reference frame.\nLatitude and longitude should be in degrees, \narm_orient_1 and arm_orient_2 should be in degrees clockwise from north. \ndebug_vectors shows the Cartesian vectors, debug_results shows latitude and longitude, \nand debug_dots shows the dot products of the Cartesian vectors.
    '''
    arm_average = (arm_orient_1+arm_orient_2)/2
    #convert to radians
    r_latitude = np.deg2rad(latitude)
    r_longitude = np.deg2rad(longitude)
    arm_average_r=np.deg2rad(arm_average)
    #create basis
    zenith_vector = [np.cos(r_longitude)*np.cos(r_latitude), np.sin(r_longitude)*np.cos(r_latitude), np.sin(r_latitude)]
    if longitude>0:
        n_longitude = r_longitude+np.pi
        n_latitude = (np.pi/2)-r_latitude
    else:
        n_longitude = r_longitude
        n_latitude = (np.pi/2)+r_latitude
    north_vector = [np.cos(n_longitude)*np.cos(n_latitude), np.sin(n_longitude)*np.cos(n_latitude), np.sin(n_latitude)]
    side_vector=np.cross(zenith_vector, north_vector)
    #find arm average components
    print(arm_average)
    
    aa_component_n=np.cos(arm_average_r)
    aa_component_side=-np.sin(arm_average_r)
    
    p1_component_side=-np.cos(arm_average_r)
    p1_component_n=-np.sin(arm_average_r)
    
    p3_component_side=np.cos(arm_average_r)
    p3_component_n=np.sin(arm_average_r)
    
    a_aa_component_n=-np.cos(arm_average_r)
    a_aa_component_side=np.sin(arm_average_r)
    aa_vector=np.zeros(3)
    p1_vector=np.zeros(3)
    p3_vector=np.zeros(3)
    a_aa_vector=np.zeros(3)
    for i in [0,1,2]:
        aa_vector[i]=aa_component_n*north_vector[i]+aa_component_side*side_vector[i]
        p1_vector[i]=p1_component_n*north_vector[i]+p1_component_side*side_vector[i]
        p3_vector[i]=p3_component_n*north_vector[i]+p3_component_side*side_vector[i]
        a_aa_vector[i]=a_aa_component_n*north_vector[i]+a_aa_component_side*side_vector[i]
    aa_rotating_direction=CartesianToLatLong(aa_vector)
    p1_rotating_direction=CartesianToLatLong(p1_vector)
    p3_rotating_direction=CartesianToLatLong(p3_vector)
    a_aa_rotating_direction=CartesianToLatLong(a_aa_vector)
    if debug_vectors==True:
        print("vectors for debug")
        print("zenith vector:", zenith_vector)
        print("north vector:", north_vector)
        print("side vector:", side_vector)
        print("arm average vector:", aa_vector)
        print("perpendicular vector (Q1):", p1_vector)
        print("perpendicular vector (Q3):", p3_vector)
        print("opposite to arm average vector:", a_aa_vector)
        print("dot product", np.dot(aa_vector, p1_vector))
    if debug_results==True:
        print("results for debug")
        print(aa_rotating_direction)
        print(p1_rotating_direction)
        print(p3_rotating_direction)
        print(a_aa_rotating_direction)
    if debug_dots==True:
        print("dots for debug")
        print("quadrant vector dots (should be 0)")
        print(np.dot(aa_vector, p1_vector))
        print(np.dot(p1_vector, a_aa_vector))
        print(np.dot(a_aa_vector, p3_vector))
        print(np.dot(p3_vector, aa_vector))
        print("zenith vector dots (should be 0)")
        print(np.dot(zenith_vector, p1_vector))
        print(np.dot(zenith_vector, a_aa_vector))
        print(np.dot(zenith_vector, p3_vector))
        print(np.dot(zenith_vector, aa_vector))
        print("reference frame dots (should be 0)")
        print(np.dot(zenith_vector, north_vector))
        print(np.dot(zenith_vector, side_vector))
        print(np.dot(north_vector, side_vector))
        print("opposites dots (should be 1)")
        print(np.dot(aa_vector, a_aa_vector))
        print(np.dot(p3_vector, p1_vector))
    results_array=np.zeros((4,2))
    results_array[0][0]=aa_rotating_direction[0]
    results_array[0][1]=aa_rotating_direction[1]
    results_array[1][0]=p1_rotating_direction[0]
    results_array[1][1]=p1_rotating_direction[1]
    results_array[2][0]=p3_rotating_direction[0]
    results_array[2][1]=p3_rotating_direction[1]
    results_array[3][0]=a_aa_rotating_direction[0]
    results_array[3][1]=a_aa_rotating_direction[1]
    return results_array

def getCoordsAtJDT (coords_rotating, timeJulian, epoch=2451545.0, zero_at_epoch=280.17, debug=False):
    '''
    Converts coordinates from an Earth-centric rotating frame to an Earth-centric fixed frame. \ncoords_rotating should be in degrees in [latitude, longitude] format. \ntimeJulian should be the Julian date of the event. \nDefault epoch is J2000. \nzero_at_epoch is the right ascension of the projection of the Earth's prime meridian onto the celestial sphere at the epoch.
    '''
    #at J2000, the right ascension of the prime meridian was +280.17 degrees, to within 15 arcseconds
    timeJulianSinceEpoch=timeJulian-epoch
    SiderealDaysJ=timeJulianSinceEpoch*24.065709824279/24
    SiderealRotation=(SiderealDaysJ-(int(SiderealDaysJ)))*360
    SiderealRotationDegrees=SiderealRotation%360
    latitude=coords_rotating[0]
    longitude_rotating=coords_rotating[1]
    #moves from west to east
    zero_at_time=zero_at_epoch+SiderealRotationDegrees
    longitude=(zero_at_time+longitude_rotating)%360
    coords_at_time = [latitude, longitude]
    if debug==True:
        print("timeJulianSinceEpoch:", timeJulianSinceEpoch)
        print("SiderealDaysJ:", SiderealDaysJ)
        print("SiderealRotation:", SiderealRotation)
        print("SiderealRotationDegrees:", SiderealRotationDegrees)
        print("latitude:", latitude)
        print("longitude_rotating:", longitude_rotating)
        print("zero_at_time:", zero_at_time)
        print("longitude:", longitude)
    return coords_at_time

def getMinimaAtJDTArray (detector_stats, timeJulian, epoch=2451545.0, zero_at_epoch=280.17):
    '''
    Creates an array of the coordinates of detector minima at a given time. \ndetector_stats should be in [latitude, longitude, arm_orient_1, arm_orient_2] \n where latitude and longitude are in degrees and arm_orients should be in degrees clockwise from due north.
    '''
    latitude=detector_stats[0]
    longitude=detector_stats[1]
    arm_orient_1=detector_stats[2]
    arm_orient_2=detector_stats[3]
    RF_array=getMinimaRotatingFrame(latitude, longitude, arm_orient_1, arm_orient_2)
    c_array = np.zeros((4,2))
    for i in range(len(RF_array)):
        s_coords=[0,0]
        s_coords=getCoordsAtJDT(RF_array[i], timeJulian, epoch=epoch, zero_at_epoch=zero_at_epoch)
        c_array[i][0]=s_coords[0]
        c_array[i][1]=s_coords[1]
    return c_array

def CoordsTuple (coords_array, dispersion=3.0):
    '''
    Converts an array of four [latitude, longitude] coordinates into an astropy PointsTuple.
    '''
    len_c = len(coords_array)
    pts1=PointsTuple(
        [
            (coords_array[0][1], coords_array[0][0], dispersion),
            (coords_array[1][1], coords_array[1][0], dispersion),
            (coords_array[2][1], coords_array[2][0], dispersion),
            (coords_array[3][1], coords_array[3][0], dispersion),
        ]
    )
    return pts1

def CoordsScatter (coords_array, skymap_m, filename, dispersion=3.0):
    '''
    Places a PointsTuple of an array of four [latitude, longitude] coordinates on a skymap.
    '''
    len_c = len(coords_array)
    skymap_m.plot(CoordsTuple(coords_array, dispersion))
    plt.savefig(filename+".png")

H1_detector_stats = [46.5, -119.5, -36, -126]
L1_detector_stats = [30.5, -90.75, 162, 252]
V1_detector_stats = [43.5, 10.5, 19, -71]
J2000=2451545.0
GW170817_JDT=2457983.028518



#Hanford
#H1RotatingNorm=getMinimaRotatingFrame(46.5, -119.5, -36, -126)
#Livingston
#L1RotatingNorm=getMinimaRotatingFrame(30.5, -90.75, 162, 252)
#Virgo
#V1RotatingNorm=getMinimaRotatingFrame(43.5, 10.5, 19, -71)

#print("170817")
#print(GW170817_JDT)
#for i in range(4):
#    print("H1 minimum")
#    print(getCoordsAtJDT(H1RotatingNorm[i], GW170817_JDT))
#    print("L1 minimum")
#    print(getCoordsAtJDT(L1RotatingNorm[i], GW170817_JDT))
#    print("V1 minimum")
#    print(getCoordsAtJDT(V1RotatingNorm[i], GW170817_JDT))
#
#print(getMinimaAtJDTArray(H1_detector_stats, GW170817_JDT))
#print(getMinimaAtJDTArray(L1_detector_stats, GW170817_JDT))
#print(getMinimaAtJDTArray(V1_detector_stats, GW170817_JDT))
