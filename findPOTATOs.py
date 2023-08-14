#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from datetime import datetime
from astropy.time import Time
from sklearn.neighbors import BallTree
from os.path import exists
from astropy.coordinates import SkyCoord, Angle, Distance
from linking_library import *

#  C.R. Nugent and N. Tan
#  August 2023

########## PARAMETERS ##########
input_filename='image_triplets_20011120.csv'
input_directory='../NEAT_reprocessing/output/'
max_speed = 0.05 #maximum speed an asteroid can travel to be detected, in arcseconds/second
#you don't want this more than ~1/5th of the size of the frame, anything
#faster is both unlikely and undetectable as it will leave the frame before 
#three detections can be made
# angle between a-b-c is variable min_tracklet_angle (degrees)
min_tracklet_angle= 135 #degrees
timing_uncertainty= 5 #seconds
max_mag_variance= 2 #the maximum amount brightness can vary across a tracklet, in mag
# will pick the biggest of these to determine radius of which to search
Maximum_residual = 0.3 #arcseconds #This is the maximum residual allowed after orbfit fit
astrometric_accuracy=3 #arcseconds
###############################

get_night_id=input_filename.split('.')
get_night_id=get_night_id[0].split('_')
night=get_night_id[2]
image_triplets_list=pd.read_csv(input_filename)

for m in np.arange(len(image_triplets_list)):
    file_a='o_sources'+image_triplets_list.filea[m]+'.csv'
    file_b='o_sources'+image_triplets_list.fileb[m]+'.csv'
    file_c='o_sources'+image_triplets_list.filec[m]+'.csv'
    #night='2002010102'
    #file_a='o_sources20020101020555c.csv'
    #file_b='o_sources20020101023551c.csv'
    #file_c='o_sources20020101030745c.csv'

    #file_a='o_sources20020104101451c.csv'
    #file_b='o_sources20020104093953c.csv'
    #file_c='o_sources20020104090911c.csv'

    a=pd.read_csv(input_directory+file_a)
    b=pd.read_csv(input_directory+file_b)
    c=pd.read_csv(input_directory+file_c)
    a_time=Time(a.mjd[0].astype(float), format='mjd', scale='utc')
    b_time=Time(b.mjd[0].astype(float), format='mjd', scale='utc')
    c_time=Time(c.mjd[0].astype(float), format='mjd', scale='utc')

    # # Remove Stationary Sources
    # The nearest neighbors code (balltree, haversine metric)
    # needs RA and Dec in radians. This code is used a lot here
    # so we'll just create the needed columns once right now.
    a['ra_rad'] = np.radians(a['RA'])
    a['dec_rad'] = np.radians(a['Dec'])
    b['ra_rad'] = np.radians(b['RA'])
    b['dec_rad'] = np.radians(b['Dec'])
    c['ra_rad'] = np.radians(c['RA'])
    c['dec_rad'] = np.radians(c['Dec'])
    a_moving, b_moving, c_moving=remove_stationary_sources(a, b, c, astrometric_accuracy) #threshold of 2 arsec is rougly our accuracy near the edges

    #In the following, distances stored in dataframes are in radians, because that's
    #what python likes. RA and DEC are stored in degrees.
    start_time=datetime.now()
    # After removing stationary sources, tree of b comparing to a.

    # minimum speed an asteroid can travel to be detected, in arcseconds/second
    # this is calculated based on our astromemtric accuracy and time between frames
    time_interval_s=(b_time-a_time).sec
    print('time interval',time_interval_s)

    max_dist_rad=np.radians(max_speed*time_interval_s/3600)
    min_dist_rad=np.radians(astrometric_accuracy/3600)
    print('max_dist_rad, in degrees', np.degrees(max_dist_rad))
    #print ('TEMPORARY FOR TROUBLESHOOTING RESETTING MAX DIST')
    #max_dist_rad=np.radians(2)

    tree_a = BallTree(a_moving[['ra_rad', 'dec_rad']],leaf_size=5, metric='haversine')
    indicies_b, distances_b = tree_a.query_radius(b_moving[['ra_rad', 'dec_rad']], r=max_dist_rad, return_distance=True)
    #print(indicies_b,distances_b)
    pair_id=[]
    a_source_index=[]
    b_source_index=[]
    ab_dist=[]
    dec_a=[]
    dec_b=[]
    ra_a=[]
    ra_b=[]
    mag_a=[]
    mag_b=[]
    observatory_code=[]
    band=[]
    tracklet_count=0
    #We want every source in a that is within a certain radius of each source in b, but not too slow. 
    for i in range(len(indicies_b)):
        for j in range(len(indicies_b[i])):
            #print(distances_b[i][j],min_dist_rad)
            if(distances_b[i][j] > min_dist_rad):
                pair_id.append(tracklet_count)
                a_source_index.append(indicies_b[i][j])
                b_source_index.append(i)
                ab_dist.append(distances_b[i][j])
                dec_a.append(a_moving['Dec'][indicies_b[i][j]])
                dec_b.append(b_moving['Dec'][i])
                ra_a.append(a_moving['RA'][indicies_b[i][j]])
                ra_b.append(b_moving['RA'][i])
                mag_a.append(a_moving['magnitude'][indicies_b[i][j]])
                mag_b.append(b_moving['magnitude'][i])
                observatory_code.append(b_moving['observatory_code'][i])
                band.append(b_moving['band'][i])
                tracklet_count += 1
            else:
                print("Warning: Distance less than min_dist_rad. Stationary source removal not working properly.")        

    #to pandas df
    candidate_tracklet=pd.DataFrame(pair_id)
    candidate_tracklet['point_a']=a_source_index
    candidate_tracklet['point_b']=b_source_index
    candidate_tracklet['ab_dist']=ab_dist #don't forget this is in radians
    candidate_tracklet['dec_a']=dec_a
    candidate_tracklet['dec_b']=dec_b
    candidate_tracklet['ra_a']=ra_a
    candidate_tracklet['ra_b']=ra_b
    candidate_tracklet['mag_a']=mag_a
    candidate_tracklet['mag_b']=mag_b
    candidate_tracklet['observatory_code']=observatory_code
    candidate_tracklet['band']=band

    time_interval2_s=(c_time-a_time).sec
    #print('time_interval2_s',time_interval2_s)

    # get predicted RA and Dec of where point c should be
    pred_dist=(candidate_tracklet['ab_dist']*time_interval2_s)/time_interval_s #distance likely to travel between frame b and c
    #print('pred',pred_dist,candidate_tracklet['ab_dist'],time_interval_s)
    r_due_to_timing=timing_uncertainty*(candidate_tracklet['ab_dist']/time_interval_s) #radians
    r_due_to_angle=np.arctan(np.radians(180-min_tracklet_angle)*pred_dist)

    #print('r values',r_due_to_timing,r_due_to_angle)
    r_to_search_c=pd.concat([r_due_to_timing, r_due_to_angle], axis=1).max(axis=1)

    #r_to_search_c=np.max(r_due_to_timing.to_numpy,r_due_to_angle.to_numpy)
    #print(r_to_search_c)
    candidate_tracklet['pred_dist']=pred_dist #radians
    candidate_tracklet['r_due_to_angle']=r_due_to_angle #radians

    if len(candidate_tracklet) == 0:
        print("No pairs found.")
        break

    # this next part is the distance between point a, and second point
    # that is projected at ra_a, dec_b
    # not precisely correct, but close enough
    delta_dec=np.radians(candidate_tracklet['dec_b']-candidate_tracklet['dec_a'])
    candidate_tracklet['delta_dec']=delta_dec


    #predict where object will be in frame c
    # you'll be searching around this point
    delta_ra=np.radians(candidate_tracklet.ra_b)-np.radians(candidate_tracklet.ra_a)
    delta_dec=np.radians(candidate_tracklet.dec_b)-np.radians(candidate_tracklet.dec_a)
    pred_c_pos_ra=np.radians(candidate_tracklet.ra_b)+delta_ra
    pred_c_pos_dec=np.radians(candidate_tracklet.dec_b)+delta_dec
    #print(len(pred_c_pos_ra))


    # use yet another ball tree to see if there's any  
    # detections in expected radius from predicted position
    # in frame c
    pred_c = pd.DataFrame(pred_c_pos_ra, columns=['ra_rad'])
    pred_c['dec_rad']=pred_c_pos_dec

    # see if anything is around the predicted position in c
    tree_c = BallTree(c_moving[['ra_rad', 'dec_rad']], leaf_size=5,metric='haversine')
    indicies_c, distances_c = tree_c.query_radius(pred_c[['ra_rad', 'dec_rad']], r=r_to_search_c, return_distance=True)

    # if detection(s) are around precdicted position in c, then add to tracklet list.
    point_c=[]
    dec_c=[]
    ra_c=[]
    mag_c=[]
    bc_dist=[]

    # I keep thinking there's a cleaner way to do this, but 
    # this works so who cares.
    # We're moving the existing tracklet information into
    # a numpy matrix so we can easily add rows and ignore others 
    np_tracklets=candidate_tracklet.to_numpy()
    new_tracklets =[]

    for i in range(len(indicies_c)):
        if indicies_c.size > 0:
            for j in range(len(indicies_c[i])):
                new_tracklets.append(np_tracklets[:][i])
                point_c.append(indicies_c[i][j])
                dec_c.append(c_moving['Dec'][indicies_c[i][j]])
                ra_c.append(c_moving['RA'][indicies_c[i][j]])
                mag_c.append(c_moving['magnitude'][indicies_c[i][j]])
                bc_dist.append(distances_c[i][j])
        else:
            "No length 3 tracklets found in these frames."

    # Reassemble that dataframe 
    complete_tracklets=df = pd.DataFrame(new_tracklets, columns = ['index','point_a','point_b','ab_dist','dec_a','dec_b','ra_a','ra_b','mag_a','mag_b','observatory_code','band','pred_dist','r_due_to_angle','delta_dec'])
    complete_tracklets['point_c']=point_c
    complete_tracklets['dec_c']=dec_c
    complete_tracklets['ra_c']=ra_c
    complete_tracklets['mag_c']=mag_c
    complete_tracklets['bc_dist']=bc_dist 

    if len(complete_tracklets) == 0:
        print("No tracklets found.")
        break

    # # Tracklet screening
    # A slow moving tracklet has a relatively large search radious for
    # point c, meaning that in some cases the resulting tracklet might
    # have an extreme angle between a-b-c (a c is found, but behind a)
    # so do another screening.
    # this assumes an arbitrary distance to calculate angle. 
    # also screens for magnitude
    mag_min=np.min([complete_tracklets.mag_a[1],complete_tracklets.mag_b[1],complete_tracklets.mag_c[1]])
    mag_max=np.max([complete_tracklets.mag_a[1],complete_tracklets.mag_b[1],complete_tracklets.mag_c[1]])

    for i in range(len(complete_tracklets)):
        mag_min=np.min([complete_tracklets.mag_a[i],complete_tracklets.mag_b[i],complete_tracklets.mag_c[i]])
        mag_max=np.max([complete_tracklets.mag_a[i],complete_tracklets.mag_b[i],complete_tracklets.mag_c[i]])
        if (mag_max-mag_min) <max_mag_variance:    
            coordA = SkyCoord(ra=complete_tracklets.ra_a[i],dec= complete_tracklets.dec_a[i],unit=(u.deg, u.deg), distance=70*u.kpc)
            coordB = SkyCoord(ra=complete_tracklets.ra_b[i],dec= complete_tracklets.dec_b[i],unit=(u.deg, u.deg), distance=70*u.kpc)
            coordC = SkyCoord(ra=complete_tracklets.ra_c[i],dec= complete_tracklets.dec_c[i],unit=(u.deg, u.deg), distance=70*u.kpc)
            lenAB = coordA.separation_3d(coordB)
            lenBC = coordB.separation_3d(coordC)
            lenCA = coordC.separation_3d(coordA)
            cosine_angle = ((lenAB ** 2) + (lenBC ** 2) - (lenCA ** 2)) / (2 * lenAB * lenBC)
            angle = np.degrees(np.arccos(cosine_angle))

        # print("Angle between points A, B, and C:", angle, "degrees")
            if angle.value < min_tracklet_angle:
                complete_tracklets.drop(index=[i],inplace=True)
                #print("dropped")
        else:
            complete_tracklets.drop(index=[i],inplace=True)
            
    complete_tracklets.reset_index(inplace=True)


    # now filter based on findorb
    trackletfilename="tracklets_"+night+'.txt'
    for i in range(len(complete_tracklets)):
        findOrbTxt = open(os.path.expanduser("~/.find_orb/fo.txt"),"w")

        tracklet_id='cn0000'+str(i)
        decimal_time_a=str(a_time).split('.')
        decimal_time_b=str(b_time).split('.')
        decimal_time_c=str(c_time).split('.')

        coordA = SkyCoord(ra=complete_tracklets.ra_a[i],dec= complete_tracklets.dec_a[i],unit=(u.deg, u.deg), distance=70*u.kpc)
        coordB = SkyCoord(ra=complete_tracklets.ra_b[i],dec= complete_tracklets.dec_b[i],unit=(u.deg, u.deg), distance=70*u.kpc)
        coordC = SkyCoord(ra=complete_tracklets.ra_c[i],dec= complete_tracklets.dec_c[i],unit=(u.deg, u.deg), distance=70*u.kpc)

        formatted_data = "     "
        formatted_data += "{}".format(tracklet_id)+'  C'
        formatted_data += "{}".format(a_time.strftime('%Y %m %d'))+'.'
        formatted_data += "{:1}".format(decimal_time_a[1][:5])+' '
        formatted_data += coordA.to_string(style='hmsdms',pad=True,sep=' ',precision=2)+'         '
        formatted_data += "{:.1f}".format(complete_tracklets.mag_a[i])+'   '
        formatted_data += complete_tracklets.band[i]+'    '+str(complete_tracklets.observatory_code[i])+'\n'
        
        formatted_data += "     "
        formatted_data += "{}".format(tracklet_id)+'  C'
        formatted_data += "{}".format(b_time.strftime('%Y %m %d'))+'.'
        formatted_data += "{:1}".format(decimal_time_b[1][:5])+' '
        formatted_data += coordB.to_string(style='hmsdms',pad=True,sep=' ',precision=2)+'         '
        formatted_data += "{:.1f}".format(complete_tracklets.mag_b[i])+'   '
        formatted_data += complete_tracklets.band[i]+'    '+str(complete_tracklets.observatory_code[i])+'\n'
        
        formatted_data += "     "
        formatted_data += "{}".format(tracklet_id)+'  C'
        formatted_data += "{}".format(c_time.strftime('%Y %m %d'))+'.'
        formatted_data += "{:1}".format(decimal_time_c[1][:5])+' '
        formatted_data += coordC.to_string(style='hmsdms',pad=True,sep=' ',precision=2)+'         '
        formatted_data += "{:.1f}".format(complete_tracklets.mag_c[i])+'   '
        formatted_data += complete_tracklets.band[i]+'    '+str(complete_tracklets.observatory_code[i])+'\n'

        #print(formatted_data)
        findOrbTxt.writelines(formatted_data)

        findOrbTxt.close()
        trackletFound = find_orb(Maximum_residual, nullResid = True, MOIDLim = True)
        if trackletFound == True:
            if exists(trackletfilename):
                with open(trackletfilename, 'a', encoding="utf-8") as f:
                    f.write(formatted_data)
                    f.close
            else:
                with open(trackletfilename, 'x', encoding="utf-8") as f:
                    f.write(formatted_data)
                    f.close
        else: #drop it
            complete_tracklets.drop(index=[i],inplace=True)

    #save stats
    now = datetime.now() 
    yearmonthday = now.strftime("%Y%m%d")
    outputname="o_linking_"+yearmonthday+".csv"
    run_time=str(now-start_time)
    num_sources=str(len(a)+len(b)+len(c))
    num_tracklets_prescreen=str(len(complete_tracklets))
                
    if exists(outputname):
        with open(outputname, 'a', encoding="utf-8") as f:
            f.write(file_a +","+file_b+","+file_c+","+datetime.today().strftime('%Y-%m-%d')+','+                datetime.today().strftime('%H:%M:%S') +","+ run_time +","+ num_sources +","+ num_tracklets_prescreen+'\n')
        f.close
    else:
        with open(outputname, 'x', encoding="utf-8") as f:
            f.write("filea,fileb,filec,date_corrected,time_corrected,run_time_s,num_sources,num_tracklets_prescreen\n")
            f.write(file_a +","+file_b+","+file_c+","+datetime.today().strftime('%Y-%m-%d')+','+                datetime.today().strftime('%H:%M:%S') +","+ run_time +","+ num_sources +","+ num_tracklets_prescreen+'\n')
        f.close
