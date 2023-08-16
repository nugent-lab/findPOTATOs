import astropy.units as u
import os 
from subprocess import Popen # used to call Find_Orb
import subprocess 
import re # regular expressions, used to search for mean residuals in Find_orb output files
from time import sleep
import pandas as pd
import numpy as np
from datetime import datetime
from astropy.time import Time
from sklearn.neighbors import BallTree



def find_orb(maxResidual, nullResid = True, MOIDLim = False):
    """
    Feeds observations in MPC format that are located ~/.find_orb/fo.txt to
    the non-interactive version of find orb, fo. find_orb stores orbital
    elements in  ~/.find_orb/elements.txt, which this function will read to
    find the mean residual to the orbital fit. If the mean residual is less
    than maxResidual (specified in ") and all observations in
    ~/.find_orb/fo.txt was used to generate the orbital fit, then the
    function will return True. In other cases (e.g. find_orb doesn't run;
    mean residual greater than maxResidual; not all observations in
    ~/.find_orb/fo.txt used), the function will return False.
    """
    elements_path="~/.find_orb/elements.txt" #for mac 
    if os.path.exists(os.path.expanduser(elements_path)):
        os.remove(os.path.expanduser(elements_path))
    # this line works on mac & some unix installs but not the MGHPCC
    #sp = Popen(['cd ~/.find_orb\n~/find_orb/find_orb/fo fo.txt -c'], shell=True)
    # this line is for the MGHPCC. Either way, you need the directory where your fo files are
    # the subprocess module reacts poorly to the supercomputer.
    #os.system('fo fo.txt -c')
    sp = subprocess.call(['fo fo.txt -c'], shell=True)
    totSleep = 0
    # wait for find_orb to create elements.txt. If it takes longer than 20 seconds
    # then find_orb probably can't find an orbit.
    while not os.path.exists(os.path.expanduser(elements_path)):
        sleep(0.2)
        totSleep = totSleep + 0.2
        if totSleep > 20:
            break
    if os.path.exists(os.path.expanduser(elements_path)):
        if os.path.getsize(os.path.expanduser(elements_path)) == 0:
            sleep(0.2)
        #numObs = sum(1 for line in open(os.path.expanduser("~/.find_orb/fo.txt")))
        numObs = sum(1 for line in open(os.path.expanduser("/projectnb/ct-ast/findPOTATOs/fo.txt")))


        # save all inputs to find_orb
        #open("outputs/AllPotentialTracklets.txt", "a+").writelines([l for l in open(os.path.expanduser("~/.find_orb/fo.txt")).readlines()])
        open("outputs/AllPotentialTracklets.txt", "a+").writelines([l for l in open(os.path.expanduser("/projectnb/ct-ast/findPOTATOs/fo.txt")).readlines()])
        for line in open(os.path.expanduser(elements_path)):
            li=line.strip()
            if not li.startswith("#"):
                open("outputs/AllPotentialTracklets.txt", "a").writelines(line.rstrip())
                open("outputs/AllPotentialTracklets.txt", "a").writelines("\n")
        open("outputs/AllPotentialTracklets.txt", "a").writelines("\n\n")

        resCheck = False
        for line in open(os.path.expanduser(elements_path)):
            match = re.search('mean residual (\d+)".(\d+)', line)
            match2 = re.search('MOID: (\d+).(\d+)', line)
            if match:
                res = int(match.group(1)) + float(('0.'+match.group(2)))
                if nullResid:
                    if (res < maxResidual) & (res > 0): # maxResidual in "
                        resCheck = True
                    else:
                        resCheck = False
                else:
                    if (res < maxResidual): # maxResidual in "
                        resCheck = True
                    else:
                        resCheck = False
            if (match2):
                if MOIDLim:
                    MOID = int(match2.group(1)) + float(('0.'+match2.group(2)))
                    if MOID > MOIDLim:
                        print("MOID:",MOID," exceeds MOIDLim:", MOIDLIM)
                        break
        if  resCheck:
            return True
        else:
            print("Residuals,",res," exceed maxResidual:", maxResidual)
            return False
    else:
        return False


def remove_stationary_sources(df1, df2, df3, thresh):
    """
    Compares three dataframes, removes sources that are at same coordinate 
    location to within threshold. Returns cleaned dataframes that 
    consist of transitory sources.
    
    The for loop through the matches is a bit unpythonic, could be
    improved but works for now.
    
    Args:
        df1: first dataset to be considered, dataframe. Needs columns ra_rad, dec_rad (ra and dec in radians).
        df2: second dataset that will be compared to first.  Needs columns ra_rad, dec_rad (ra and dec in radians).
        df3: third datset that will be compared to first.  Needs columns ra_rad, dec_rad (ra and dec in radians).
        thresh: threshold that we should consider stationary sources, in arcseconds. 

    Returns:
        df1_moving: just the transitory sources in the first dataframe
        df2_moving: just the transitory sources in the second dataframe
        df3_moving: just the transitory sources in the third dataframe
    """
    #convert threshold (arsec) to degrees, then to radians
    thresh_rad=np.radians(thresh/3600)
    #print("thresh_rad",thresh_rad)

    #intialize output dataframes
    #gonna delete the duplicates from the df?_moving dataframes
    df1_moving=df1.copy()
    df2_moving=df2.copy()
    df3_moving=df3.copy()
    
    #compare 1 to 2
    tree1 = BallTree(df1[['ra_rad', 'dec_rad']],leaf_size=5, metric='haversine')
    indicies = tree1.query_radius(df2[['ra_rad', 'dec_rad']], r=thresh_rad)
    for i in range(len(indicies)):
        for j in range(len(indicies[i])):
            try: #using a try, since there's gonna be duplicates between comparing a to b, b to c, etc
                df1_moving.drop(index=indicies[i][j], inplace=True)
            except:
                pass 
    
    #compare 1 to 3
    indicies = tree1.query_radius(df3[['ra_rad', 'dec_rad']], r=thresh_rad)
    for i in range(len(indicies)):
        for j in range(len(indicies[i])):
            try:
                df1_moving.drop(index=indicies[i][j], inplace=True)
            except:
                pass 
    
    #compare 2 to 1
    tree2 = BallTree(df2[['ra_rad', 'dec_rad']],leaf_size=5,  metric='haversine')
    indicies = tree2.query_radius(df1[['ra_rad', 'dec_rad']], r=thresh_rad)
    for i in range(len(indicies)):
        for j in range(len(indicies[i])):
            try:
                df2_moving.drop(index=indicies[i][j], inplace=True)
            except:
                pass
        
    #compare 2 to 3
    indicies = tree2.query_radius(df3[['ra_rad', 'dec_rad']], r=thresh_rad)
    for i in range(len(indicies)):
        for j in range(len(indicies[i])):
            try:
                df2_moving.drop(index=indicies[i][j], inplace=True)
            except:
                pass
    
    #compare 3 to 1
    tree3 = BallTree(df3[['ra_rad', 'dec_rad']],leaf_size=5,  metric='haversine')
    indicies = tree3.query_radius(df1[['ra_rad', 'dec_rad']], r=thresh_rad)
    for i in range(len(indicies)):
        for j in range(len(indicies[i])):
            try:
                df3_moving.drop(index=indicies[i][j], inplace=True)
            except:
                pass
        
    #compare 3 to 2
    indicies = tree3.query_radius(df2[['ra_rad', 'dec_rad']], r=thresh_rad)
    for i in range(len(indicies)):
        for j in range(len(indicies[i])):
            try:
                df3_moving.drop(index=indicies[i][j], inplace=True)
            except:
                pass
    
    # Reset index after dropping rows
    df1_moving.reset_index(inplace=True)
    df2_moving.reset_index(inplace=True)
    df3_moving.reset_index(inplace=True)
    
    print("Percentage of sources remaining from original dataframes df1, df2, df3:", 
        round(len(df1_moving) / len(df1), 2),
        round(len(df2_moving) / len(df2), 2),
        round(len(df3_moving) / len(df3), 2))
    return df1_moving, df2_moving, df3_moving
