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
from ades import *
import sys
import fitsio  # if you want thumbnails of sources

#  C.R. Nugent and N. Tan
#  August 2023

########## PARAMETERS ##########
# input filename provided at call
input_directory = "../NEAT_reprocessing/output/"

print_thumbs = (
    "y"  # turn this on (='y') if you want to save thumbnails of the detections
)
# in all tracklets.
fits_path_1 = "../tricam/data/p" # location of fits files
fits_path_2="/obsdata/processed/"  # the night is sandwitched between these two

export_ades = (
    "y"  # turn this on (='y') if you want observations exported in XML ADES format
)
# more on ADES here: https://minorplanetcenter.net/iau/info/ADES.html

# input_directory='sources/'
max_speed = (
    0.05  # maximum speed an asteroid can travel to be detected, in arcseconds/second
)
# you don't want this more than ~1/5th of the size of the frame, anything
# faster is both unlikely and undetectable as it will leave the frame before
# three detections can be made
# angle between a-b-c is variable min_tracklet_angle (degrees)
min_tracklet_angle = 135  # degrees
timing_uncertainty = 5  # seconds
max_mag_variance = 3  # the maximum amount brightness can vary across a tracklet, in mag
# will pick the biggest of these to determine radius of which to search
Maximum_residual = (
    0.9  # arcseconds #This is the maximum residual allowed after orbfit fit
)
astrometric_accuracy = 1.5  # arcseconds
findorb_check = (
    "n"  # if =='y', check tracklets using Bill Gray's Find Orb for accuracy.
)
exposure_correction = (
    10  # seconds. This code takes input as time at beginning of exposure.
)
# The MPC wants the time of the midpoint of exposure. Exposure times are 20 seconds, so
# this code requires a exposure correction of +10 seconds.

########## ADES PARAMETERS ##########
# header information. None of these will be changed by the
# following code.
ades_dict = {
    "mpcCode": "644",  # MPC-assigned observatory code
    "observatoryName": "Palomar Mountain/NEAT",
    "submitter": "C. Nugent",
    "observers1": "K. Lawrence",
    "observers2": "E. Helin",
    "measurers": "C. Nugent",
    "coinvestigators1": "J. (G.) Bauer",
    "coinvestigators2": "Y. Kim",
    "telescope_design": "reflector",
    "telescope_aperture": "1.1",
    "telescope_detector": "CCD",
    "fundingSource": "NASA",
    "comment": "None",
}
# observation information. Some of these are dummy values that will
# be updated later.
ades_obs_dict = {
    # various codes can be found here:
    # https://www.minorplanetcenter.net/iau/info/ADESFieldValues.html
    #'permID': '04933',#IAU permanent designation
    #'provID': '2022 LB1',#MPC provisional designation (in unpacked form)
    # for unnumbered objects.
    "trkSub": "None",  # Observer-assigned tracklet identifier
    "mode": "CCD",  # Mode of instrumentation (probably CCD)
    "stn": "644",  # Observatory code assigned by the MPC
    # UTC date and time of the observation, ISO 8601 exended format,
    # i.e. yyyy-mm-ddThh:mm:ss.sssZ.
    # The reported time precision should be appropriate for the
    # astrometric accuracy, but no more than 6 digits are permitted
    # after the decimal. The trailing Z indicates UTC and is required.
    "obsTime": "2016-08-29T12:23:34.12Z",
    #'rmsTime': '3' #Random uncertainty in obsTime in seconds as estimated by the observer
    "ra": "215.10254",  # decimal degrees in the J2000.0 reference frame
    "dec": "-12.547",  # decimal degrees in the J2000.0 reference frame
    # For ra-dec and deltaRA- deltaDec observations, the random component
    # of the RA*COS(DEC) and DEC uncertainty (1σ) in arcsec as estimated
    # by the observer as part of the image processing and astrometric reduction.
    "rmsRA": "1.5",
    "rmsDec": "1.5",
    # Correlation between RA and DEC or between distance and PA that may
    # result from the astrometric reduction.
    #'rmsCorr': '-0.214',
    "astCat": "UBSC",  # Star catalog used for the astrometric reduction
    # ‘UNK’, will be used for some archival observations to indicate that
    # the astrometric catalog is unknown.
    "mag": "21.91",  # Apparent magnitude in specified band.
    "rmsMag": "0.15",  # Apparent magnitude uncertainty, 1-sigma
    "band": "g",  # Passband designation for photometry.
    "photCat": "Gaia3",  # Star catalog used for the photometric reduction.
    # full list here: https://www.minorplanetcenter.net/iau/info/ADESFieldValues.html
    #'photAp': '13.3', #Photometric aperture radius in arcsec.
    #'logSNR': '0.78', #The log10 of the signal-to-noise ratio of the source
    # in the image integrated on the entire aperture used for the astrometric centroid.
    #'seeing': '0.8', #Size of seeing disc in arcsec, measured at Full-Width,
    # Half-Max (FWHM) of target point spread function (PSF).
    "exp": "20.0",  # Exposure time in seconds.
    #'remarks': 'None' #A comment provided by the observer. This field can be
    # used to report additional information that is not reportable in the notes
    # field, but that may be of relevance for interpretation of the observations.
    # Should be used sparingly by major producers.
}
#################################


input_filename = sys.argv[1]

get_night_id = input_filename.split(".")
get_night_id = get_night_id[0].split("_")
night = get_night_id[2]
fits_path=fits_path_1+night+fits_path_2
image_triplets_list = pd.read_csv(input_filename)
tracklet_num = 0



if export_ades == "y":
    xml_filename = "outputs/tracklets_ADES_" + night + ".xml"
    xml_tracklet_found = "n"


for m in np.arange(len(image_triplets_list)):
    file_a = "o_sources" + image_triplets_list.filea[m] + ".csv"
    file_b = "o_sources" + image_triplets_list.fileb[m] + ".csv"
    file_c = "o_sources" + image_triplets_list.filec[m] + ".csv"
    print("Checking file triplet number", m, "consisting of:", file_a, file_b, file_c)

    # Put frames in exposure order, so that frame a is first, b is second, and c is third.
    init_a = pd.read_csv(input_directory + file_a)
    init_b = pd.read_csv(input_directory + file_b)
    init_c = pd.read_csv(input_directory + file_c)

    init_a_time = Time(init_a.mjd[0].astype(float), format="mjd", scale="utc")
    init_b_time = Time(init_b.mjd[0].astype(float), format="mjd", scale="utc")
    init_c_time = Time(init_c.mjd[0].astype(float), format="mjd", scale="utc")

    # put frames in order
    order_frames = pd.DataFrame(
        {
            "names": [file_a, file_b, file_c],
            "times": [init_a_time, init_b_time, init_c_time],
        }
    )

    order_frames.sort_values(by=["times"], inplace=True)
    order_frames.reset_index(inplace=True)
    # print(order_frames)

    a = pd.read_csv(input_directory + order_frames.names[0])
    b = pd.read_csv(input_directory + order_frames.names[1])
    c = pd.read_csv(input_directory + order_frames.names[2])
    # correct for time- MPC wants midponint of expsure, but this code takes
    # the input as beginning of exposure time
    exposure_correction_mjd = exposure_correction / (
        24 * 60 * 60
    )  # convert from seconds to fraction of days
    a_time = Time(a.mjd[0].astype(float), format="mjd", scale="utc")
    b_time = Time(b.mjd[0].astype(float), format="mjd", scale="utc")
    c_time = Time(c.mjd[0].astype(float), format="mjd", scale="utc")
    a_time += exposure_correction_mjd
    b_time += exposure_correction_mjd
    c_time += exposure_correction_mjd
    decimal_time_a = str(a_time).split(".")
    decimal_time_b = str(b_time).split(".")
    decimal_time_c = str(c_time).split(".")
    # print("frames,", order_frames.names[0],order_frames.names[1],order_frames.names[2])
    # print(a_time,b_time,c_time)

    # # Remove Stationary Sources
    # The nearest neighbors code (balltree, haversine metric)
    # needs RA and Dec in radians. This code is used a lot here
    # so we'll just create the needed columns once right now.
    a["ra_rad"] = np.radians(a["RA"])
    a["dec_rad"] = np.radians(a["Dec"])
    b["ra_rad"] = np.radians(b["RA"])
    b["dec_rad"] = np.radians(b["Dec"])
    c["ra_rad"] = np.radians(c["RA"])
    c["dec_rad"] = np.radians(c["Dec"])

    if print_thumbs == "y":
        if not os.path.exists("thumbs/"):
            try:
                os.system("mkdir thumbs")
            except:
                print("Could not make thumbnail directory.")

        # derive filenames of original fits files
        # then load images with fitsio
        fits_frame_a = order_frames.names[0].split("s")
        fits_frame_a = fits_frame_a[2]
        fits_frame_a = fits_frame_a[:-1] + "fit.fz"

        fits_frame_b = order_frames.names[1].split("s")
        fits_frame_b = fits_frame_b[2]
        fits_frame_b = fits_frame_b[:-1] + "fit.fz"

        fits_frame_c = order_frames.names[2].split("s")
        fits_frame_c = fits_frame_c[2]
        fits_frame_c = fits_frame_c[:-1] + "fit.fz"
        print("images", fits_frame_a, fits_frame_b, fits_frame_c)

        telescope_image_a = fitsio.read(fits_path + fits_frame_a)
        telescope_image_b = fitsio.read(fits_path + fits_frame_b)
        telescope_image_c = fitsio.read(fits_path + fits_frame_c)



    a_moving, b_moving, c_moving = remove_stationary_sources(
        a, b, c, astrometric_accuracy
    )  # threshold of 2 arsec is rougly our accuracy near the edges

    # In the following, distances stored in dataframes are in radians, because that's
    # what python likes. RA and DEC are stored in degrees.
    start_time = datetime.now()
    # After removing stationary sources, tree of b comparing to a.

    # minimum speed an asteroid can travel to be detected, in arcseconds/second
    # this is calculated based on our astromemtric accuracy and time between frames
    time_interval_s = (b_time - a_time).sec
    # print('time interval',time_interval_s, 'seconds')

    max_dist_rad = np.radians(max_speed * time_interval_s / 3600)
    min_dist_rad = np.radians(astrometric_accuracy / 3600)
    # print('max_dist_rad, in degrees', np.degrees(max_dist_rad))
    # print ('TEMPORARY FOR TROUBLESHOOTING RESETTING MAX DIST')
    # max_dist_rad=np.radians(2)

    tree_a = BallTree(a_moving[["ra_rad", "dec_rad"]], leaf_size=5, metric="haversine")
    indicies_b, distances_b = tree_a.query_radius(
        b_moving[["ra_rad", "dec_rad"]], r=max_dist_rad, return_distance=True
    )
    # print(indicies_b,distances_b)
    pair_id = []
    a_source_index = []
    b_source_index = []
    ab_dist = []
    dec_a = []
    dec_b = []
    ra_a = []
    ra_b = []
    x_a = []
    x_b = []
    y_a = []
    y_b = []
    mag_a = []
    mag_b = []
    observatory_code = []
    band = []
    tracklet_count = 0
    # We want every source in a that is within a certain radius of each source in b, but not too slow.
    for i in range(len(indicies_b)):
        for j in range(len(indicies_b[i])):
            # print(distances_b[i][j],min_dist_rad)
            if distances_b[i][j] > min_dist_rad:
                pair_id.append(tracklet_count)
                a_source_index.append(indicies_b[i][j])
                b_source_index.append(i)
                ab_dist.append(distances_b[i][j])
                dec_a.append(a_moving["Dec"][indicies_b[i][j]])
                dec_b.append(b_moving["Dec"][i])
                ra_a.append(a_moving["RA"][indicies_b[i][j]])
                ra_b.append(b_moving["RA"][i])
                x_a.append(a_moving["xcentroid"][indicies_b[i][j]])
                x_b.append(b_moving["xcentroid"][i])
                y_a.append(a_moving["ycentroid"][indicies_b[i][j]])
                y_b.append(b_moving["ycentroid"][i])
                mag_a.append(a_moving["magnitude"][indicies_b[i][j]])
                mag_b.append(b_moving["magnitude"][i])
                observatory_code.append(b_moving["observatory_code"][i])
                band.append(b_moving["band"][i])
                tracklet_count += 1
            else:
                print(
                    "Warning: Distance less than min_dist_rad. Stationary source removal not working properly."
                )

    # to pandas df
    candidate_tracklet = pd.DataFrame(pair_id)
    candidate_tracklet["point_a"] = a_source_index
    candidate_tracklet["point_b"] = b_source_index
    candidate_tracklet["ab_dist"] = ab_dist  # don't forget this is in radians
    candidate_tracklet["dec_a"] = dec_a
    candidate_tracklet["dec_b"] = dec_b
    candidate_tracklet["ra_a"] = ra_a
    candidate_tracklet["ra_b"] = ra_b
    candidate_tracklet["x_a"] = x_a
    candidate_tracklet["y_a"] = y_a
    candidate_tracklet["x_b"] = x_b
    candidate_tracklet["y_b"] = y_b
    candidate_tracklet["mag_a"] = mag_a
    candidate_tracklet["mag_b"] = mag_b
    candidate_tracklet["observatory_code"] = observatory_code
    candidate_tracklet["band"] = band

    time_interval2_s = (c_time - a_time).sec
    # print('time_interval2_s',time_interval2_s)

    # get predicted RA and Dec of where point c should be
    pred_dist = (
        candidate_tracklet["ab_dist"] * time_interval2_s
    ) / time_interval_s  # distance likely to travel between frame b and c
    # print('pred',pred_dist,candidate_tracklet['ab_dist'],time_interval_s)
    r_due_to_timing = timing_uncertainty * (
        candidate_tracklet["ab_dist"] / time_interval_s
    )  # radians
    r_due_to_angle = np.arctan(np.radians(180 - min_tracklet_angle) * pred_dist)

    # print('r values',r_due_to_timing,r_due_to_angle)
    r_to_search_c = pd.concat([r_due_to_timing, r_due_to_angle], axis=1).max(axis=1)

    # r_to_search_c=np.max(r_due_to_timing.to_numpy,r_due_to_angle.to_numpy)
    # print(r_to_search_c)
    candidate_tracklet["pred_dist"] = pred_dist  # radians
    candidate_tracklet["r_due_to_angle"] = r_due_to_angle  # radians

    if len(candidate_tracklet) == 0:
        print("No pairs found.")
        continue  # skip to next frame triplet
    # this next part is the distance between point a, and second point
    # that is projected at ra_a, dec_b
    # not precisely correct, but close enough
    delta_dec = np.radians(candidate_tracklet["dec_b"] - candidate_tracklet["dec_a"])
    candidate_tracklet["delta_dec"] = delta_dec

    # predict where object will be in frame c
    # you'll be searching around this point
    delta_ra = np.radians(candidate_tracklet.ra_b) - np.radians(candidate_tracklet.ra_a)
    delta_dec = np.radians(candidate_tracklet.dec_b) - np.radians(
        candidate_tracklet.dec_a
    )
    pred_c_pos_ra = np.radians(candidate_tracklet.ra_b) + delta_ra
    pred_c_pos_dec = np.radians(candidate_tracklet.dec_b) + delta_dec

    # use yet another ball tree to see if there's any
    # detections in expected radius from predicted position
    # in frame c
    pred_c = pd.DataFrame(pred_c_pos_ra, columns=["ra_rad"])
    pred_c["dec_rad"] = pred_c_pos_dec

    # see if anything is around the predicted position in c
    tree_c = BallTree(c_moving[["ra_rad", "dec_rad"]], leaf_size=5, metric="haversine")
    indicies_c, distances_c = tree_c.query_radius(
        pred_c[["ra_rad", "dec_rad"]], r=r_to_search_c, return_distance=True
    )

    # if detection(s) are around precdicted position in c, then add to tracklet list.
    point_c = []
    dec_c = []
    ra_c = []
    x_c = []
    y_c = []
    mag_c = []
    bc_dist = []

    # I keep thinking there's a cleaner way to do this, but
    # this works so who cares.
    # We're moving the existing tracklet information into
    # a numpy matrix so we can easily add rows and ignore others
    np_tracklets = candidate_tracklet.to_numpy()
    new_tracklets = []

    for i in range(len(indicies_c)):
        if indicies_c.size > 0:
            for j in range(len(indicies_c[i])):
                new_tracklets.append(np_tracklets[:][i])
                point_c.append(indicies_c[i][j])
                dec_c.append(c_moving["Dec"][indicies_c[i][j]])
                ra_c.append(c_moving["RA"][indicies_c[i][j]])
                x_c.append(c_moving["xcentroid"][indicies_c[i][j]])
                y_c.append(c_moving["ycentroid"][indicies_c[i][j]])
                mag_c.append(c_moving["magnitude"][indicies_c[i][j]])
                bc_dist.append(distances_c[i][j])
        else:
            "No length 3 tracklets found in these frames."
    # Reassemble that dataframe
    complete_tracklets = pd.DataFrame(
        new_tracklets,
        columns=[
            "index",
            "point_a",
            "point_b",
            "ab_dist",
            "dec_a",
            "dec_b",
            "ra_a",
            "ra_b",
            "x_a",
            "y_a",
            "x_b",
            "y_b",
            "mag_a",
            "mag_b",
            "observatory_code",
            "band",
            "pred_dist",
            "r_due_to_angle",
            "delta_dec",
        ],
    )
    complete_tracklets["point_c"] = point_c
    complete_tracklets["dec_c"] = dec_c
    complete_tracklets["ra_c"] = ra_c
    complete_tracklets["x_c"] = x_c
    complete_tracklets["y_c"] = y_c
    complete_tracklets["mag_c"] = mag_c
    complete_tracklets["bc_dist"] = bc_dist

    if len(complete_tracklets) == 0:
        print("No tracklets found.")
        continue  # skip to next frame triplet

    # # Tracklet screening
    # A slow moving tracklet has a relatively large search radious for
    # point c, meaning that in some cases the resulting tracklet might
    # have an extreme angle between a-b-c (a c is found, but behind a)
    # so do another screening.
    # this assumes an arbitrary distance to calculate angle.
    # also screens for magnitude
    angle_array = []
    mag_array = []
    mag_min_array = []
    for i in range(len(complete_tracklets)):
        mag_min = np.min(
            [
                complete_tracklets.mag_a[i],
                complete_tracklets.mag_b[i],
                complete_tracklets.mag_c[i],
            ]
        )
        mag_max = np.max(
            [
                complete_tracklets.mag_a[i],
                complete_tracklets.mag_b[i],
                complete_tracklets.mag_c[i],
            ]
        )

        if (mag_max - mag_min) < max_mag_variance:
            coordA = SkyCoord(
                ra=complete_tracklets.ra_a[i],
                dec=complete_tracklets.dec_a[i],
                unit=(u.deg, u.deg),
                distance=70 * u.kpc,
            )
            coordB = SkyCoord(
                ra=complete_tracklets.ra_b[i],
                dec=complete_tracklets.dec_b[i],
                unit=(u.deg, u.deg),
                distance=70 * u.kpc,
            )
            coordC = SkyCoord(
                ra=complete_tracklets.ra_c[i],
                dec=complete_tracklets.dec_c[i],
                unit=(u.deg, u.deg),
                distance=70 * u.kpc,
            )
            # print("coords", coordA,coordB,coordC)
            lenAB = coordA.separation_3d(coordB)
            lenBC = coordB.separation_3d(coordC)
            lenCA = coordC.separation_3d(coordA)
            cosine_angle = ((lenAB**2) + (lenBC**2) - (lenCA**2)) / (
                2 * lenAB * lenBC
            )
            angle = np.degrees(np.arccos(cosine_angle))

            # print("Angle between points A, B, and C:", angle, "degrees")
            if angle.value < min_tracklet_angle:
                complete_tracklets.drop(index=[i], inplace=True)
            else:
                angle_array.append(angle)
                mag_array.append(mag_max - mag_min)
                mag_min_array.append(
                    mag_max
                )  # because the "minimum" mag you want is the faintest one

        else:
            complete_tracklets.drop(index=[i], inplace=True)

    complete_tracklets.reset_index(inplace=True)
    complete_tracklets["angle"] = angle_array
    complete_tracklets["mag_diff"] = mag_array
    complete_tracklets["mag_min"] = mag_min_array

    print("Initial tracklet screening of", len(complete_tracklets), "complete.")
    sys.stdout.flush()  # print out everything before running FindOrb
    # now filter based on findorb
    tracklet_features = "outputs/tracklet_features" + night + ".txt"
    trackletfilename = "outputs/tracklets_" + night + ".txt"

    for i in range(len(complete_tracklets)):
        tracklet_id = "cn" + str(tracklet_num).rjust(5, "0")
        tracklet_num += 1

        # ratio of velocity between points a-b and points b-c
        ab_bc_vratio = (complete_tracklets.ab_dist[i] / (b_time - a_time)) / (
            complete_tracklets.bc_dist[i] / (c_time - b_time)
        )

        coordA = SkyCoord(
            ra=complete_tracklets.ra_a[i],
            dec=complete_tracklets.dec_a[i],
            unit=(u.deg, u.deg),
            distance=70 * u.kpc,
        )
        coordB = SkyCoord(
            ra=complete_tracklets.ra_b[i],
            dec=complete_tracklets.dec_b[i],
            unit=(u.deg, u.deg),
            distance=70 * u.kpc,
        )
        coordC = SkyCoord(
            ra=complete_tracklets.ra_c[i],
            dec=complete_tracklets.dec_c[i],
            unit=(u.deg, u.deg),
            distance=70 * u.kpc,
        )
        sky_sep = coordA.separation(coordC)
        sky_sep_arcs = sky_sep.arcsecond

        formatted_data = "     "
        formatted_data += "{}".format(tracklet_id) + "  C"
        formatted_data += "{}".format(a_time.strftime("%Y %m %d")) + "."
        formatted_data += "{:1}".format(decimal_time_a[1][:5]) + " "
        formatted_data += (
            coordA.to_string(style="hmsdms", pad=True, sep=" ", precision=2)
            + "         "
        )
        formatted_data += "{:.1f}".format(complete_tracklets.mag_a[i]) + "   "
        formatted_data += (
            complete_tracklets.band[i]
            + "    "
            + str(complete_tracklets.observatory_code[i])
            + "\n"
        )

        formatted_data += "     "
        formatted_data += "{}".format(tracklet_id) + "  C"
        formatted_data += "{}".format(b_time.strftime("%Y %m %d")) + "."
        formatted_data += "{:1}".format(decimal_time_b[1][:5]) + " "
        formatted_data += (
            coordB.to_string(style="hmsdms", pad=True, sep=" ", precision=2)
            + "         "
        )
        formatted_data += "{:.1f}".format(complete_tracklets.mag_b[i]) + "   "
        formatted_data += (
            complete_tracklets.band[i]
            + "    "
            + str(complete_tracklets.observatory_code[i])
            + "\n"
        )

        formatted_data += "     "
        formatted_data += "{}".format(tracklet_id) + "  C"
        formatted_data += "{}".format(c_time.strftime("%Y %m %d")) + "."
        formatted_data += "{:1}".format(decimal_time_c[1][:5]) + " "
        formatted_data += (
            coordC.to_string(style="hmsdms", pad=True, sep=" ", precision=2)
            + "         "
        )
        formatted_data += "{:.1f}".format(complete_tracklets.mag_c[i]) + "   "
        formatted_data += (
            complete_tracklets.band[i]
            + "    "
            + str(complete_tracklets.observatory_code[i])
            + "\n"
        )
        # print(formatted_data)

        trackletFound = "n"  # this will change later if you use findorb and the tracklet is less than
        # the residual
        res = "NaN"  # this is the findorb residual; will stay 'NaN' if you don't run findorb.
        if findorb_check == "y":
            findOrbTxt = open("/projectnb/ct-ast/findPOTATOs/fo.txt", "w")
            findOrbTxt.writelines(formatted_data)
            findOrbTxt.close()

            trackletFound, res = find_orb(
                Maximum_residual, nullResid=True, MOIDLim=True
            )
            # print("results of find_orb:",trackletFound,res)

        if (findorb_check == "y" and trackletFound == "y") or findorb_check == "n":
            print("Candidate tracklet found!\n", formatted_data)
            if exists(trackletfilename):
                with open(trackletfilename, "a", encoding="utf-8") as f:
                    f.write(formatted_data)
                    f.close
            else:
                with open(trackletfilename, "x", encoding="utf-8") as f:
                    f.write(formatted_data)
                    f.close

            if exists(tracklet_features):
                with open(tracklet_features, "a", encoding="utf-8") as f:
                    f.write(
                        tracklet_id
                        + ","
                        + str(complete_tracklets.angle[i].value)
                        + ","
                        + str(complete_tracklets.mag_diff[i])
                        + ","
                        + str(complete_tracklets.mag_min[i])
                        + ","
                        + str(res)
                        + ","
                        + str(sky_sep_arcs)
                        + ","
                        + str(ab_bc_vratio)
                        + "\n"
                    )
                    f.close
            else:
                with open(tracklet_features, "x", encoding="utf-8") as f:
                    f.write(
                        "tracklet_id,angle_deg,mag_diff,mag_min,residual,sky_sep,ab_bc_vratio\n"
                    )
                    f.write(
                        tracklet_id
                        + ","
                        + str(complete_tracklets.angle[i].value)
                        + ","
                        + str(complete_tracklets.mag_diff[i])
                        + ","
                        + str(complete_tracklets.mag_min[i])
                        + ","
                        + str(res)
                        + ","
                        + str(sky_sep_arcs)
                        + ","
                        + str(ab_bc_vratio)
                        + "\n"
                    )
                    f.close
            if print_thumbs == "y":
                # print("Saving thumbnails for tracklet:"tracklet_id)
                save_thumbnails_ml(
                    fits_frame_a,
                    tracklet_id,
                    "a",
                    complete_tracklets.x_a[i],
                    complete_tracklets.y_a[i],
                    telescope_image_a,
                )
                save_thumbnails_ml(
                    fits_frame_b,
                    tracklet_id,
                    "b",
                    complete_tracklets.x_b[i],
                    complete_tracklets.y_b[i],
                    telescope_image_b,
                )
                save_thumbnails_ml(
                    fits_frame_c,
                    tracklet_id,
                    "c",
                    complete_tracklets.x_c[i],
                    complete_tracklets.y_c[i],
                    telescope_image_c,
                )

            if export_ades == "y":
                # update keys in obsData dictionary for all three points
                # Coordinate A
                ades_obs_dict["trkSub"] = tracklet_id
                ades_obs_dict["obsTime"] = str(a_time.isot) + "Z"
                ades_obs_dict["ra"] = coordA.ra.deg
                ades_obs_dict["dec"] = coordA.dec.deg
                ades_obs_dict["mag"] = "{:.1f}".format(complete_tracklets.mag_a[i])

                if xml_tracklet_found == "n":  # first tracklet found, write header
                    xml_tracklet_found = "y"
                    XMLElement, ades, obsData = generate_xml(
                        xml_filename, ades_dict, ades_obs_dict
                    )

                else:  # update the existing xml
                    XMLElement, ades, obsData = update_xml(
                        XMLElement, ades, obsData, ades_obs_dict
                    )

                # Coordinate B
                ades_obs_dict["obsTime"] = str(b_time.isot) + "Z"
                ades_obs_dict["ra"] = coordB.ra.deg
                ades_obs_dict["dec"] = coordB.dec.deg
                ades_obs_dict["mag"] = "{:.1f}".format(complete_tracklets.mag_b[i])
                XMLElement, ades, obsData = update_xml(
                    XMLElement, ades, obsData, ades_obs_dict
                )

                # Coordinate C
                ades_obs_dict["obsTime"] = str(c_time.isot) + "Z"
                ades_obs_dict["ra"] = coordC.ra.deg
                ades_obs_dict["dec"] = coordC.dec.deg
                ades_obs_dict["mag"] = "{:.1f}".format(complete_tracklets.mag_c[i])
                XMLElement, ades, obsData = update_xml(
                    XMLElement, ades, obsData, ades_obs_dict
                )

        if findorb_check == "y" and trackletFound == "n":  # drop it
            print("tracklet rejected")
            complete_tracklets.drop(index=[i], inplace=True)

    # save stats
    now = datetime.now()
    yearmonthday = now.strftime("%Y%m%d")
    outputname = "outputs/o_linking_" + yearmonthday + ".csv"
    run_time = str(now - start_time)
    num_sources = str(len(a) + len(b) + len(c))
    num_tracklets_prescreen = str(len(complete_tracklets))

    if not exists(outputname):
        f = open(outputname, "x", encoding="utf-8")
        f.write(
            "filea,fileb,filec,date_corrected,time_corrected,run_time_s,num_sources,num_tracklets_prescreen\n"
        )
    else:
        f = open(outputname, "a", encoding="utf-8")
    f.write(
        file_a
        + ","
        + file_b
        + ","
        + file_c
        + ","
        + datetime.today().strftime("%Y-%m-%d")
        + ","
        + datetime.today().strftime("%H:%M:%S")
        + ","
        + run_time
        + ","
        + num_sources
        + ","
        + num_tracklets_prescreen
        + "\n"
    )
    f.close

if export_ades == "y":
    # write the ADES xml to file
    tree = XMLElement.ElementTree(ades)
    xml_string = minidom.parseString(XMLElement.tostring(ades)).toprettyxml()
    with open(xml_filename, "w", encoding="UTF-8") as files:
        files.write(xml_string)
