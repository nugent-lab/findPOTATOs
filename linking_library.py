import os
import subprocess
import re  # regular expressions, used to search for mean residuals in Find_orb output files
from time import sleep
import pandas as pd
import numpy as np
import matplotlib as plt

# from datetime import datetime
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, Angle, Distance
from sklearn.neighbors import BallTree
from PIL import Image


def find_orb(maxResidual, nullResid=True, MOIDLim=False):
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

    Args:
        maxResidual: float, maximum residual allowed for tracklet to be approved

    Returns:
        trackletFound: string, equals 'yes' if passed
        res: associated residual with found tracklet
    """
    trackletFound = "n"
    elements_path = "~/.find_orb/elements.txt"  # for mac
    if os.path.exists(os.path.expanduser(elements_path)):
        os.remove(os.path.expanduser(elements_path))
    # this line works on mac & some unix installs but not the MGHPCC
    # sp = Popen(['cd ~/.find_orb\n~/find_orb/find_orb/fo fo.txt -c'], shell=True)
    # this line is for the MGHPCC. Either way, you need the directory where your fo files are
    # the subprocess module reacts poorly to the supercomputer.
    # os.system('fo fo.txt -c')
    findorb_call = "fo fo.txt -c"
    sp = subprocess.call(findorb_call, shell=True)
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
        # numObs = sum(1 for line in open(os.path.expanduser("~/.find_orb/fo.txt")))
        numObs = sum(
            1
            for line in open(os.path.expanduser("/projectnb/ct-ast/findPOTATOs/fo.txt"))
        )

        # save all inputs to find_orb
        open("outputs/AllPotentialTracklets.txt", "a+").writelines(
            [l for l in open(
                    os.path.expanduser("/projectnb/ct-ast/findPOTATOs/fo.txt")
                ).readlines()
            ]
        )
        for line in open(os.path.expanduser(elements_path)):
            li = line.strip()
            if not li.startswith("#"):
                open("outputs/AllPotentialTracklets.txt", "a").writelines(line.rstrip())
                open("outputs/AllPotentialTracklets.txt", "a").writelines("\n")
        open("outputs/AllPotentialTracklets.txt", "a").writelines("\n\n")

        resCheck = False
        for line in open(os.path.expanduser(elements_path)):
            match = re.search('mean residual (\d+)".(\d+)', line)
            match2 = re.search("MOID: (\d+).(\d+)", line)
            if match:
                res = int(match.group(1)) + float(("0." + match.group(2)))
                if nullResid:
                    if (res < maxResidual) & (res > 0):  # maxResidual in "
                        resCheck = True
                    else:
                        resCheck = False
                else:
                    if res < maxResidual:  # maxResidual in "
                        resCheck = True
                    else:
                        resCheck = False
            if match2:
                if MOIDLim:
                    MOID = int(match2.group(1)) + float(("0." + match2.group(2)))
                    if MOID > MOIDLim:
                        print("MOID:", MOID, " exceeds MOIDLim:", MOIDLIM)
                        break
        if resCheck:
            trackletFound = "y"

    else:
        print("Could not open file", os.path.expanduser(elements_path))
    return trackletFound, res


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
    # convert threshold (arsec) to degrees, then to radians
    thresh_rad = np.radians(thresh / 3600)
    # print("thresh_rad",thresh_rad)

    # intialize output dataframes
    # gonna delete the duplicates from the df?_moving dataframes
    df1_moving = df1.copy()
    df2_moving = df2.copy()
    df3_moving = df3.copy()

    # compare 1 to 2
    tree1 = BallTree(df1[["ra_rad", "dec_rad"]], leaf_size=5, metric="haversine")
    indicies = tree1.query_radius(df2[["ra_rad", "dec_rad"]], r=thresh_rad)
    for i in range(len(indicies)):
        for j in range(len(indicies[i])):
            try:  # using a try, since there's gonna be duplicates between comparing a to b, b to c, etc
                df1_moving.drop(index=indicies[i][j], inplace=True)
            except:
                pass

    # compare 1 to 3
    indicies = tree1.query_radius(df3[["ra_rad", "dec_rad"]], r=thresh_rad)
    for i in range(len(indicies)):
        for j in range(len(indicies[i])):
            try:
                df1_moving.drop(index=indicies[i][j], inplace=True)
            except:
                pass

    # compare 2 to 1
    tree2 = BallTree(df2[["ra_rad", "dec_rad"]], leaf_size=5, metric="haversine")
    indicies = tree2.query_radius(df1[["ra_rad", "dec_rad"]], r=thresh_rad)
    for i in range(len(indicies)):
        for j in range(len(indicies[i])):
            try:
                df2_moving.drop(index=indicies[i][j], inplace=True)
            except:
                pass

    # compare 2 to 3
    indicies = tree2.query_radius(df3[["ra_rad", "dec_rad"]], r=thresh_rad)
    for i in range(len(indicies)):
        for j in range(len(indicies[i])):
            try:
                df2_moving.drop(index=indicies[i][j], inplace=True)
            except:
                pass

    # compare 3 to 1
    tree3 = BallTree(df3[["ra_rad", "dec_rad"]], leaf_size=5, metric="haversine")
    indicies = tree3.query_radius(df1[["ra_rad", "dec_rad"]], r=thresh_rad)
    for i in range(len(indicies)):
        for j in range(len(indicies[i])):
            try:
                df3_moving.drop(index=indicies[i][j], inplace=True)
            except:
                pass

    # compare 3 to 2
    indicies = tree3.query_radius(df2[["ra_rad", "dec_rad"]], r=thresh_rad)
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

    print(
        "Percentage of sources remaining from original dataframes df1, df2, df3:",
        round(len(df1_moving) / len(df1), 2),
        round(len(df2_moving) / len(df2), 2),
        round(len(df3_moving) / len(df3), 2),
    )
    return df1_moving, df2_moving, df3_moving


def mpc_reader(file_name):
    """
    Reads in 80-char observation file.
    Args:
        file_name: string of filename to be read
    Returns:
        dat: dataframe of file
    """
    txt_file = open(file_name)
    track_id = []
    date = []
    ra = []
    dec = []
    mag = []
    # run through until file ends
    while 1:
        # pull 50 lines from memory at a time
        lines = txt_file.readlines(50)
        # break if at the end of file
        if not lines:
            break
        # look through each line
        for line in lines:
            # print(line)
            track_id.append(line[:12].strip(" "))
            # get the date into the right format for astropy
            prelim_date = line[15:25].replace(" ", "-")
            t = Time(prelim_date, format="isot", scale="utc")
            frac = float(line[25:31])
            newt = Time(t.mjd + frac, format="mjd", scale="utc")
            date.append(newt.mjd)
            c = SkyCoord(line[32:43], line[45:55], unit=(u.hourangle, u.deg))
            ra.append(c.ra.degree)
            dec.append(c.dec.degree)
            mag.append(line[65:70])
            # print(newt.mjd, c.ra.degree, c.dec.degree)
    d = {"tracklet_id": track_id, "date_mjd": date, "RA": ra, "Dec": dec, "mag": mag}
    dat = pd.DataFrame(data=d)
    return dat


def save_thumbnails(fits_frame, tracklet_id, abc, x_pos, y_pos, telescope_image):
    """
    Saves thumbnails in human-readable format.
    Args:
        fits_frame: string, name of image
        tracklet_id: string, unique id of tracklet
        abc: string, 'a', 'b', 'c' to indicate order of source in tracklet
        x_pos: float, RA of source you want thumbnail of
        y_pos: float, Dec of source you want thumbnail of
        telescope_image: array, actual image
    Returns:
        none
    """
    buffer = 20
    image_name = fits_frame.split(".")
    if y_pos < 0:
        y_pos = buffer
    if x_pos < 0:
        x_pos = buffer

    fig, ax = plt.subplots()
    m, s = np.mean(telescope_image), np.std(telescope_image)
    plt.xlim([int(x_pos) - (buffer), int(x_pos) + (buffer)])
    plt.ylim([int(y_pos) - (buffer), int(y_pos) + (buffer)])
    plt.axis("off")
    plt.margins(0, 0)
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    im = ax.imshow(
        telescope_image,
        interpolation="nearest",
        cmap="gray",
        vmin=m - s,
        vmax=m + s,
        origin="lower",
    )
    plt.savefig(
        "thumbs/thumb_" + image_name[0] + "_" + tracklet_id + "_" + abc + ".png",
        format="png",
    )
    plt.close("all")
    return


def save_thumbnails_ml(fits_frame, tracklet_id, abc, x_pos, y_pos, telescope_image):
    """
    Saves thumbnails for machine learning. Thumbs
    are very small and scaled from 0 to 255.
    Args:
        fits_frame: string, name of image
        tracklet_id: string, unique id of tracklet
        abc: string, 'a', 'b', 'c' to indicate order of source in tracklet
        x_pos: float, RA of source you want thumbnail of
        y_pos: float, Dec of source you want thumbnail of
        telescope_image: array, actual image
    Returns:
        none
    """
    buffer = 9
    image_name = fits_frame.split(".")
    left = int(x_pos) - (buffer - 1)
    right = int(x_pos) + (buffer)
    upper = int(y_pos) + (buffer - 1)
    lower = int(y_pos) - (buffer)

    if left < 0:
        left = 0
    elif upper < 0:
        upper = 0
    elif right > telescope_image.shape[1]:
        right = telescope_image.shape[1]
    elif lower >= telescope_image.shape[0]:
        lower = telescope_image.shape[0]

    # Crop the image and scale it to use all values.
    cropped_array = telescope_image[lower:upper, left:right]
    try:
        scaled_array = (
            (cropped_array - np.min(cropped_array))
            / (np.max(cropped_array) - np.min(cropped_array))
            * 255
        ).astype(np.uint8)
        cropped_image = Image.fromarray(scaled_array, mode="L")
        cropped_image.save(
            "thumbs/mlthumb_" + image_name[0] + "_" + tracklet_id + "_" + abc + ".png"
        )
    except ValueError:
        print("Error rescaling, here's the cropped_array", cropped_array)
    return
