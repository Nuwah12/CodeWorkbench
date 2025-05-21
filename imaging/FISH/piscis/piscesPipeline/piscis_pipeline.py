import os
import numpy as np
from sklearn.neighbors import NearestNeighbors

from timeit import default_timer as timer
from datetime import timedelta

import nd2
import tifffile

import piscis

import argparse
import yaml
import logging

logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__) # init the logging object

def _parse_args():
    """
    Parse input argument and the settings.yml file.
    """
    parser = argparse.ArgumentParser(prog="Piscis Pipeline",
                                     description="Call spots in smFISH data with Piscis")
    parser.add_argument("settingsFile", type=str, help="Path to settings.yml") # only one argument - the settings file
   
    args = parser.parse_args()
    
    try:
        settings = yaml.safe_load(open(args.settingsFile)) # load the yaml
    except FileNotFoundError as e:
        logger.error(f"{e}: Couldn't find specified settings file {args.settingsFile}. Check that it exists and re-specify.")
        exit(1)
    
    return settings

def _read_img(path, img_type):
    """
    Read in an image file with the correct library.
    Arguments:
        path (str)              - path to the image file
        img_type (str)          - image filetype ["tif", "tiff", "nd2"]
    """
    if img_type == "tif" or img_type == "tiff":
        return tifffile.imread(path)
    elif img_type == "nd2":
        if ".nd2" in path:
            return nd2.imread(path)
        else:
            raise ValueError("File does not end in .nd2, check filetype and extension.")

def _get_channel_dim(img, numChannels):
    """
    Estimate which dimension along the ndarray contains the channel information.
    This is done by simply looking for the dimension which has size equal to the number of dimensions provided in settings.yml
    Arguments:
        img (ndarray)           - img to be analyzed
        numChannels (int)       - the total number of channels in img
    """
    for i, s in enumerate(img.shape):
        if s != numChannels:
            continue
        else:
            logger.info(f"Inferring dim. {i} (size = {s}) is the channel dimension")
            return i

def _checks(args):
    """
    Sanity checks for arguments passed via settings.yml
    Arguments:
        args (dict)             - dictionary of arguments parsed from settings.yml
    """
    # img dir exists
    if not os.path.exists(args["image_dir"]):
        raise FileNotFoundError(f"Couldn't find image directory {args['image_dir']}")
        exit(1)
    # channel to call spots on is in list of channels
    if not args["spot_channel"] in args["channels"]:
        raise ValueError(f"Specified channel to call spots on ({args["spot_channel"]}) is not present in provided channel names ({args["channels"]}). Please check the provided names and re-specify.")
        exit(1)

def _call_spots_piscis(piscis_obj, img, threshold):
    """
    Call spots on an image using the Piscis model.
    Arguments:
        piscis_obj (piscis)     - an instance of a piscis object
        img - (ndarray)         - the image to be analyzed
        threshold (int/float)   - piscis threshold model parameter
    """
    logger.info("Starting spot detection...")
    start = timer()
    pred = piscis_obj.predict(img, threshold=1)

    elapsed = timedelta(seconds=timer()-start) # time
    logger.info(f"Finished calling spots; {elapsed} elapsed")

    return pred

def _dedup_spots(spots, img, threshold):
    """
    Deduplicate called spots via a nearest neighbors method
    Arguments:
        spots (list of lists)   - list of called spots for image img, with 1 sublist for each z slice
        img (ndarray)           - image to be deduplicated
        threshold (int/float)   - diameter in pixels for finding neighborhoods
    """

def main():
    """
    Main function. Gathers image files and arguments to perform spot calling.
    """
    settings = _parse_args() # returns settings as dict {key: value}
    
    img_files = [os.path.join(settings["image_dir"], f) for f in os.listdir(settings["image_dir"])]
    imgtype = settings["image_type"]
    callChannel = settings["spot_channel"]

    _checks(settings) # run checks on supplied argument values
    
    model = piscis.Piscis(model_name=settings["model"]) # make piscis obj to reuse
    threshold = settings["piscis_thresh"] # piscis threshold parameter

    # Main loop
    for i in img_files:
        logger.info(f"Starting image {i}")
        j = _read_img(i, imgtype) # read image to np nd array
        
        channelDim = _get_channel_dim(j, len(settings["channels"])) # guesstimate the channel dimension
        channelIdx = settings["channels"].index(settings["spot_channel"]) # get the index of the spot calling channel

        index = [slice(None)] * j.ndim
        index[channelDim] = channelIdx
        
        j = j[tuple(index)] # index the array to just the channel for spot calling
        
        pred_spots = _call_spots_piscis(model, j, threshold)
        #print(j[tuple(index)].shape) 
        


if __name__ == "__main__":
    main()
