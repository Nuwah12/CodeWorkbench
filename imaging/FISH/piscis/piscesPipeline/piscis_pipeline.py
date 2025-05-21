import os
from pathlib import Path
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
logger = logging.getLogger("PiscisPipeline") # init the logging object

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
    assert Path(args["image_dir"]).exists(), f"Could not find supplied image directory '{args['image_dir']}'"
    # channel to call spots on is in list of channels
    assert args["spot_channel"] in args["channels"], f"Specified channel to call spots on ({args["spot_channel"]}) is not present in provided channel names ({args["channels"]}). Please check the provided names and re-specify."

def _max_proj_image(img):
    """
    Returns a masimum projection of the input image over all Z slices
    """
    return np.max(img, axis=0)

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
    logger.info(f"Finished calling spots; {len(np.concatenate(pred))} spots found in {elapsed}")

    return pred

def _dedup_spots(spots, img, radius):
    """
    Deduplicate called spots via a nearest neighbors method
    Arguments:
        spots (list of lists)   - list of called spots for image img, with 1 sublist for each z slice
        img (ndarray)           - image to be deduplicated
        threshold (int/float)   - diameter in pixels for finding neighborhoods
    """
    logger.info("Deduplicating called spots")
    all_spots = np.concatenate(spots) # unnest the list for all z slices

    # Nearest Neighbors 
    nn = NearestNeighbors(radius=radius)
    nn.fit(all_spots)
    neighborhoods = nn.radius_neighbors(all_spots, return_distance=False)

    # Max-project the image so we can find the brightest spots
    img_mp = _max_proj_image(img)

    # Now, select the brightest pixel in each neighborhood as the representative spot
    ## TODO: Fix deduplication scheme
    deduped_spots = []
    seen = set()

    for i, group in enumerate(neighborhoods):
        # Skip if we've already handled a neighbor in this group
        if any(j in seen for j in group):
            continue

        # Get intensities from max-projected image
        intensities = [
            img_mp[int(all_spots[j][0]), int(all_spots[j][1])] for j in group
        ]
        brightest_idx = group[np.argmax(intensities)]
        deduped_spots.append(all_spots[brightest_idx])
        seen.update(group)
    
    logging.info(f"{len(deduped_spots)} spots remain after deduplication.")
    return deduped_spots

def _interactive_plot(img, spots):
    ## TODO: implement plotly interactive plotting here
    pass

def main():
    """
    Main function. Gathers image files and arguments to perform spot calling.
    """
    settings = _parse_args() # returns settings as dict {key: value}
    _checks(settings)
    
    img_files = [os.path.join(settings["image_dir"], f) for f in os.listdir(settings["image_dir"])]
    imgtype = settings["image_type"]
    callChannel = settings["spot_channel"]
    
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

        dedup_spots = _dedup_spots(pred_spots, j, settings["dedup_radius"])
        #print(j[tuple(index)].shape) 
        exit(0)        


if __name__ == "__main__":
    main()
