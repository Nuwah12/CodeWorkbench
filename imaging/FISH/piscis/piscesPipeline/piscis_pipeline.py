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

# Normalize and convert Z-slice to uint8
def _normalize_to_uint8(slice_2d):
    p_min, p_max = np.percentile(slice_2d, (1, 99))
    norm = np.clip((slice_2d - p_min) / (p_max - p_min), 0, 1)
    return (norm * 255).astype(np.uint8)


def _interactive_plot(img, spots, mode="all", outf="."):
    """
    Create an interactive plot using plotly and save to html
    Arguments:
        img (ndarray) - the image to be plotted
        spots (array) - the spots to be plotted
        mode (str) - one of ['all', 'max', 'neighbors'], determines what is plotted
            'all' - plots all dots on their respective z slices
            'dedup' - plots deduplicated spots on a max z projection of img
            'neighbors' - plots all dots on a max projection of img with lines denoting dots that were called as one neighborhood
    """
    lines_x = []
    lines_y = []
    for i, nbrs in enumerate(neighbors):
        for j in nbrs:
            if i >= j:  # avoid duplicates
                continue
    lines_x.extend([x[i], x[j], None])
    lines_y.extend([y[i], y[j], None])
    
    img_slice = normalize_to_uint8(loc1_untr_maxProj)

    pil_img = Image.fromarray(img_slice)
    buffer = BytesIO()
    pil_img.save(buffer, format="PNG")
    encoded = base64.b64encode(buffer.getvalue()).decode()

    # Create figure with image background
    fig = go.Figure()

    # Add lines between neighbors
    fig.add_trace(go.Scatter(
        x=lines_x,
        y=lines_y,
        mode='lines',
        line=dict(color='blue', width=1),
        name='Neighbor Links'))

    # Add spot markers
    fig.add_trace(go.Scatter(
        x=x,
        y=y,
        mode='markers',
        marker=dict(color='rgba(255,0,0,1)', size=5),
        name='Spots'))

    # Overlay the image
    fig.update_layout(
        images=[dict(
            source=f'data:image/png;base64,{encoded}',
            xref="x", yref="y",
            x=0, y=0,
            sizex=img_slice.shape[1],  # X-axis size (width)
            sizey=img_slice.shape[0],  # Y-axis size (height)
            sizing="stretch",
            opacity=1.0,
            layer="below")],
        height=700,
        width=700,
        title=f"Max-Z projection with Neighbors")

    fig.update_yaxes(autorange='reversed')
    fig.write_html('../piscis/loc01_ICOS_untr_spotNeighborhoods.html')


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

    plot_all = settings["plot_all"]
    plot_dedup = settings["plot_dedup"]
    plot_neighborhoods = settings["plot_neighborhoods"]
    plot_out_dir = settings["plot_out_dir"]

    # Main loop
    for i in img_files:
        logger.info(f"Starting image {i}")
        j = _read_img(i, imgtype) # read image to np nd array
        jname = i.split("/")[len(i.split("/"))-1]
         
        channelDim = _get_channel_dim(j, len(settings["channels"])) # guesstimate the channel dimension
        channelIdx = settings["channels"].index(settings["spot_channel"]) # get the index of the spot calling channel

        index = [slice(None)] * j.ndim
        index[channelDim] = channelIdx
        
        j = j[tuple(index)] # index the array to just the channel for spot calling
        
        # call spots
        pred_spots = _call_spots_piscis(model, j, threshold)
        
        # de-duplicate the spots called over all z slices
        dedup_spots = _dedup_spots(pred_spots, j, settings["dedup_radius"])
        
        if plot_all:
            _interactive_plot(j, pred_spots, mode="all", outf=f"{plot_out_dir}/{jname}_interactivePlot_allDots_allZslices.html")
        if plot_dedup:
            _interactive_plot(_max_proj_image(j), mode="dedup", outf=f"{plot_out_dir}/{jname}_interactivePlot_dedupDots_maxProj.html")
        if plot_neighborhoods:
            _interactive_plot(_max_proj_image(j), mode="neighbors", outf=f"{plot_out_dir}/{jname}_interactivePlot_allDots_neighborhoods_maxProj.html")
        exit(0)        


if __name__ == "__main__":
    main()
