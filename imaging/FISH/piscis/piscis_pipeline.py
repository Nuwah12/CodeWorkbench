# standard libs
import os
import re
from pathlib import Path
import numpy as np
# plotting
import plotly.graph_objects as go
from PIL import Image
import base64
from io import BytesIO
# timing
from timeit import default_timer as timer
from datetime import timedelta
# image I/O
import nd2
import tifffile
# spot calling
import piscis
# argument input / logging
import argparse
import yaml
import logging

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger("piscis") # init the logging object

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
        logger.error(f" {e}: Couldn't find specified settings file {args.settingsFile}. Check that it exists and re-specify.")
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
        #print(f"{i},{s}")
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
    Returns a maximum projection of the input image over all Z slices
    """
    return np.max(img, axis=0)

def _call_spots_piscis(piscis_obj, img, threshold, max_proj=True, stack=True):
    """
    Call spots on an image using the Piscis model.
    Arguments:
        piscis_obj (piscis)     - an instance of a piscis object
        img - (ndarray)         - the image to be analyzed
        threshold (int/float)   - piscis threshold model parameter
    """
    logger.info("Starting spot detection.")
    start = timer()
    pred = piscis_obj.predict(img, threshold=1, stack=stack)
    
    elapsed = timedelta(seconds=timer()-start) # time
    if max_proj:
        num_spots = len(pred)
    else:
        num_spots = len(np.concatenate(pred))
    
    logger.info(f"Finished calling spots; {num_spots} spots found in {elapsed}")

    return pred

# Normalize and convert Z-slice to uint8
def _normalize_to_uint8(slice_2d):
    p_min, p_max = np.percentile(slice_2d, (1, 99))
    norm = np.clip((slice_2d - p_min) / (p_max - p_min), 0, 1)
    return (norm * 255).astype(np.uint8)


def _interactive_plot(img, spots, mode, outf, neighborhoods=None):
    """
    Create an interactive plot using plotly and save to html
    Arguments:
        img (ndarray) - the image to be plotted
        spots (array) - the spots to be plotted
        mode (str) - one of ["stack", "max"], determines what is plotted
    """
    logger.info(f"Plotting detected spots with mode {mode}")
    if mode == "stack":
        Z = img.shape[0]  # total number of z slices
        frames = []
        for z in range(Z):
            img_slice = _normalize_to_uint8(img[z]) # get the current z slice and normalize it to [0, 255]

            coords_z = spots[spots[:, 0] == z]  # filter spots at z
            y = coords_z[:, 1]   # y coord
            x = coords_z[:, 2]   # x coord

            frame = go.Frame(
                data=[
                    go.Heatmap(  # show image
                        z=img_slice,
                        colorscale='gray',
                        showscale=False),
                    go.Scatter(  # show spots
                        x=x,
                        y=y,
                        mode='markers',
                        marker=dict(color='red', size=5),
                        name='Spots')],
                name=str(z))

            frames.append(frame)
        
        fig = go.Figure(
            data=[
                go.Heatmap(
                    z=_normalize_to_uint8(img[z]),
                    colorscale='gray',
                    showscale=False),
                go.Scatter(
                    x=coords_z[:, 2],  # X
                    y=coords_z[:, 1],  # Y
                    mode="markers",
                    marker=dict(color='red', size=5),
                    name='Spots')],
            frames=frames)

        # Add slider and play buttons
        fig.update_layout(
            sliders=[{
            "steps": [
            {"method": "animate", "args": [[str(z)], {"mode": "immediate"}], "label": f"Z={z+1}"} for z in range(Z)],
                "currentvalue": {"prefix": "Slice: "}}],
            height=700,
            width=700,
            title="Z-stack Spot Viewer")

        fig.update_yaxes(autorange="reversed")  # Important for image-style orientation
        fig.write_html(outf)

    if mode == "max":
        vmin, vmax = np.percentile(img, (1, 99))  # or (0.5, 99.5) for more aggressive stretch
        img_slice = _normalize_to_uint8(img)
        y = spots[:, 1]
        x = spots[:, 2]

        fig = go.Figure()

        fig.add_trace(go.Heatmap(
            z=img,
            colorscale='gray',
            showscale=False,
            zmin=vmin,
            zmax=vmax
        ))

        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='markers',
            marker=dict(color='red', size=5),
            name='Spots'
        ))

        fig.update_layout(
            title="All spots",
            height=700,
            width=700
        )
        fig.update_yaxes(autorange="reversed")
        fig.write_html(outf)

def _pair_mask_img(img_paths, mask_paths):
    """
    Pair image and mask files.
    """
    mapping = {}
    # First map Location_XX : mask filename
    # Then loop through image files, split string to get Location_XX, and index mask dict with it to get path
    # Add to dictionary {image_filename : mask_filename}
    loc_map = {}
    for m in mask_paths:
        f = m.split("/")[len(m.split("/"))-1]
        loc = f.split("_")
        loc_map[f"{loc[0]}_{loc[1]}"] = m # map ID : mask path
    
    for i in img_paths:
        f = i.split("/")[len(i.split("/"))-1]
        loc = f.split("_")
        mapping[i] = loc_map[f"{loc[0]}_{loc[1]}"]

    return mapping

def _mask_spots(mask, spots):
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
    channels = settings["channels"]

    mask_dir = settings["mask_dir"]
    mask_files = [os.path.join(mask_dir, f) for f in os.listdir(mask_dir)]

    model = piscis.Piscis(model_name=settings["model"]) # make piscis obj to reuse
    stack = settings["stack"]
    threshold = settings["piscis_thresh"] # piscis threshold parameter
    scale = settings["piscis_scale"]
    min_dist = settings["piscis_min_distance"]

    plot_max = settings["plot_max"]
    plot_z = settings["plot_z"]
    plot_out_dir = settings["plot_out_dir"]
    
    mask_mapping = _pair_mask_img(img_files, mask_files)
    # Main loop
    for n,i in enumerate(img_files): # TODO: Now, we want to loop througn the {image : mask} dictionary
        logger.info(f"Starting image {i}")
        #logger.info(f"Using mask file {mask_files[n]}")
        exit(0)
        j = _read_img(i, imgtype) # read image to np nd array
        jname = i.split("/")[len(i.split("/"))-1]
        
        if len(channels) > 1: # if there is more than one channel, subset the image with the one we want 
            channelDim = _get_channel_dim(j, len(channels)) # guesstimate the channel dimension
            channelIdx = channels.index(callChannel) # get the index of the spot calling channel

            index = [slice(None)] * j.ndim
            index[channelDim] = channelIdx
        
            j = j[tuple(index)] # index the array to just the channel for spot calling
    
        # call spots
        pred_spots = _call_spots_piscis(model, j, threshold, stack=stack)

        # plot
        if plot_max:
            _interactive_plot(_max_proj_image(j), pred_spots, mode="max", outf=f"{plot_out_dir}/{jname}_interactivePlot_maxProj_allSpots.html")
        if plot_z:
            _interactive_plot(j, pred_spots, mode="stack", outf=f"{plot_out_dir}/{jname}_interactivePlot_zStack_allSpots.html")
        
        # output
        if not Path(settings["spot_out"]).is_dir():
            Path(settings["spot_out"]).mkdir()
        np.savetxt(f"{settings["spot_out"]}/{jname}_allCalledSpots.tsv", pred_spots, delimiter="\t")
      
if __name__ == "__main__":
    main()
