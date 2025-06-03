### Run Piscis in a Pipeline
The python script `piscis_pipeline.py` will run the [piscis](https://github.com/zjniu/Piscis) spot-calling algorithm in a high-throughput manner givn a dataset of smFISH microscopy images in `nd2` or `tiff` format. \
Spots can be called on multi or single channel images as well as stacks of images or flat images. \

#### File naming
The pipeline expects both image and mask files to contain the prefix `Location_XX`, i.e. `Location_01`, `Location_02`, etc. The numbers need not be increasing or in numerical order, but this prefix **_must_** match between image files and their corresponding mask files. \
Example:
```
Location_01_DAPI_MaxProj_cp_masks.tif
Location_02_DAPI_MaxProj_cp_masks.tif
Location_03_DAPI_MaxProj_cp_masks.tif
...
```

#### Usage
Before running the pipeline, I suggest installing and activating the conda environment in the above directory - `imagingEnv.yml`. It contains all packages in correct versions needed to run it. \
To install, copy the `yml` file to your local machine and run: 
```
conda env create -f imgagingEnv.yml
conda activate bigfish
```
To execute the pipeline, run
```
python3 piscis_pipeline.py [-h] settings.yml
```
The script takes one argument, the `settings.yml` file. All parameters related to spot calling and image processing are in this file - no other rcommand line arguments are needed. \
The available arguments include:
##### Image arguments
* `image_dir`: Path to directory with image files. Please ensure the files to be analyzed are in one directory.
* `image_type`: Type of image to be analyzed. Currently supported options are `tif`/`tiff` and `nd2`
* `channels`: An array with one entry for each fluoresence channel in the image. The dimension containing channel information is inferred based on the size of this list.
  * **NOTE:** Please ensure the ordering of channels here is the same as it is in your image file, i.e. in the below example, channel 0 is "DAPI", channel 1 is "Cy3", and channel 2 is "A647".
##### Spot calling / de-duplication
* `piscis_thresh`: The `threshold` parameter for the piscis model. It should not be changed much. Check [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC10862914/) for more information.
* `piscis_scale`: Piscis `scale` parameter
* `piscis_min_distance`: Piscis `minimum_distance` parameter
* `model`: The (string) representing the piscis model to be used. Only change if results are not looking good.
* `spot_channel`: The name of the channel from the above list to be used to call spots.
##### Plotting
* `plot_max`: Boolean; Make an interactive plot of all dots on the max-projected image.
* `plot_z`: Boolean; Interactive plot of all dots on their respective Z-slices.
* `plot_out_dir`: Directory to output interactive plots to in `html` format.

A properly formatted `settings.yml` would look like:
```
image_dir: "/path/to/image.tif"              
image_type: "tif"             
spot_out: "./spots"   

channels:                   
 - "A647"
 - "DAPI"
 - "Cy3"
spot_channel: "A647"        
stack: true                   
piscis_thresh: 1              
piscis_scale: 1                
piscis_min_distance: 1        

model: "20230905"         

plot_max: true             
plot_z: true                    
plot_out_dir: "."        
```
