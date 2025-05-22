### Run Piscis in a Pipeline
The python script `piscis_pipeline.py` will run the [piscis](https://github.com/zjniu/Piscis) spot-calling algorithm in a high-throughput manner givn a dataset of smFISH microscopy images in `nd2` or `tiff` format. \
By default, spots are called on each Z-stack slice independently. Afterwards, clusters of spots are de-duplicated using sklearn's `NearestNeighbors` module to identify the clusters, and the brightest spot from each neighborhood is retained.

#### Usage
```
python3 piscis_pipeline [-h] settings.yml
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
* `model`: The (string) representing the piscis model to be used. Only change if results are not looking good.
* `dedup_radius`: Radius (in pixels) for neighborhood identification in the de-duplication step.
* `spot_channel`: The name of the channel from the above list to be used to call spots.
* `call_projection`: **TBA**
##### Plotting
* `plot_all`: Boolean; Make an interactive plot of all dots on each Z-slice independently.
* `plot_neighborhoods`: Boolean; Interactive plot of all spots on a Z-max-projected image with lines showing neighborhood membership.
* `plot_dedup`: Boolean; Interactive plot of de-duplicated spots on max-projected image.
* `plot_out_dir`: Directory to output interactive plots to in `html` format.

A properly formatted `settings.yml` would look like:
```
image_dir: "/path/to/images"                          
image_type: "nd2"             

channels:                     
 - "DAPI"
 - "Cy3"
 - "A647"
spot_channel: "Cy3"    
call_projection: false  # TBA!      
piscis_thresh: 1              

model: "20230905"      

dedup_radius: 1          

plot_all: true        
plot_neighborhoods: true    
plot_dedup: true             
plot_out_dir: "."
```
