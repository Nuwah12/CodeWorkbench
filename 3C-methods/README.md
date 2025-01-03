## 3C Analysis methods for Faryabi Lab
Hosted here are analysis / processing methods for Chromatin Conformtion Capture (3C) experiments.
* [Faryabi Lab DockerHub](https://hub.docker.com/)
* [Juicer](https://github.com/aidenlab/juicer)
* [Cooler](https://github.com/open2c/cooler)
* [hic2cool](https://github.com/4dn-dcic/hic2cool)
* [Homer](http://homer.ucsd.edu/homer/interactions/)
* [GENOVA](https://github.com/robinweide/GENOVA)

### Overview of methods included in this repository
#### `/processing`
* **Workflows**:
  * 3C processing (`./Ingest_3C`) - uses `Juicer` for initial alignment/processing/normalization, then converts to `.(m)cool` using `hic2cool`. First step to processing 3C data.
* `make_homer_input.sh`: Convert Juicer's `merged_nodup.txt` (or, really any file in [AidenLab's "long" format](https://github.com/aidenlab/juicer/wiki/Pre#long-format)) to a [Homer-compatible format](http://homer.ucsd.edu/homer/interactions/HiCtagDirectory.html). This new input is then used to create a Tag Directory via Homer's `makeTagDirectory`.
#### `/visualization`
* A series of R scripts mainly using the `GENOVA` package to do basic analyses and plot 3C data. \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  `genova_loadCooler.R` - Load a .cool/.mcool file as a GENOVA-compatible object for plotting. \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `genova_viewMatrix.R` - Preconfigured GENOVA commands for viewing single & differential matrics in a variety of formats. \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `genova_aggregatePeakAnalysis.R` - Calculate APA values and plot. 
#### `/quantify_matrix`
* `FetchPixelFromMatrixGivenBedpe.py`: Given a matrix in `.cool/.mcool` form and a `BEDPE` file of interacting coordinates, extract raw and normalized contact counts.
