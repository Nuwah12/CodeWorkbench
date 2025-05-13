### Spot Calling with BigFISH
Copy this repository to your local machine with `git clone`.\
Steps to start Jupyter notebook:
1. Navigate to this directory
2. On **Aries**, activate the bigfish virtual envornment: 
    * `source /opt/miniforge3/bin/activate` 
    * `mamba activate bigfish` 
    * If you see an error like: 
   ```
   'mamba' is running as a subprocess and can't modify the parent shell.
    Thus you must initialize your shell before using activate and deactivate.
   ```
    * Run `eval "$(mamba shell hook --shell bash)"`, and then activate the bigfish environment
4. Start the Notebook (jupyter lab) with `jupyter lab --ip <Aries IP> --port 8888`
5. Navigate to `<Aries IP>:8888` to access notebook. 

Method-specific steps/information is in the notebook.

#### Note on File Naming:
BigFISH offers the ability to build a 'recipe' when reading files from multiple different channels / FOVs. It is expected that, defined by the recipe, files will only contain unique identifiers with the exception of one 'optional' part of the name that is shared amongst the files. For example, if your files are named like: \
`Experiment_FOV1_DAPI.tif`, `Experiment_FOV2_DAPI.tif` ... `Experiment_FOV10_DAPI.tif` \
the recipe would look like:
```
recipe = {
    "fov": ["FOV1,"FOV2",...,"FOV10"],
    "c": "DAPI",
    "opt": "Experiment",
    "ext": "tif",
    "pattern": "opt_fov_c.ext"}
```
BigFISH will search for files with the pattern specified by the `pattern` key. 
