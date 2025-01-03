## 3C Assay preprocessing pipeline
This workflow uses [juicer](https://github.com/aidenlab/juicer) to align, filter, and normalize the data, followed by [hic2cool](https://github.com/4dn-dcic/hic2cool) for conversion into a `.cool` file.
### Usage:
This workflow is written in `WDL`, and is configured to use `Cromwell` to run. Instructions on running a Cromwell workflow can be found in the README of this repo: [Dockerized Workflows](https://github.com/faryabiLab/dockerize-workflows).

Rundown of the workflow parameters found in `juicer_cooler.json`:
* `sampleList` - Sample sheet containing the prefixes of all samples to be processed.
* `ChromSizes` - A standard BED-style Chromosome Sizes file.
* `TopDir` - Directory that the workflow will be run from. Fastq files should be stored here in a directory called `fastq`.
* `GenomeAssembly` - Path to a `.fasta` genome sequence file - **the directory should also contain all BWA Index files**
* `Enzyme` - Path to a file containing the restriction sites of the enzyme used in this assay ('none' if no res. enzyme used, 'HindIII' and 'MboI' are included by default.)

### Output:
Output location can be specified via the `options.json` file under key `final_workflow_outputs_dir`, and final log locations can be similarly specified under key `final_call_logs_dir`. 
