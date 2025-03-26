import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import cooler
import bioframe
import cooltools
import bbi
from scipy.stats import linregress
import cooltools.lib.plotting
import argparse
import sys
import datetime

parser = argparse.ArgumentParser()

parser.add_argument('cool', type=str, help="File path to cool/mcool file")
parser.add_argument('genome', type=str, help="Genome build used to fetch chromosome arm sizes ('hg19' or 'hg38')")
parser.add_argument('peaks', type=str, help="File path to bed (or bedpe) style file with ChIP peaks to pile 3C contacts around")
parser.add_argument('filetype', type=str, help="'bed' or 'bedpe', depending on what type of file the peaks are represented in")
parser.add_argument('map_out', type=str, help="File path to output pileup heatmap")
parser.add_argument('--filter', type=str, help="Use only top quartile of peaks by score and enrichment. If on (T), bigWig must be supplied via --bigwig (default = F)", default = "F", choices = ["T", "F"])
parser.add_argument('--bigwig', type=str, help="File path to bigWig file corresponding to ChIP input. Used to extract enrichment information. (default = None)", default = None)
parser.add_argument('--chip_flank', type=int, help="Extend input ChIP peaks by this amount, in bp (default = 250)", default=250)
parser.add_argument('--ctrl', type=str, help=".cool/.mcool file to use as background when calculating bin enrichment")
parser.add_argument('--map_flank', type=int, help="Area to be shown in pileup plot, in bp from the center (default = 100000)", default=100000)
parser.add_argument('--res', type=int, help="Resolution to use when piling up (default = 10000)", default=10000)
parser.add_argument('--title', type=str, help="Title of output plieup plot (default = 'Pileup')", default = "Pileup")
parser.add_argument('--nproc', type=int, help="Number of cores to use (default = 10)", default=10)

args = parser.parse_args()

# If input is bedpe, split into two dataframes

print("{}: Fetching chromosome arm sizes . . .".format(datetime.datetime.now()))
## Get chromsizes: NOTE this is hard coded to match with distiller
hg19_chromsizes = pd.read_table("/mnt/data0/tfriedr/Genomes/hg38.chrom_red_test.sizes", header = None)
hg19_chromsizes.columns = ["chrom","length"]
hg19_cens = bioframe.fetch_centromeres(args.genome)
hg19_arms = bioframe.make_chromarms(hg19_chromsizes, hg19_cens)

print("{}: Loading interaction matrix . . .".format(datetime.datetime.now()))
#Load cool files
if ".cool" in args.cool:
    cooler_obj = cooler.Cooler(args.cool)
elif ".mcool" in args.cool:
    cooler_obj = cooler.Cooler("{}::/resolutions/{}".format(args.cool, args.res))

print("{}: Calculating expected interactions . . .".format(datetime.datetime.now()))
#Calculate expected interactions within chromosome arms
#Check if using different file for backgrond. If yes, load it into cooler obj
if args.ctrl is not None:
    if ".mcool" in args.ctrl:
        ctrl_cool = cooler.Cooler("{}::/{}".format(args.ctrl, args.res))
    elif ".cool" in args.ctrl:
        ctrl_cool = cooler.Cooler("{}::/{}".format(args.ctrl, args.res))
    expected =  cooltools.expected_cis(ctrl_cool,
                                   view_df = hg19_arms,
                                   nproc = args.nproc,
                                   chunksize = 1_000_000)
else:
    expected = cooltools.expected_cis(cooler_obj,
                                   view_df = hg19_arms,
                                   nproc = args.nproc,
                                   chunksize = 1_000_000)


print("{}: Loading ChIP peaks . . .".format(datetime.datetime.now()))
#Load anchor features (chip peaks)
if args.filetype == 'bed':
    chip = bioframe.read_table(args.peaks, schema='bed') #
    chip['mid'] = (chip.end + chip.start)//2 # Use midpoint of peak region as feature
elif args.filetype == 'bedpe':
    paired_sites = bioframe.read_table(args.peaks, header = 0)
    #paired_sites = paired_sites.set_axis(['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'], axis=1)
    paired_sites[['start1','end1','start2','end2']] = paired_sites[['start1','end1','start2','end2']].apply(pd.to_numeric) 
    print(paired_sites.dtypes)
############################################ Optional filtering step #####################################
if args.filter == "T":
    print("{}: Fetching ChIP peak enrichment . . .".format(datetime.datetime.now()))
    chip_signal = bbi.stackup( #query signal from bigWig file 
        args.bigwig,
        chip.chrom,
        chip.mid - args.chip_flank,
        chip.mid + args.chip_flank,
        bins=1)
    chip['FC_score'] = chip_signal

    print("{}: Fetching highest scoring peaks . . .".format(datetime.datetime.now()))
    #Divide score and FC into quartiles
    chip['quartile_score'] = pd.qcut(chip['score'], 4, labels = False) + 1
    chip['quartile_FC_score'] = pd.qcut(chip['FC_score'], 4, labels = False) + 1

    #Assign each peak to a category based on quartiles
    chip['peaks_importance'] = chip.apply(
        lambda x: 'Top by both scores' if x.quartile_score==4 and x.quartile_FC_score==4 else
        'Top by Motif score' if x.quartile_score==4 else
        'Top by FC score' if x.quartile_FC_score==4 else 'Ordinary peaks', axis = 1)

    #Subset chip peaks for 'top by both scores' peaks
    chip = chip[chip['peaks_importance']=='Top by both scores']\
            .sort_values('FC_score', ascending=False)\
            .reset_index(drop=True)
    #Some chip sites are too close together, will effect analysis
    #Collapes sites falling into same size genomic bins as cool res (10,000)
    chip = bioframe.cluster(chip, min_dist=args.res)\
        .drop_duplicates('cluster')\
        .reset_index(drop=True)

elif args.filter == "F":
    print("{}: Skipping ChIP peak filtering!".format(datetime.datetime.now()))
############################################################################################################

#Pair ChIP anchors and get midpoints - skip if bedpe was input
if args.filetype == 'bed':
    paired_sites = bioframe.pair_by_distance(chip, min_sep=200_000, max_sep=1_000_000, suffixes=('1', '2'))
    paired_sites.loc[:, 'mid1'] = (paired_sites['start1'] + paired_sites['end1'])//2
    paired_sites.loc[:, 'mid2'] = (paired_sites['start2'] + paired_sites['end2'])//2

print("{}: Pileup and average snips . . .".format(datetime.datetime.now()))
#Make snips:
stack = cooltools.pileup(cooler_obj, paired_sites, view_df=hg19_arms, expected_df=expected, flank=args.map_flank)
#Aggregate snips:
mtx = np.nanmean(stack, axis=2)

print("{}: Generating visualization at resolution {} . . .".format(datetime.datetime.now(), args.res))
#Visualize matrix:
plt.imshow(
    np.log2(mtx),
    vmax = 1,
    vmin = -1,
    cmap='BluYl')

plt.colorbar(label = 'log2 mean obs/exp')
ticks_pixels = np.linspace(0, args.map_flank*2//args.res,5)
ticks_kbp = ((ticks_pixels-ticks_pixels[-1]/2)*args.res//1000).astype(int)
plt.xticks(ticks_pixels, ticks_kbp)
plt.yticks(ticks_pixels, ticks_kbp)
plt.xlabel('relative position, kbp')
plt.ylabel('relative position, kbp')
plt.title(args.title)

plt.savefig(args.map_out)

print("{}: DONE!".format(datetime.datetime.now()))
