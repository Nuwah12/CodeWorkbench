import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
from matplotlib.ticker import EngFormatter
from matplotlib.lines import Line2D
import argparse
import sys
import cooler
from cooltools.lib.numutils import fill_diag
import datetime
import pandas as pd
import numpy as np
import bioframe
import os.path

parser = argparse.ArgumentParser()
parser.add_argument('--view', '-v', type=str, help="Region to view on map, ex. chr1,150000000,160000000", required = True)
parser.add_argument('--dots', '-d', type=str, help="Path to the file containing the called dots", required = True, nargs="+")
parser.add_argument('--out', '-o', type=str, help="Output DIRECTORY to store matrices w/ dots", required = True)
parser.add_argument('--cool', '-c', type=str, help="FILE path to the .mcool/.cool file", required = True)
parser.add_argument('--res', '-r', type=int, help="Resolution to plot matrix in", required = True)

args = parser.parse_args()

#Plotting function
def rectangles_around_dots(dots_df, region, loc="upper", lw=1, ec="blue", fc="none"):
    """
    yield a series of rectangles around called dots in a given region
    """
    # select dots from the region:
    df_reg = bioframe.select(
        bioframe.select(dots_df, region, cols=("chrom1","start1","end1")),
        region,
        cols=("chrom2","start2","end2"),
    )
    rectangle_kwargs = dict(lw=lw, ec=ec, fc=fc)
    # draw rectangular "boxes" around pixels called as dots in the "region":
    for s1, s2, e1, e2 in df_reg[["start1", "start2", "end1", "end2"]].itertuples(index=False):
        width1 = e1 - s1
        width2 = e2 - s2
        if loc == "upper":
            yield patches.Rectangle((s2, s1), width1, width2, **rectangle_kwargs)
        elif loc == "lower":
            yield patches.Rectangle((s1, s2), width1, width2, **rectangle_kwargs)
        else:
            raise ValueError("loc has to be uppper or lower")

if ".mcool" in args.cool:
    clr = cooler.Cooler("{}::/resolutions/{}".format(args.cool, args.res))
elif ".cool" in args.cool:
    clr = cooler.Cooler(args.cool)
else:
    sys.exit("Incorrect cooler format . . . try again")

for i in args.dots:
    print("{}: Plotting dots on interaction heatmap, dot location = {} . . .".format(datetime.datetime.now(), i))
    dots = pd.read_csv(i, sep = "\t")
    print("{} dots read".format(dots.shape[0]))
    #Define region to view dots
    view = args.view.split(",")
    start = int(view[1])
    end = int(view[2])
    region = (view[0], start, end)

    matshow_kwargs = dict(
        cmap='YlOrBr',
        norm=LogNorm(vmax=0.5),
        extent=(start, end, end, start))

    colorbar_kwargs = dict(fraction = 0.046, label = 'Correct Frequencies')

    #Heatmap
    region_matrix = clr.matrix(balance = True).fetch(region)
    for diag in [-1,0,1]:
        region_matrix = fill_diag(region_matrix, np.nan, i = diag)

    f, ax = plt.subplots(figsize = (7, 7))
    im = ax.matshow(region_matrix, **matshow_kwargs)
    
    plt.colorbar(im, ax=ax, **colorbar_kwargs)

    for box in rectangles_around_dots(dots, region, lw=1.5):
       ax.add_patch(box)
    
    #output file name
    outf = "{}/{}_dotsOnMatrix_{}_res{}.png".format(args.out, os.path.basename(i), args.view, args.res)
    plt.savefig(outf)
