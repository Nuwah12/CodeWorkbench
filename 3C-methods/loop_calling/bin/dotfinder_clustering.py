import cooler
import pandas as pd
import numpy as np
from sklearn.cluster import Birch
import argparse
import logging
import os.path

# Implement cooltools dotfinder clustering algorithm for the final unioned loopset in the tri-method 

parser = argparse.ArgumentParser()

parser.add_argument("--dots", type=str, required=True)
parser.add_argument("--radius", type=str, required=True)

args = parser.parse_args()

def clust_2D_pixels(
    pixels_df,
    threshold_cluster=2,
    bin1_id_name="bin1_id",
    bin2_id_name="bin2_id",
    clust_label_name="c_label",
    clust_size_name="c_size"):

    """
    Parameters
    ----------
    pixels_df : dataframe with pixel coordinates that must have AT LEAST 2 COLUMNS that are 'bin1_id' (row) and 'bin2_ID' (column)
    threshold_cluster : clustering radius for Birch clustering **derived from 40kb radius** (eg. 2 = 40kb/2 = 20kb clustering radius)
    bin1_id_name : Name of 'bin1_id' column, by default 'bin1_id'
    bin2_id_name : Name of 'bin_2id' column, by default 'bin2_id'
    clust_label_name : Name of clusters to label (default = "c_label")
    clust_size_name : Name of the cluster of pixels size (default = "c_size")
    """
    
    pixels = pixels_df[[bin1_id_name, bin2_id_name]].values.astype(np.float64)
    pixel_idxs = pixels_df.index

    # Perform Birch clustering of pixels
    brc = Birch(
            n_clusters=None,
            threshold=threshold_cluster,
            compute_labels=True)

    brc.fit(pixels)
    
    clustered_labels = brc.labels_

    #centrid coordinates
    clustered_centroids = brc.subcluster_centers_

    #count unique labels and get their continuous indices 
    uniq_labels, inverse_idx, uniq_counts = np.unique(clustered_labels,
                                                     return_inverse=True,
                                                     return_counts=True)
    #cluster sizes of matching labels
    cluster_sizes = uniq_counts[inverse_idx]

    #take centroids cooresponding to labels
    centroids_per_pixel = np.take(clustered_centroids, clustered_labels, axis=0)

    logging.info(
        f"detected {uniq_counts.size} clusters of {uniq_counts.mean():.2f}+/-{uniq_counts.std():.2f} size"
    )

    #Make output dataframe
    centroids_n_labels_df = pd.DataFrame(
                                        centroids_per_pixel,
                                        index=pixel_idxs,
                                        columns = ["c" + bin1_id_name, "c" + bin2_id_name])

    #add labels per pixel
    centroids_n_labels_df[clust_label_name] = clustered_labels.astype(np.int64)

    #add cluster sizes
    centroids_n_labels_df[clust_size_name] = cluster_sizes.astype(np.int64)

    return centroids_n_labels_df

    
#Main function
def clustering_step(
    scored_df,
    dots_clustering_radius,
    assigned_regions_name="region",
    obs_raw_name="count",
):

    """
    Parameters
    ----------
    scored_df : dataframe with enriched pixels that are ready to be clustered and are annotated with their genomic coordinates
    dots_clustering_radius : Birch clustering radius threshold
    assigned_regions_name : Name of the column in scored_df to use for for grouping pixels before clstering - NONE = FULL-CHROMOSOME CLUSTERING (DO THIS)
    obs_raw_name : name of the column with raw observed pixel counts

    """
    #Read in scored_df to pandas dataframe (since we will be using this out-of-the-box)
    scored_df = pd.read_csv(scored_df, sep = "\t")

    #Check if these column names exist in dataframe - if not, raise error
    if not {"chrom1", "chrom2", "start1", "start2", obs_raw_name}.issubset(scored_df):
        raise ValueError("Scored pixels provided for clustering are not properly annotated.")

    #Assigning a region to cluster independently 
    scored_df = scored_df.copy()
    if not assigned_regions_name in scored_df.columns:
        logging.warning("No regions assigned to the scored pixels before clustering, using chromosomes")
        scored_df[assigned_regions_name] = np.where(scored_df["chrom1"] == scored_df["chrom2"], scored_df["chrom1"], np.nan)
    
    #Group dataframe by the group name assigned above
    pixel_clust_list = []
    scored_pixels_by_region = scored_df.groupby(assigned_regions_name, observed = True)

    #For each region, cluster the dots in that region
    for region, _df, in scored_pixels_by_region:
        logging.info("Clustering enriched pixels in region: {}".format(region))

        #Using genomic coordinates for clustering, not bin_id
        pixel_clust = clust_2D_pixels(_df, 
                                     threshold_cluster=int(dots_clustering_radius),
                                     bin1_id_name="start1",
                                     bin2_id_name="start2")
        pixel_clust_list.append(pixel_clust)
    logging.info("Finished clustering!")

    #concatenate clustering results
    #first, check if no clusters were found and return an empty dataframe
    if not pixel_clust_list:
        logging.warning("No clusters found fo any region, output will be empty")
        empty_output = pd.DataFrame(
                                    [],
                                    columns=list(scored_df.columns)
                                    +[
                                    assigned_regions_name + "1",
                                    assigned_regions_name + "2",
                                    "c_label",
                                    "c_size",
                                    "cstart1",
                                    "cstart2",
                                    ],)
        return empty_output

    else:
        #Concatenate clustering results from different regions
        pixel_clust_df = pd.concat(
                                  pixel_clust_list, 
                                  ignore_index=False)

    #Merge pixel_clust_df with scored_df
    df = pd.merge(
                 scored_df,
                 pixel_clust_df,
                 how="left",
                 left_index=True,
                 right_index=True)

    df[assigned_regions_name+"1"] = df[assigned_regions_name].astype(str)
    df[assigned_regions_name+"2"] = df[assigned_regions_name].astype(str)

    #Report only the centroids with highest observed in each cluster
    chrom_clust_group = df.groupby(
                                  [assigned_regions_name+"1", assigned_regions_name+"2", "c_label"],
                                  observed=True)

    centroids = df.loc[chrom_clust_group[obs_raw_name].idxmax()]

    return centroids


###main
out = clustering_step(args.dots, args.radius)

out.to_csv("{}_clustered.dots".format(args.dots), sep = "\t")
























