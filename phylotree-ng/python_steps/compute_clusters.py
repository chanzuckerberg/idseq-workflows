import argparse
import json
import math

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import cluster

matplotlib.use('Agg')

class TooDivergentError(Exception):
    pass

def main(ska_distances: str, trim_height: float):
    # we have observed some strange parsing behavior of this file, ensure it works with end to end testing
    # we may need to use a regex separator
    df = pd.read_csv(ska_distances, sep='\t')
    df = pd.concat([
        pd.DataFrame(dict(zip(df.columns, [df.iloc[0][0], df.iloc[0][0]] + [0] * 6)), index=[0]),
        df,
        pd.DataFrame(dict(zip(df.columns, [df.iloc[-1][0], df.iloc[-1][0]] + [0] * 6)), index=[0]),
    ])
    df.reset_index(drop=True, inplace=True)

    # long dataframe to wide
    df2 = df.pivot_table(index=['Sample 1'], columns='Sample 2', values='Mash-like distance')

    # fill lower triangle of the matrix (currently NA)
    npdf = df2.to_numpy()
    i_lower = np.tril_indices(len(df2.index), -1)
    npdf[i_lower] = npdf.T[i_lower]
    df3 = pd.DataFrame(npdf)
    df3.columns = df2.columns
    df3.index = df2.index
    df3.fillna(0, inplace=True)

    Z = cluster.hierarchy.linkage(1-df3, method='complete')
    cutree = cluster.hierarchy.cut_tree(Z, height=float(trim_height))
    ordered_clusterids = [i[0] for i in cutree]
    cluster_assignments = dict(zip(df3.index, ordered_clusterids))
    n_clusters = len(set(ordered_clusterids))

    cluster_sets = dict(zip(list(set(ordered_clusterids)), [[] for _ in range(n_clusters)]))

    for i in cluster_assignments.keys():
        cluster_sets[cluster_assignments[i]].append(i)

    # write cluster contents to files for future processing
    n_clusters_out = 0
    for c, s in cluster_sets.items():
        if(len(s)) > 2: # only output files where there are > 2 samples
            n_clusters_out += 1
            filenames = '\n'.join(cluster_sets[c])
            with open("./cluster_files/cluster_" + str(c), "w") as text_file:
                text_file.write(filenames)

    if n_clusters_out > 1:
        exit(json.dumps(dict(
            wdl_error_message=True,
            error="TooDivergentError",
            cause="Sequences are too divergent to create a single phylo tree",
        )))

    color_list = sns.color_palette("Dark2", 8)
    long_color_list = color_list*math.ceil(len(set(ordered_clusterids))/len(color_list))
    col_colors = [long_color_list[i] for i in ordered_clusterids]
    sns.clustermap(df3, cmap='coolwarm_r', vmin = 0, vmax = 0.15, col_linkage = Z, col_colors = col_colors, figsize=(15,15))
    plt.savefig('clustermap.png', bbox_inches='tight')
    plt.savefig('clustermap.svg', bbox_inches='tight')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ska-distances")
    parser.add_argument("--cut-height", type=float)
    args = parser.parse_args()
    main(args.ska_distances, args.cut_height)