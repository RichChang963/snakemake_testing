# coding: utf-8

import os
import pandas as pd
from _helpers import FilePath, configure_logging


def separate_nodes(geo_path: FilePath):
    if not os.path.exists(f"{geo_path}/pickles"):
        os.makedirs(f"{geo_path}/pickles")

    node_df = pd.read_json(geo_path + "/node_wildcards.json")
    for node in node_df.id:
        small_df = node_df.loc[node_df.id == node]
        small_df = small_df[["id", "node"]]
        small_df.to_pickle(
            f"{geo_path}/pickles/node_{node}.pkl"
            )    

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("transform_nodewise_timeseries")
    configure_logging(snakemake)

    geo_path = snakemake.input.geo_path
   
    separate_nodes(geo_path=geo_path)


