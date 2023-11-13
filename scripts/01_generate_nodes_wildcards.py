# coding: utf-8

import os
# extra pre-setting for geopandas to exclude unclosed rings in polygons
os.environ["USE_PYGEOS"] = "0"
os.environ["OGR_GEOMETRY_ACCEPT_UNCLOSED_RING"] = "NO"
import geopandas as gpd
from _helpers import FilePath, configure_logging


def clean_node(nuts: FilePath):
    """Save only region name and node_ids into pd.DataFrame

    Parameters
    ----------
    nuts : FilePath
        Input geojson file of nuts level
    excl_iso_list: List, optional
        Exclude specific countries. Default by None
    country: string, optional
        Filter specific country. Default by None
    """
    nodes = gpd.read_file(nuts)

    nodes = nodes[["id", "na"]].rename(columns={"na": "node"})

    nodes.sort_values(by="id", ascending=True).sort_index().to_json(snakemake.output.node_wildcards)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("generate_wildcards")
    configure_logging(snakemake)
    
    clean_node(nuts=snakemake.input.nodes)


