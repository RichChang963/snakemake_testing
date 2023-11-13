# coding: utf-8

import logging

import pandas as pd
import numpy as np
import pypsa
import warnings
from _helpers import FilePath, configure_logging

logger = logging.getLogger(__name__)

warnings.simplefilter(action="ignore", category=FutureWarning)



if __name__ == "__main__":                
    if "snakemake" not in globals():
        from _helpers import mock_snakemake
        
        snakemake = mock_snakemake("solve_network", nodes="DEU_131")
    configure_logging(snakemake)

    node_id = snakemake.wildcards.nodes
    print(node_id)
    df = pd.read_pickle(snakemake.input.nodes)
    df.to_csv(snakemake.output.network)
