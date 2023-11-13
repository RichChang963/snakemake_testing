# coding: utf-8

from os.path import normpath, exists
from shutil import copyfile, move
import json

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: "config.yaml"

wildcard_constraints:
    nodes="[-+a-zA-Z0-9\.\s]*"

RDIR = config["results_dir"] + config["output_folder"]
CIDIR = config["input_dir"] + "/" + config["country_input_dir"]
INDIR = config["input_dir"]

def open_wildcards():
    try:
        with open(f"{INDIR}/node_wildcards.json", "r") as f:
            node_config = json.load(f)
        return node_config["id"].values()
    except:
        print(f"No json file available in {INDIR}.")


rule generate_wildcards:
    input:
        nodes = f"{INDIR}/four_nodes.geojson",
    output:
        node_wildcards = f"{INDIR}/node_wildcards.json",
    log: "logs/generate_wildcards.log"
    benchmark: "benchmarks/generate_wildcards"
    threads: 2
    script: "scripts/01_generate_nodes_wildcards.py"


rule separate_nodes:
    input:
        geo_path = f"{INDIR}",
    log: "logs/02_separate_nodes.log"
    benchmark: "benchmarks/02_separate_nodes"
    threads: 2
    script: "scripts/02_separate_nodes.py"


rule test_all_wildcards:
    input:
        expand(RDIR + "networks/lcoh_{nodes}.nc",
        nodes=open_wildcards())
        

rule test_wildcard:
    input:
        nodes = CIDIR + "pickles/node_{nodes}.pkl",
    output:
        network = RDIR + "networks/lcoh_{nodes}.nc",
    log: "logs/solve_network_lcoh_{nodes}.log"
    benchmark: "benchmarks/solve_network_lcoh_{nodes}"
    threads: 2
    script: "scripts/03_test_wildcard.py"
