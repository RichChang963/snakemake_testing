# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2022-2023 The PyPSA-Eur and Agora Energiewende Authors

# SPDX-License-Identifier: MIT

# coding: utf-8

import os
from os import PathLike
from pathlib import Path
from typing import Union
import logging

FilePath = Union[str, "PathLike[str]"]


def configure_logging(snakemake, skip_handlers=False):
    """Configure the basic behaviour for the logging module.
    Note: Must only be called once from the __main__ section of a script.
    The setup includes printing log messages to STDERR and to a log file defined
    by either (in priority order): snakemake.log.python, snakemake.log[0] or "logs/{rulename}.log".
    Additional keywords from logging.basicConfig are accepted via the snakemake configuration
    file under snakemake.config.logging.

    Parameters
    ----------
    snakemake : Snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : boolean, optional
        Do (not) skip the default handlers created for redirecting output to STDERR and file.
        Default by False
    """
    kwargs = snakemake.config.get("logging", dict())
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath(
            "..", "logs", f"{snakemake.rule}.log"
        )
        logfile = snakemake.log.get(
            "python", snakemake.log[0] if snakemake.log else fallback_path
        )
        kwargs.update(
            {
                "handlers": [
                    # Prefer the 'python' log, otherwise take the first log for each
                    # Snakemake rule
                    logging.FileHandler(logfile),
                    logging.StreamHandler(),
                ]
            }
        )
    logging.basicConfig(**kwargs)


def mock_snakemake(rulename, **wildcards):
    """This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.
    If a rule has wildcards, you have to specify them in **wildcards.

    Parameters
    ----------
    rulename : str
        Name of the rule for which the snakemake object should be generated
    **wildcards:
        Keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.

    Returns
    -------
    Snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    """
    import os
    import snakemake as sm
    from packaging.version import Version, parse
    from pypsa.descriptors import Dict
    from snakemake.script import Snakemake

    script_dir = Path(__file__).parent.resolve()
    assert (
        Path.cwd().resolve() == script_dir
    ), f"mock_snakemake has to be run from the repository scripts directory {script_dir}"
    os.chdir(script_dir.parent)
    for p in sm.SNAKEFILE_CHOICES:
        if os.path.exists(p):
            snakefile = p
            break
    kwargs = dict(rerun_triggers=[]) if parse(sm.__version__) > Version("7.7.0") else {}
    workflow = sm.Workflow(snakefile, overwrite_configfiles=[], **kwargs)
    workflow.include(snakefile)
    workflow.global_resources = {}
    rule = workflow.get_rule(rulename)
    dag = sm.dag.DAG(workflow, rules=[rule])
    wc = Dict(wildcards)
    job = sm.jobs.Job(rule, dag, wc)

    def make_accessable(*ios):
        for io in ios:
            for i in range(len(io)):
                io[i] = os.path.abspath(io[i])

    make_accessable(job.input, job.output, job.log)
    snakemake = Snakemake(
        job.input,
        job.output,
        job.params,
        job.wildcards,
        job.threads,
        job.resources,
        job.log,
        job.dag.workflow.config,
        job.rule.name,
        None,
    )
    # create log and output dir if not existent
    for path in list(snakemake.log) + list(snakemake.output):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    os.chdir(script_dir)

    return snakemake