#!/bin/bash

snakemake --cluster qsub --jobscript sge_jobscript.sh --jobs 60
