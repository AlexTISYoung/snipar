#!/usr/bin/env bash

set -e

source ${snipar_env}/activate

ibd.py \
	--bgen ${BGEN_PATH}/ukb_hap_chr@_v2_filtered \
	--chr_range 1-22 \
	--king ${KING_OUTPUT_PATH}/ukb_all_ancestry.kin0 \
	--agesex ${PHEN_PATH}/agesex.txt \
	--out ${IBD_PATH}/all_ancestries/chr@ \
	--ibdmatrix \
	--ld_out \
	--threads 15
