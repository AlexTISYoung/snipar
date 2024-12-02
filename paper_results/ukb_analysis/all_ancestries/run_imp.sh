#!/usr/bin/env bash

set -e

source ${snipar_env}/bin/activate

impute.py \
	--ibd ${IBD_PATH}/all_ancestries/chr@.ibd \
	--bgen ${BGEN_PATH}/all_ancestries/ukb_hap_chr@_v2_filtered \
	--king ${KING_OUTPUT_PATH}/ukb_all_ancestry.kin0 \
	--agesex ${PHEN_PATH}/agesex.txt \
	--out ${IMP_PATH}/all_ancestries/chr@ \
	--chr_range 1-22 \
	--threads 10
