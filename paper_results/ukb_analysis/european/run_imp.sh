#!/usr/bin/env bash


source ${snipar_env}/activate

impute.py \
	--ibd ${IBD_PATH}/european/v2chr_@.ibd \
	--bgen ${BGEN_PATH}/european/ukb_hap_chr#_v2_filtered \
	--king ${KING_OUTPUT_PATH}/ukb_v2.kin0 \
	--agesex ${PHEN_PATH}/agesex.txt \
	--out ${IMP_PATH}/european/chr_@ \
	--chr_range 1-22 \
	--threads 5
