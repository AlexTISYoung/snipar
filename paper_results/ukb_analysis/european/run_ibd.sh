#!/usr/bin/env bash

ENV=/homes/nber/jguan/SNIPar_tidy/bin

source $ENV/activate

ibd.py \
	--bgen ${BGEN_PATH}/european/ukb_hap_chr@_v2_filtered \
	--king ${KING_OUTPUT_PATH}/ukb_v2.kin0 \
	--agesex ${PHEN_PATH}/agesex.txt \
	--out ${IBD_PATH}/european/v2 \
	--ibdmatrix \
	--ld_out \
	--threads 20
