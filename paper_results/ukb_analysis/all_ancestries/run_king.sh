#!/usr/bin/env bash

set -e

DIR=/disk/genetics/ukb/jguan/preprocessed
PLINK_DIR=${DIR}/plink
merged=${PLINK_DIR}/ukb_merged_v2

# king
${KING_PATH}/king -b ${BED_PATH}/all_ancestries/genome_filtered.bed \
    --related \
    --degree 1 \
    --prefix ${KING_OUTPUT_PATH}/ukb_all_ancestry \
	--cpus 5
