set -e


# king
${KING_PATH}/king -b ${BED_PATH}/ukb_merged_v2.bed \
    --related \
    --degree 1 \
    --prefix ${KING_OUTPUT_PATH}/ukb_v2 \
	--cpus 5
