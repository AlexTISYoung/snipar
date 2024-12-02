set -e


${KING_PATH}/king -b ${BED_PATH}/ukb_merged_v2.bed \
    --ibdseg \
	--degree 3 \
	--cpus 20 \
    --prefix ${KING_IBDSEG_PATH}/ukb_v2_3rd
