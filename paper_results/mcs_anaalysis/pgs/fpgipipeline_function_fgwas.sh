#!/usr/bin/env bash

set -e
source ${snipar_env}/bin/activate
function withinfam_pred(){

    WTFILE=$1
    EFFECT=$2
    PHENONAME=$3
    OUTSUFFIX=$4
    BINARY=$5
    METHOD=$6
    ANCESTRY=$7

    OUTPATH="prs_analysis/mcs/nofilter/${PHENONAME}/${METHOD}/${ANCESTRY}"
    RAWPATH="mcs_data"
    COVAR="mcs_data/phen/PCs.txt"

    echo $OUTPATH/${EFFECT}${OUTSUFFIX}
    echo $WTFILE
    mkdir -p $OUTPATH


    pheno_out="$RAWPATH/phen/${PHENONAME}/${ANCESTRY}"


    bedfilepath="mcs_data/bgen/tmp/chr@.dose"
    impfilespath="mcs_data/imputed_parents/chr@"

    ## get proband and parental pgis using snipar        
    pgs.py \
        $OUTPATH/${EFFECT}${OUTSUFFIX} \
        --bed $bedfilepath \
        --imp $impfilespath \
        --beta_col "ldpred_beta" \
        --SNP "sid" \
        --A1 "nt1" \
        --A2 "nt2" \
        --weights prs_analysis/prscs_weights/nofilter/${METHOD}/${PHENONAME}/${PHENONAME}_${EFFECT}_fpgs_formatted.txt \
        --scale_pgs | tee $OUTPATH/${EFFECT}${OUTSUFFIX}.log 

    scoresout="$OUTPATH/${EFFECT}${OUTSUFFIX}.pgs.txt"
    fpgs_out="prs_analysis/mcs/nofilter/fpgs/${PHENONAME}/${METHOD}/${ANCESTRY}"
    mkdir -p $fpgs_out

    ## run fPGI regression
    echo "Run fPGI regression..."
    echo "Output to ${fpgs_out}/${EFFECT}${OUTSUFFIX}_PCadjusted"
    pgs.py ${fpgs_out}/${EFFECT}${OUTSUFFIX}_PCadjusted \
        --pgs ${scoresout} \
        --covar $COVAR \
        --phenofile ${pheno_out}/pheno_withPCs.pheno \
        --scale_pgs \
        --scale_phen | tee "prs_analysis/mcs/nofilter/fpgs/logs/${PHENONAME}_${EFFECT}${OUTSUFFIX}_PCadjusted_${ANCESTRY}_full.reg.log"

}

function main(){

    PHENONAME=$1
    OUTSUFFIX=$2
    BINARY=$3
    METHOD=$4
    ANCESTRY=$5
    POPULATION=$6

    RAWPATH="/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed"
    direct_weights="/disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/prscs_weights/nofilter/${METHOD}/${PHENONAME}/direct/weights/meta_weights.snpRes"
    pheno_out="$RAWPATH/phen/${PHENONAME}/${ANCESTRY}"

    # main prediction -- direct effect pgi
    withinfam_pred $direct_weights \
        "direct" "$PHENONAME" \
        "$OUTSUFFIX" "$BINARY" "$METHOD" "$ANCESTRY"

    if [ "$POPULATION" == "dir_pop" ]; then
        # population effect pgi
        population_weights="/disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/prscs_weights/nofilter/${METHOD}/${PHENONAME}/population/weights/meta_weights.snpRes"
        withinfam_pred $population_weights \
            "population" "$PHENONAME" \
            "$OUTSUFFIX" "$BINARY" "$METHOD" "$ANCESTRY"
    fi
    
}