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
    COVAR="mcs_data/phen/PCs_sas.txt"

    
    # mkdir -p $RAWPATH/phen/${PHENONAME}
    echo $OUTPATH/pop_pgs/${EFFECT}${OUTSUFFIX}
    echo $WTFILE
    mkdir -p $OUTPATH/pop_pgs
    mkdir -p prs_analysis/mcs/nofilter/pop_pgs/logs

    pheno_out="$RAWPATH/phen/${PHENONAME}/${ANCESTRY}"


    bedfilepath="mcs_data/bgen/SAS/tmp/chr@.dose"


    # get proband and parental pgis using snipar        
    python pgs_population_eff.py \
        $OUTPATH/pop_pgs/${EFFECT}${OUTSUFFIX} \
        --bed $bedfilepath \
        --beta_col "ldpred_beta" \
        --SNP "sid" \
        --A1 "nt1" \
        --A2 "nt2" \
        --weights prs_analysis/prscs_weights/nofilter/${METHOD}/${PHENONAME}/${PHENONAME}_${EFFECT}_fpgs_formatted.txt \
        --scale_pgs


   
}

function main(){

    PHENONAME=$1
    OUTSUFFIX=$2
    BINARY=$3
    METHOD=$4
    ANCESTRY=$5
    POPULATION=$6


    RAWPATH="mcs_data"
    direct_weights="prs_analysis/prscs_weights/nofilter/${METHOD}/${PHENONAME}/direct/weights/meta_weights.snpRes"
    pheno_out="$RAWPATH/phen/${PHENONAME}/${ANCESTRY}"

    withinfam_pred $direct_weights \
        "direct" "$PHENONAME" \
        "$OUTSUFFIX" "$BINARY" "$METHOD" "$ANCESTRY"

    if [ "$POPULATION" == "dir_pop" ]; then

        population_weights="prs_analysis/prscs_weights/nofilter/${METHOD}/${PHENONAME}/population/weights/meta_weights.snpRes"
        withinfam_pred $population_weights \
            "population" "$PHENONAME" \
            "$OUTSUFFIX" "$BINARY" "$METHOD" "$ANCESTRY"
    fi
    
}