#!/usr/bin/env bash


set -e
source ~/snipar_release/bin/activate
function withinfam_pred(){

    WTFILE=$1
    EFFECT=$2
    PHENONAME=$3
    OUTSUFFIX=$4
    BINARY=$5
    METHOD=$6
    ANCESTRY=$7

    OUTPATH="/disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/mcs/nofilter/${PHENONAME}/${METHOD}/${ANCESTRY}"
    RAWPATH="/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed"
    COVAR="/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/phen/PCs_sas.txt"
    # PHENOFILE="/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/phen/phenotypes_eur.txt"
    
    # mkdir -p $RAWPATH/phen/${PHENONAME}
    echo $OUTPATH/pop_pgs/${EFFECT}${OUTSUFFIX}
    echo $WTFILE
    mkdir -p $OUTPATH/pop_pgs
    mkdir -p /disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/mcs/nofilter/pop_pgs/logs

    pheno_out="$RAWPATH/phen/${PHENONAME}/${ANCESTRY}"


    bedfilepath="/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/bgen/SAS/tmp/chr@.dose"
    impfilespath="/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/imputed_parents/chr@"

    # get proband and parental pgis using snipar        
    python pgs_population_eff.py \
        $OUTPATH/pop_pgs/${EFFECT}${OUTSUFFIX} \
        --bed $bedfilepath \
        --beta_col "ldpred_beta" \
        --SNP "sid" \
        --A1 "nt1" \
        --A2 "nt2" \
        --weights /disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/prscs_weights/nofilter/${METHOD}/${PHENONAME}/${PHENONAME}_${EFFECT}_fpgs_formatted.txt \
        --scale_pgs


    scoresout="$OUTPATH/pop_pgs/${EFFECT}${OUTSUFFIX}.pgs.txt"
    pop_pgs_out="/disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/mcs/nofilter/pop_pgs/${PHENONAME}/${METHOD}/${ANCESTRY}"


    echo "Output to ${pop_pgs_out}/${EFFECT}${OUTSUFFIX}_PCadjusted"
    pgs.py ${pop_pgs_out}/${EFFECT}${OUTSUFFIX}_PCadjusted \
    --pgs ${scoresout} \
    --covar $COVAR \
    --gen_models 1 \
    --phenofile ${pheno_out}/pheno_withPCs.pheno \
    --scale_pgs \
    --scale_phen | tee "/disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/mcs/nofilter/pop_pgs/logs/${PHENONAME}_${EFFECT}${OUTSUFFIX}_PCadjusted_${ANCESTRY}_full.reg.log"

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

    withinfam_pred $direct_weights \
        "direct" "$PHENONAME" \
        "$OUTSUFFIX" "$BINARY" "$METHOD" "$ANCESTRY"

    if [ "$POPULATION" == "dir_pop" ]; then

        population_weights="/disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/prscs_weights/nofilter/${METHOD}/${PHENONAME}/population/weights/meta_weights.snpRes"
        withinfam_pred $population_weights \
            "population" "$PHENONAME" \
            "$OUTSUFFIX" "$BINARY" "$METHOD" "$ANCESTRY"
    fi
    
}