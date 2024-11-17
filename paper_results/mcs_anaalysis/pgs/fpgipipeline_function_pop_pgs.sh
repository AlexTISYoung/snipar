#!/usr/bin/env bash

# within_family_path="/var/genetics/proj/within_family/within_family_project"
# snipar_path="/var/genetics/proj/within_family/snipar_simulate"
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
    COVAR="/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/phen/PCs.txt"
    # PHENOFILE="/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/phen/phenotypes_eur.txt"
    
    # mkdir -p $RAWPATH/phen/${PHENONAME}
    echo $OUTPATH/pop_pgs/${EFFECT}${OUTSUFFIX}
    echo xxxxxx
    # echo $WTFILE
    mkdir -p $OUTPATH/pop_pgs
    mkdir -p /disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/mcs/nofilter/pop_pgs/logs
    # mkdir -p $within_family_path/processed/fpgs/${PHENONAME}
    
       ## format weight files for fpgs                                             
    # python format_weights.py \
    #     $WTFILE \
    #     --chr 0 --pos 2 --rsid 1 --a1 3 --a2 4 --beta 5 \
    #     --sep "delim_whitespace" \
    #     --outfileprefix /disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/prscs_weights/${METHOD}/${PHENONAME}/${PHENONAME}_${EFFECT}_fpgs_formatted \
    #     --sid-as-chrpos \
    #     --prscs

    # generate pheno file
    pheno_out="$RAWPATH/phen/${PHENONAME}/${ANCESTRY}"
    # mkdir -p ${pheno_out}
    # python $within_family_path/scripts/fpgs/format_pheno.py \
    #     $PHENOFILE \
    #     --iid IID --fid FID --phenocol $PHENONAME \
    #     --outprefix ${pheno_out}/pheno  \
    #     --sep "delim_whitespace" \
    #     --binary $BINARY

    bedfilepath="/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/bgen/tmp/chr@.dose"
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
        --scale_pgs #| tee $OUTPATH/${EFFECT}${OUTSUFFIX}.log
    echo ssss

    scoresout="$OUTPATH/pop_pgs/${EFFECT}${OUTSUFFIX}.pgs.txt"
    pop_pgs_out="/disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/mcs/nofilter/pop_pgs/${PHENONAME}/${METHOD}/${ANCESTRY}"
    # mkdir -p $pop_pgs_out

    # ## run fPGI regression
    # echo "Run fPGI regression..."
    # python pgs_population_eff.py ${pop_pgs_out}/${EFFECT}${OUTSUFFIX} \
    #     --pgs ${scoresout} \
    #     --phenofile ${pheno_out}/pheno.pheno \
    #     --scale_pgs \
    #     --scale_phen | tee "/disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/mcs/pop_pgs/logs/${PHENONAME}_${EFFECT}${OUTSUFFIX}_${ANCESTRY}_full.reg.log"
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

    # processed_dir="/var/genetics/data/mcs/private/latest/processed/pgs/fpgs/${PHENONAME}/${METHOD}/${ANCESTRY}"
    RAWPATH="/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed"
    direct_weights="/disk/genetics/ukb/jguan/ukb_analysis/output/prs_analysis/prscs_weights/nofilter/${METHOD}/${PHENONAME}/direct/weights/meta_weights.snpRes"
    pheno_out="$RAWPATH/phen/${PHENONAME}/${ANCESTRY}"
    # covar_fid="/var/genetics/data/mcs/private/latest/raw/downloaded/NCDS_SFTP_1TB_1/imputed/phen/covar_pedigfid.txt"

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