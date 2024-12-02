#!/usr/bin/env bash

## --------------------------------------------------------------------------------------------------------
## args:
## PHENONAME=$1. phenotype name
## OUTSUFFIX=$2. suffix for output file (usually blank)
## BINARY=$3. 1 for binary, 0 for continuous outcome phenotype
## METHOD=$4. method used to produce sumstats (unified, robust, sibdiff, young)
## ANCESTRY=$5. ancestry for MCS validation sample. previously we ran for sas and eur, but now just eur.
## POPULATION=$6. "dir_pop" for both direct and population effect pgi, empty otherwise
## --------------------------------------------------------------------------------------------------------

set -e
source ${snipar_env}/bin/activate
source fpgipipeline_function_pop_pgs_sas.sh

# unified
main "bmi" "" "0" "unified" "sas" "dir_pop" 
main "height" "" "0" "unified" "sas" "dir_pop" 
main "ea" "" "0" "unified" "sas" "dir_pop" 

## new robust
main "bmi" "" "0" "new_robust" "sas"
main "height" "" "0" "new_robust" "sas"
main "ea" "" "0" "new_robust" "sas"

## sibdiff
main "bmi" "" "0" "sibdiff" "sas"
main "height" "" "0" "sibdiff" "sas"
main "ea" "" "0" "sibdiff" "sas"

## young
main "bmi" "" "0" "young" "sas"
main "height" "" "0" "young" "sas"
main "ea" "" "0" "young" "sas"

