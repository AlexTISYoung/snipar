
phen='/var/genetics/data/mcs/private/latest/processed/phen'
DIR='/var/genetics/data/mcs/private/latest/raw/gen/NCDS_SFTP_1TB_1'


for maf in f721 nof21; do
    for type in europe_filter_maf05 non_europe_filter_maf05 filter_maf05; do
        python fGWAS.py "${DIR}/genotype/${type}/imputed_parents/pcs_${maf}"  "${phen}/simulated_stratification_phenotype.txt" --outprefix "${phen}/simulated_stratification_phenotype_${type}_pcs_${maf}" --bed "${DIR}/genotype/${type}/by_chrom/chr21" --covar "${DIR}/genotype/whole_genome_pc.csv" > "${phen}/simulated_stratification_phenotype_${type}_pcs_${maf}.log"
    done
done

for maf in f721 nof21; do
    for type in europe_filter_maf05 non_europe_filter_maf05 filter_maf05; do
        echo =========================================================
        echo =========================================================
        echo =========================================================
        echo "${phen}/simulated_stratification_phenotype_${type}_pcs_${maf}.log"
        echo =========================================================
        echo =========================================================
        echo =========================================================
        cat "${phen}/simulated_stratification_phenotype_${type}_pcs_${maf}.log"
    done
done




import h5py as hf
import numpy as np
addresses = {
'europe_nof':'/var/genetics/data/mcs/private/latest/processed/phen/simulated_stratification_phenotype_europe_filter_maf05_pcs_nof21.sumstats.hdf5',
'europe_f':'/var/genetics/data/mcs/private/latest/processed/phen/simulated_stratification_phenotype_europe_filter_maf05_pcs_f721.sumstats.hdf5',
'non_europe_nof':'/var/genetics/data/mcs/private/latest/processed/phen/simulated_stratification_phenotype_non_europe_filter_maf05_pcs_nof21.sumstats.hdf5',
'non_europe_f':'/var/genetics/data/mcs/private/latest/processed/phen/simulated_stratification_phenotype_non_europe_filter_maf05_pcs_f721.sumstats.hdf5',
'whole_sample_nof':'/var/genetics/data/mcs/private/latest/processed/phen/simulated_stratification_phenotype_filter_maf05_pcs_nof21.sumstats.hdf5',
'whole_sample_f':'/var/genetics/data/mcs/private/latest/processed/phen/simulated_stratification_phenotype_filter_maf05_pcs_f721.sumstats.hdf5',
}
files = {key:hf.File(val) for key,val in addresses.items()}
ses = {key:np.array(val['estimate_ses']) for key,val in files.items()}
effects = {key:np.array(val['estimate']) for key,val in files.items()}
zs = {key:effects[key]/ses[key] for key in files}
zvars = {key:np.nanvar(zs[key], 0) for key in files}
for key in zvars:
    print(key, zvars[key])




europe_nof [0.9568446 0.9142382 1.01344  ]
europe_f [0.9560069  0.91048557 1.0115039 ]
non_europe_nof [0.98201025 1.0787051  0.92175025]
non_europe_f [0.9806894  1.0359663  0.90396047]
whole_sample_nof [0.9627659  0.93808347 1.0760192 ]
whole_sample_f [0.9708536 1.009584  1.0296088]