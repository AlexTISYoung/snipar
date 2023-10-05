import h5py
import numpy as np 
from snipar.utilities import *

phased_dir = '/disk/genetics/ukb/jguan/ukb_analysis/output/parent_imputed/v2/'
unphased_dir = '/disk/genetics/ukb/jguan/ukb_analysis/output/parent_imputed/v4/'
outdir = '/disk/genetics/ukb/alextisyoung/imputed/'

snp_keys = ['estimated_genotyping_error','imputed_par_gts',
            'mendelian_error_ratio','non_duplicates','parent_ratio-backup',
            'pos','ratio_ibd0','sib_ratio-backup']

non_snp_keys = ['families','parental_status','pedigree']

for chrom in range(1,22):
    hdf5_phased = h5py.File(phased_dir+'chr_'+str(chrom)+'.hdf5', 'r')
    hdf5_unphased = h5py.File(unphased_dir+'chr_'+str(chrom)+'.hdf5', 'r')
    hdf5_out = h5py.File(outdir+'chr_'+str(chrom)+'.hdf5', 'w')
    # Merge
    phased_bim = convert_str_array(hdf5_phased['bim_values'])
    unphased_bim = convert_str_array(hdf5_unphased['bim_values'])
    phased_bim_alleles = np.array([str.split(x,',') for x in phased_bim[:,5]])
    combined_bim = np.vstack((unphased_bim[:,[0,1,3,4,5]],np.hstack((phased_bim[:,[2,0,3]],phased_bim_alleles))))
    hdf5_out['bim_columns'] = ['chrom','rsid','pos','allele1','allele2']
    # Remove duplicates from combined bim
    unique_snps = np.unique(combined_bim[:,1],return_index=True)
    unique_snps_pos = np.array(combined_bim[unique_snps[1],2],dtype=int)
    pos_sort = np.argsort(unique_snps_pos)
    unique_merge_bim = np.vstack((np.repeat(chrom,unique_snps[0].shape[0]),
                                    unique_snps[0][pos_sort],unique_snps_pos[pos_sort])).T 
    merge_alleles = np.zeros((unique_snps[0].shape[0],2),dtype=unique_merge_bim.dtype)
    # Map to merged bim
    merge_bim_dict = make_id_dict(unique_merge_bim,1)
    unphased_map = np.array([merge_bim_dict[x] for x in unphased_bim[:,1]])
    unphased_map_sort = np.argsort(unphased_map)
    unphased_map = unphased_map[unphased_map_sort]
    phased_map = np.array([merge_bim_dict[x] for x in phased_bim[:,1]])
    phased_map_sort = np.argsort(phased_map)
    phased_map = phased_map[phased_map_sort]
    # Alleles
    merge_alleles[unphased_map,:] = unphased_bim[unphased_map_sort,4:6]
    merge_alleles[phased_map,:] = phased_bim_alleles[phased_map_sort,:]
    unique_merge_bim = np.hstack((unique_merge_bim,merge_alleles))
    hdf5_out['bim_values'] = encode_str_array(unique_merge_bim)
    # Fams
    phased_fams = convert_str_array(hdf5_phased['families'])
    unphased_fams = convert_str_array(hdf5_unphased['families'])
    if not np.array_equal(phased_fams,unphased_fams):
        raise ValueError('Fam IDs do not match')
    for key in hdf5_unphased.keys():
        print(key)
        if key in snp_keys:
            key_shape = hdf5_phased[key].shape
            if len(key_shape)==1:
                hdf5_out.create_dataset(key, shape=(unique_snps[0].shape[0]), dtype=hdf5_phased[key].dtype)
                hdf5_out[key][unphased_map] = np.array(hdf5_unphased[key])[unphased_map_sort]
                hdf5_out[key][phased_map] = np.array(hdf5_phased[key])[phased_map_sort]
            else:
                hdf5_out.create_dataset(key, shape=(key_shape[0],unique_snps[0].shape[0]), dtype=hdf5_phased[key].dtype)
                hdf5_out[key][:,unphased_map] = np.array(hdf5_unphased[key])[:,unphased_map_sort]
                hdf5_out[key][:,phased_map] = np.array(hdf5_phased[key])[:,phased_map_sort] 
        if key in non_snp_keys:
            hdf5_out.create_dataset(key, data=hdf5_phased[key])
    hdf5_out.close()
    print(chrom) 