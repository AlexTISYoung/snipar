from numba import njit, prange
import numpy as np
from os import path
import snipar
from snipar.utilities import make_id_dict
from snipar.gtarray import gtarray

@njit
def pos_to_cM(pos,boundaries, cM_pos):
    cM_out = np.zeros((pos.shape[0]), dtype=np.float_)
    cM_out[...] = np.nan
    current_seg = 0
    for i in range(pos.shape[0]):
        if boundaries[cM_pos.shape[0]] > pos[i] >= boundaries[0]:
            unmapped = True
            while unmapped:
                if boundaries[current_seg + 1] > pos[i] >= boundaries[current_seg]:
                    cM_out[i] = cM_pos[current_seg]
                    unmapped = False
                else:
                    current_seg += 1
    return cM_out

def decode_map_from_pos(chrom,pos):
    decode_map_path = path.join(path.dirname(snipar.__file__), f'util_data/decode_map/chr_{chrom}.gz')
    map = np.loadtxt(decode_map_path, dtype=float, skiprows=1)
    boundaries = np.hstack((np.array(map[0, 0], dtype=np.int_),np.array(map[:, 1], dtype=np.int_)))
    return pos_to_cM(pos, boundaries, map[:, 2])


# Read header of mapfile
def get_map_positions(mapfile,gts,min_map_prop = 0.5):
    map_file = open(mapfile,'r')
    map_header = map_file.readline()
    map_header = np.array(map_header.split(' '))
    map_header[len(map_header)-1] = map_header[len(map_header)-1].split('\n')[0]
    map_file.close()
    if 'pposition' in map_header and 'gposition' in map_header:
        bp_pos = np.loadtxt(mapfile,usecols = np.where(map_header=='pposition')[0][0], dtype=int, skiprows =1)
        pos_dict = make_id_dict(bp_pos)
        cm_pos = np.loadtxt(mapfile,usecols = np.where(map_header=='gposition')[0][0], dtype=float, skiprows =1)
        # Check for NAs
        if np.sum(np.isnan(cm_pos)) > 0:
            raise (ValueError('Map cannot have NAs'))
        if np.min(cm_pos) < 0:
            raise (ValueError('Map file cannot have negative values'))
        if np.var(cm_pos) == 0:
            raise (ValueError('Map file has no variation'))
        # Check ordering
        ordered_map = np.sort(cm_pos)
        if np.array_equal(cm_pos, ordered_map):
            pass
        else:
            raise (ValueError('Map not monotonic. Please make sure input is ordered correctly'))
        # Check scale
        if np.max(cm_pos) > 5000:
            raise (ValueError('Maximum value of map too large'))
        # Find positions of SNPs in map file
        map = np.zeros((gts.shape[1]),dtype=float)
        map[:] = np.nan
        in_map = np.array([x in pos_dict for x in gts.pos])
        # Check if we have at least 50% of SNPs in map
        prop_in_map = np.mean(in_map)
        if prop_in_map < min_map_prop:
            raise(ValueError('Only '+str(round(100*prop_in_map))+'% of SNPs have genetic positions in '+mapfile+'. Need at least '+str(round(100*min_map_prop))+'%'))
        print('Found genetic map positions for '+str(round(100*prop_in_map))+'% of SNPs in '+mapfile)
        # Fill in map values
        map[in_map] = cm_pos[[pos_dict[x] for x in gts.pos[in_map]]]
        # Linearly interpolate map
        if prop_in_map < 1:
            print('Linearly interpolating genetic map for SNPs not in input map')
            map = np.interp(gts.pos, gts.pos[in_map], map[in_map])
        return map
    else:
        raise(ValueError('Map file must contain columns pposition and gposition'))

def map_from_bed(bedfile, chrom):
    bimfile = bedfile.split('.bed')[0]+'.bim'
    print('Attempting to read map from ' + bimfile)
    bim = np.loadtxt(bimfile, usecols=(1,2,3), dtype=str)
    bim_map = np.array(bim[:,1],dtype=float)
    bim_pos = np.array(bim[:,2],dtype=int)
    # Check for NAs
    if np.var(bim_map) == 0:
        print('Map information not found in bim file.')
        print('Using default map (decode sex averaged map on Hg19 coordinates)')
        map = decode_map_from_pos(chrom, bim_pos)
        pc_mapped = str(round(100*(1-np.mean(np.isnan(map))),2))
        print('Found map positions for '+str(pc_mapped)+'% of SNPs')
    else:
        if np.sum(np.isnan(bim_map)) > 0:
            raise (ValueError('Map cannot have NAs'))
        if np.min(bim_map) < 0:
            raise (ValueError('Map file cannot have negative values'))
        # Check ordering
        ordered_map = np.sort(bim_map)
        if np.array_equal(bim_map, ordered_map):
            pass
        else:
            raise (ValueError('Map not monotonic. Please make sure input is ordered correctly'))
        # Check scale
        if np.max(bim_map) > 500:
            raise (ValueError('Maximum value of map too large'))
        map = bim_map
        print('Read map from bim file')
    return bim[:,0], map