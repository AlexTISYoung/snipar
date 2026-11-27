from numba import njit
import numpy as np
from os import path
import snipar

@njit
def pos_to_cM(pos, boundaries, cM_pos):
    indices = np.searchsorted(boundaries, pos, side="right") - 1
    output = np.full(pos.size, np.nan)
    for i in range(pos.size):
        if 0 <= indices[i] < cM_pos.size:
            output[i] = cM_pos[indices[i]]
    return output

def decode_map_from_pos(chrom,pos):
    decode_map_path = path.join(path.dirname(snipar.__file__), f'util_data/decode_map/chr_{chrom}.gz')
    map = np.loadtxt(decode_map_path, dtype=float, skiprows=1)
    boundaries = np.hstack((np.array(map[0, 0], dtype=np.int_),np.array(map[:, 1], dtype=np.int_)))
    return pos_to_cM(pos, boundaries, map[:, 2])


def get_map_positions(mapfile, gts, min_map_prop=0.5):
    with open(mapfile) as stream:
        header = stream.readline().split()
    if "pposition" not in header or "gposition" not in header:
        raise ValueError("Map file must contain columns pposition and gposition")

    bp_pos, cm_pos = np.loadtxt(
        mapfile, usecols=(header.index("pposition"), header.index("gposition")),
        skiprows=1, ndmin=2, unpack=True,
    )
    if not cm_pos.size:
        raise ValueError("Map file must contain at least one position")
    if not np.isfinite(cm_pos).all():
        raise ValueError("Map cannot have non-finite positions")
    if (np.diff(bp_pos) <= 0).any():
        raise ValueError("Physical positions in the map must be strictly increasing")
    if cm_pos[0] < 0:
        raise ValueError("Map file cannot have negative values")
    if (np.diff(cm_pos) < 0).any():
        raise ValueError("Map not monotonic. Please make sure input is ordered correctly")
    if cm_pos.size > 1 and cm_pos[0] == cm_pos[-1]:
        raise ValueError("Map file has no variation")
    if cm_pos[-1] > 5000:
        raise ValueError("Maximum value of map too large")

    in_map = np.isin(gts.pos, bp_pos)
    prop_in_map = np.mean(in_map) if in_map.size else 0
    if prop_in_map < min_map_prop:
        raise ValueError(f"Only {100 * prop_in_map:.0f}% of SNPs have genetic positions in "
                         f"{mapfile}. Need at least {100 * min_map_prop:.0f}%")
    print(f"Found genetic map positions for {100 * prop_in_map:.0f}% of SNPs in {mapfile}")
    if prop_in_map < 1:
        print("Linearly interpolating genetic map for SNPs not in input map")
    # Retain interpolation between the common SNPs.
    shared_positions = np.unique(gts.pos[in_map])
    return np.interp(gts.pos, shared_positions, cm_pos[np.searchsorted(bp_pos, shared_positions)])

def map_from_bed(bedfile, chrom):
    bimfile = bedfile.removesuffix(".bed")+'.bim'
    print('Attempting to read map from ' + bimfile)
    bim = np.loadtxt(bimfile, usecols=(1, 2, 3), dtype=str, ndmin=2)
    bim_map = np.array(bim[:,1],dtype=float)
    bim_pos = np.array(bim[:,2],dtype=int)
    # Check for NAs
    if not bim_map.size:
        raise ValueError('BIM file does not contain any SNPs')
    if np.all(bim_map == 0) or (bim_map.size > 1 and np.var(bim_map) == 0):
        print('Map information not found in bim file.')
        print('Using default map (decode sex averaged map on GRCh38 coordinates)')
        map = decode_map_from_pos(chrom, bim_pos)
        pc_mapped = str(round(100*(1-np.mean(np.isnan(map))),2))
        print('Found map positions for '+str(pc_mapped)+'% of SNPs')
    else:
        if not np.isfinite(bim_map).all():
            raise ValueError('Map cannot have non-finite positions')
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
