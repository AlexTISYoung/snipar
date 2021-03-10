import pandas as pd
import numpy as np
import glob
import re
import argparse


def index_containing_substring(the_list, substring):
    for i, s in enumerate(the_list):
        if substring in s:
              return i
    return -1


def read_log(path, runname):

    print(f"Run Name: {runname}")

    with open(path, 'r') as f:
        results = f.readlines()
        
    start_idx = index_containing_substring(results, "Start time:")
        
    nobs_final = int(re.findall(r'\d+', results[start_idx+3])[0])
    transform_effect = results[start_idx+4].split(':')[1].replace(r' ', '')[:-1]
    
    estimates = eval(re.sub('Estimates:', '', results[start_idx+6])[:-1])
    h2_dir = estimates['v1']
    h2_indir = estimates['v2']
    r = estimates['r']
    
    se = re.findall(r'(\d+.\d+|nan)', 
                re.sub('Standard Errors:', '', results[start_idx+7])[:-1]
                .strip()
                .strip('[')
                .strip(']'))
    h2_dir_se = float(se[0])
    h2_indir_se = float(se[1])
    r_se = float(se[2])
    
    
    solver_output = ''
    for output in results[16:25]:
        solver_output += output
        
    estimation_time_prejk = re.sub('Estimation time:', '', results[start_idx+17][:-1]).strip().split(':')
    estimation_time_prejk = f"{estimation_time_prejk[0]} hours, {estimation_time_prejk[1]} minutes, {estimation_time_prejk[2]} seconds"
    
    # If jkse is done
    if len(results) > start_idx + 18:

        jk_blocksizes = int(re.findall(r'\d+', results[start_idx + 18])[0])
        jk_ncores = int(re.findall(r'\d+', results[start_idx + 19])[0])
        
        jkse = re.findall(r'(\d+.\d+|nan)', 
                re.sub('Standard Errors:', '', results[start_idx + 20])[:-1]
                .strip()
                .strip('[')
                .strip(']'))
        h2_dir_jkse = float(jkse[0])
        h2_indir_jkse = float(jkse[1])
        r_jkse = float(jkse[2])
        
        estimation_time_postjk = re.findall(r'\d+', results[start_idx + 21])
        estimation_time_postjk = f"{estimation_time_postjk[0]} hours, {estimation_time_postjk[1]} minutes, {estimation_time_postjk[2]} seconds"
    else:
        jk_blocksizes = ""
        jk_ncores = ""
    
        h2_dir_jkse = ""
        h2_indir_jkse = ""
        r_jkse = ""
        
        estimation_time_postjk = ""
        
    outdf = pd.DataFrame(
        {
            'run' : runname,

            'solver_output' : solver_output,

            'nobs_final' : nobs_final,
            'transform_effect' : transform_effect,

            'h2_dir' : h2_dir,
            'h2_indir' : h2_indir,
            'r' : r,

            'h2_dir_se' : h2_dir_se,
            'h2_indir_se' : h2_indir_se,
            'r_se' : r_se,

            'estimation_time_prejk' : estimation_time_prejk,

            'jk_blocksizes' : jk_blocksizes,
            'jk_ncores' : jk_ncores,
            'h2_dir_jkse' : h2_dir_jkse,
            'h2_indir_jkse' : h2_indir_jkse,
            'r_jkse' : r_jkse,

            'estimation_time_postjk' : estimation_time_postjk

        },

        index = [0]
    )
        
    return outdf
    
    
def read_all_files(files, runnames):

    dfout = pd.DataFrame(
        columns = ['solver_output', 'nobs_final', 'transform_effect',
                    'h2_dir', 'h2_indir', 'r',
                    'h2_dir_se', 'h2_indir_se',  'r_se',
                    'estimation_time_prejk', 'jk_blocksizes', 'jk_ncores',
                    'h2_dir_jkse',  'h2_indir_jkse', 'r_jkse', 'estimation_time_postjk']
    )
    paths = glob.glob(files)
    for path, name in zip(paths, runnames):
        df_path = read_log(path, name)
        dfout = dfout.append(df_path)
        
    return dfout.reset_index(drop = True)

    

if __name__ == '__main__':
    
    parser=argparse.ArgumentParser()
    parser.add_argument('filenames', type = str,
                       help = 'Input log file names. Format should be glob like')
    parser.add_argument('--outpath', type = str,
                       help = 'Full path of outputted csv file. Include .csv in name')
    parser.add_argument('--runnames', type = str,
                       help = '''Comma delimited list of run names. Eg: "1,2,3,4,5"
                       If nothing is provided or if length is not equal to
                       length of filenames, the path is the runname.
                       ''')
    
    args = parser.parse_args()
    
    if args.runnames is not None:
        runnames = [item.strip() for item in args.runnames.split(',')]
    else:
        runnames = glob.glob(args.filenames)

    
    dfout = read_all_files(args.filenames,
                         runnames)
    
    dfout.to_csv(args.outpath,
                index = False)
    