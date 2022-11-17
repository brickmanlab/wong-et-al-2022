#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd

if __name__ == '__main__':
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, help="Input folder for bed files", required=True)
    args = parser.parse_args()

    bed_files = glob.glob(f'{args.input}/*.bed')
    
    for f in bed_files:
        df = pd.read_csv(f, index_col=0, sep='\t', header=None)
        classes = df[4].unique()

        if len(classes) == 1:
            print(f'Skipping splitting file {os.path.basename(f)}')
        else:
            for klass in classes:
                new_filename = f.replace('.bed', f'_c{klass}.bed')
                df[df[4] == klass].to_csv(new_filename, sep='\t', header=None)
