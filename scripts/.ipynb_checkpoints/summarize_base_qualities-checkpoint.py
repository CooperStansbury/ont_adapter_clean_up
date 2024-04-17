import sys
import pandas as pd
import pysam
import os
import numpy as np

print(f"{pd.__version__=}")
print(f"{pd.__file__=}")

if __name__ == "__main__":
    file_path = sys.argv[1]  
    outpath = sys.argv[2]

    basename = os.path.basename(file_path)

    # load the dataframe
    df = pd.read_parquet(file_path)

    # summarize
    summary = {
        'basename' : basename,
        'n_reads' : df['read_name'].nunique(),
        'mean_read_length' : df['read_length'].mean(),
        'median_read_length' : df['read_length'].median(),
        'std_read_length' : df['read_length'].std(),
        'mean_mean_base_quality' : df['mean_base_quality'].mean(),
        'std_mean_base_quality' : df['mean_base_quality'].std(),
    }
    
    # create a single record dataframe and store
    summary = pd.DataFrame([summary])
    summary.to_parquet(outpath, index=False)
    
            


    