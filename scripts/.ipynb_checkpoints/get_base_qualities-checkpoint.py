import sys
import pandas as pd
import pysam
import numpy as np


if __name__ == "__main__":
    fastq_path = sys.argv[1]  
    outpath = sys.argv[2]
    
    records = []
    with pysam.FastxFile(fastq_path, "r") as f:
        for read in f:
            read_name = read.name
            quals = read.get_quality_array()
    
            row = {
                'read_name' : read_name,
                'read_length' : len(read.sequence),
                'mean_base_quality' : np.mean(quals),
                'std_base_quality' : np.std(quals),
                'median_base_quality' : np.median(quals),
                'min_base_quality' : np.min(quals),
                'max_base_quality' : np.max(quals),
            }
            records.append(row)
    
    df = pd.DataFrame(records)
    df.to_parquet(outpath, index=False)
    
        


    