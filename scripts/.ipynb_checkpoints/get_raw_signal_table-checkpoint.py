import pandas as pd
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pod5 as p5
import pysam


if __name__ == "__main__":
    pod5_path = sys.argv[1]  
    outpath = sys.argv[2]
    
    records = []
    
    with p5.Reader(pod5_path) as reader:
        for read in reader.reads():
            # Get the signal data and sample rate
            sample_rate = read.run_info.sample_rate
            signal = read.signal
    
            signal_row = {
                'read_name' : str(read.read_id),
                'sample_rate' : float(sample_rate),
                'signal_length' : len(signal),
                'acquisition_start_time' : read.run_info.acquisition_start_time,
                'experiment_name' : str(read.run_info.experiment_name),
                'signal' : ";".join([str(x) for x in signal])
            }
    
            records.append(signal_row)
    
    df = pd.DataFrame(records)
    df.to_parquet(outpath, index=False)
    
    
        


    