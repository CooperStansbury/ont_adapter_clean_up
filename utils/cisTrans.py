import pandas as pd
import numpy as np
import itertools
import os
import sys
import math


def getCTFast(df):
    """A 'lightweight' function to clique-expand an incidence matrix and 
    count the fraction of cis- or trans- contacts for each read
    """
    newRows = []
    
    for (read_name, read_df) in df.groupby("read_base", as_index=False):
        if len(read_df) <= 1:
            continue
        
        pairs = list(
            read_df.sort_values(["read_start"], ascending=True)
            .assign(pos_on_read=lambda x: np.arange(len(x)))
            .itertuples()
        )
        
        order = len(pairs)
        n = 0
        nCis = 0
        nTrans = 0
        
        for (align_1, align_2) in itertools.combinations(pairs, 2):
            n += 1
            if align_1.chrom == align_2.chrom:
                nCis += 1
            else:
                nTrans += 1
        
        row = {
            'order' : order, 
            'nPairs' : n,
            'nCis' : nCis,
            'nTrans' : nTrans,
            'pCis' : nCis / n,
            'pTrans' : nTrans / n,
        }
        
        newRows.append(row)
        
    return pd.DataFrame(newRows)