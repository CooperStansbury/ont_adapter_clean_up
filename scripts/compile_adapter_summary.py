import pandas as pd
import sys
import re
import os
import numpy as np


def parse_alignment_file(file_path):
    """
    Parses an alignment file into a structured pandas DataFrame.

    Args:
        file_path (str): The path to the alignment file.

    Returns:
        pandas.DataFrame: A DataFrame containing the parsed alignment data.
    """

    data = []
    pattern = re.compile(r"(.+) \(read coords: (\d+)-(\d+), identity: (\d+\.\d+)%\)")

    with open(file_path, "r") as file:
        for line in file:
            match = pattern.search(line)
            if match:
                name, start, end, identity = match.groups()
                data.append({"Name": name, "Start": int(start), "End": int(end), "Identity": float(identity)})

    df = pd.DataFrame(data)
    return df



if __name__ == "__main__":
    output_path = sys.argv[1] 
    file_list = sys.argv[2:]  
    
    res = []
    for file_path in file_list:
        basename = os.path.basename(file_path)
        tmp = parse_alignment_file(file_path)
        tmp = tmp['Name'].value_counts().reset_index()
        tmp['basename'] = basename
        res.append(tmp)

    res = pd.concat(res)
    res.to_csv(output_path, index=False)
    

    
   
   
