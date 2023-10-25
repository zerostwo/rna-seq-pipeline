import pandas as pd
import sys

def merge_counts(count_files):
    dfs = []
    for file in count_files:
        # Read in file, skip first row and set index to Geneid
        df = pd.read_csv(file, sep="\t", skiprows=1)
        df = df.set_index("Geneid")
        
        # Extract sample name from filename
        sample_name = file.split("/")[-1].split(".")[0]
        
        # Select only the count column and rename it to the sample name
        count_column = df.columns[-1]  # Selects the last column
        df = df[[count_column]]
        df.columns = [sample_name]
        
        dfs.append(df)
    
    # Merge dataframes on Geneid
    merged = pd.concat(dfs, axis=1)

    return merged

if __name__ == "__main__":
    count_files = sys.argv[1:]
    matrix = merge_counts(count_files)
    print(matrix.to_csv(sep="\t"))

