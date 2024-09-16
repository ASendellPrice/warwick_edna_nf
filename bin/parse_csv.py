# parse_csv.py
import pandas as pd
import sys

def parse_csv(csv_file):
    df = pd.read_csv(csv_file)
    for index, row in df.iterrows():
        print(f"'{row['sampleID']}' {row['forward_reads']} {row['reverse_reads']}")

if __name__ == "__main__":
    parse_csv(sys.argv[1])
