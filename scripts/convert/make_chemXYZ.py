##########
# Script for converting basic 3-column x,y,z files to proper .xyz files
##########
import os
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=str, help="Path to directory with basic 3-column xyz file(s) to be converted.")
    parser.add_argument("-a", "--atom", type=str, required=False, default='C', help="Atom type to append as first column [default = carbon (C)]")
    parser.add_argument("-j", "--header", type=str, required=False, default='A properly formatted .xyz file', help="Header to put in file [default = 'A properly formatted .xyz file']")
    parser.add_argument("-o", "--outdir", type=str, required=False, default='.', help="Output location for properly-formatted .xyz files.")

    return parser.parse_args()

def main():
    args = parse_args()

    files = os.listdir(args.dir)

    for i in files:
        file_path = os.path.join(args.dir, i)
        f = pd.read_csv(file_path, sep='\t', header=None, skiprows=[0])  # Skip first row (assuming it's atom count)
        n_atoms = len(f)

        # Insert atom column at index 0
        f.insert(0, 'Atom', [args.atom] * n_atoms)

        # Define column names to maintain order
        num_cols = f.shape[1]  # Get the number of columns
        column_names = ['Atom'] + [f'Col_{j}' for j in range(1, num_cols)]  # Generate column names
        f.columns = column_names  # Assign structured column names

        # Create headers ensuring the structure remains intact
        header = pd.DataFrame([[str(n_atoms)] + [''] * (num_cols - 1)], columns=column_names)
        header_text = pd.DataFrame([[args.header] + [''] * (num_cols - 1)], columns=column_names)

        # Concatenate headers with data while ensuring column order remains fixed
        xyz_data = pd.concat([header, header_text, f], axis=0, ignore_index=True)

        # Print for verification
        #print(xyz_data)
        #print("===============")

        xyz_data.to_csv(f"{args.outdir}/{i}.xyz", sep='\t', doublequote=False, index=False, header=False)

if __name__ == "__main__":
    main()


