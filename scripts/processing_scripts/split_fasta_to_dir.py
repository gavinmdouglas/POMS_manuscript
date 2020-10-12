#!/usr/bin/python3

import argparse
import re
import os


def main():

    parser = argparse.ArgumentParser(

            description="Reads in fasta and writes all seqs to a directory of "
                        "files: 1 per different sequence. Note that sequence headers "
                        "will be used to make file names - meaning they are assumed to be very basic.",
epilog='''Usage example:

python split_fasta_to_dir.py -i FASTA -o out_directory

''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to FASTA file", required=True)

    parser.add_argument("-o", "--outdir", metavar="PATH", type=str,
                        help="Path to output directory", required=True)                     

    args = parser.parse_args()

    # Make output directory if it doesn't exist.
    os.makedirs(args.outdir, exist_ok=True)

  
    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    seq_count = 0

    # Read in FASTA line-by-line.
    with open(args.fasta, "r") as fasta:

        for line in fasta:

            # If header-line then split by whitespace, take the first element,
            # and define the sequence name as everything after the ">".
            if line[0] == ">":

                if seq_count == 0:
                    seq_count = seq_count + 1
                else:
                    seq_out.close()

                name = line.split()[0][1:]

                seq_out = open(args.outdir + "/" + name + ".fna", "w")

            print(line, file=seq_out, end="")


    seq_out.close()
        

if __name__ == '__main__':
    main()

