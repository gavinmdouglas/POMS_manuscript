#!/usr/bin/python3

import argparse
import os
import random
import sys


def read_fasta(filename, cut_header=False):

    '''Read in FASTA file and return dictionary with each independent sequence
    id as a key and the corresponding sequence string as the value.
    '''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    # Read in FASTA line-by-line.
    with open(filename, "r") as fasta:

        for line in fasta:

            # If header-line then split by whitespace, take the first element,
            # and define the sequence name as everything after the ">".
            if line[0] == ">":

                if cut_header:
                    name = line.split()[0][1:]
                else:
                    name = line[1:]

                name = name.rstrip("\r\n")

                # Intitialize empty sequence with this id.
                seq[name] = ""

            else:
                # Remove line terminator/newline characters.
                line = line.rstrip("\r\n")

                # Add sequence to dictionary.
                seq[name] += line

    return seq



def main():

    parser = argparse.ArgumentParser(

            description="Randomize order of seqs in FASTA",

epilog='''Usage example:

python reorder_fasta_seqs.py -i IN.fasta -o OUT.fasta

''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="PATH", type=str,
                        help="Path to input FASTA", required=True)

    parser.add_argument("-s", "--seed", metavar="PATH", type=int,
                        help="Random seed", default=1174, required=False)

    parser.add_argument("-o", "--output", metavar="PATH", type=str,
                        help="Path to output FASTA", required=True)                     

    args = parser.parse_args()
  
    in_fasta = read_fasta(args.input)

    random.seed(args.seed)

    seq_ids = list(in_fasta.keys())

    random.shuffle(seq_ids)

    # Write output file.
    outfile = open(args.output, "w")
    for seq_id in seq_ids:
        print(">" + seq_id, file=outfile)
        print(in_fasta[seq_id], file=outfile)
    outfile.close()

if __name__ == '__main__':
    main()

