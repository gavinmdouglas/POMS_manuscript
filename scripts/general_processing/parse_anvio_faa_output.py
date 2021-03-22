#!/usr/bin/python3

import argparse
import os
import sys
import textwrap

def main():

    parser = argparse.ArgumentParser(

            description="Remove gap characters from FASTA, clean-up header, "
                        "and remove all bins with completeness/redundancy "
                        "outside specified range.",

epilog='''Usage example:

python parse_anvio_faa_output.py -f FASTA -s SUMMARY_FILE -d SET -o OUTPUT

''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="PATH", type=str,
                        help="Path to input FASTA", required=True)

    parser.add_argument("-s", "--summary", metavar="PATH", type=str,
                        help="Path to input anvio summary file", required=True)

    parser.add_argument("-d", "--dataset", metavar="STR", type=str,
                        help="String describing dataset to add to all " +
                             "headerlines", required=True)

    parser.add_argument("-c", "--completeness", metavar="FLOAT", type=float,
                        help="Min percent completeness per bin",
                        default=90.0)

    parser.add_argument("-r", "--redundancy", metavar="FLOAT", type=float,
                        help="Max percent redundancy per bin",
                        default=10.0)

    parser.add_argument("-o", "--output", metavar="PATH", type=str,
                        help="Path to output cleaned-up FASTA",
                        required=True)                     

    args = parser.parse_args()


    seqs = {}
    with open(args.fasta, "r") as fasta:
        for line in fasta:
            # Skip blank lines
            if len(line) == 0:
                continue
            elif line[0] == ">":
                    current_seq = args.dataset + "." + line.split()[0].replace(">", "")
                    seqs[current_seq] = ""
            else:
                seqs[current_seq] += line.rstrip().replace("-", "")


    bins2keep = set()
    summary_lc = 0

    with open(args.summary, "r") as summary:
        for line in summary:
            
            # Skip first line
            if summary_lc == 0:
                summary_lc += 1
                continue

            line_split = line.split()

            bin_name = args.dataset + "." + line_split[0]

            completeness = float(line_split[6])
            redundancy = float(line_split[7])

            if completeness >= args.completeness and redundancy <= args.redundancy:
                bins2keep.add(bin_name)


    # Write output file.
    outfile = open(args.output, "w")
    for bin_name in list(bins2keep):
        print(">" + bin_name, file=outfile)
        print(textwrap.fill(seqs[bin_name], width=50), file=outfile)
    outfile.close()

if __name__ == '__main__':
    main()

