#!/usr/bin/python3

import argparse
from collections import defaultdict
import gzip
import os


def main():

    parser = argparse.ArgumentParser(

        description="Parse MIDAS genome annotations to a single table.",
        epilog='''Usage example:

    python parse_midas_gene_counts.py -i IN_DIRECTORY -c figfam -o OUT_TABLE

''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--in_dir", metavar="PATH", type=str,
                        help="Path to MIDAS genome folder", required=True)

    parser.add_argument("-c", "--category", metavar="CATEGORY", type=str,
                        help="Functional category", required=True)

    parser.add_argument("-o", "--output", metavar="PATH", type=str,
                        help="Path to output table", required=True)                     

    args = parser.parse_args()

    genome_folders = os.listdir(args.in_dir)

    genome_annot = defaultdict(dict)

    all_func = set()

    for genome_id in genome_folders:

        genome_annot[genome_id] = defaultdict(int)

        features_file = os.path.join(args.in_dir, genome_id,
                                     "centroid_functions.txt.gz")

        with gzip.open(features_file,'rt') as f:

            for line in f:
                
                line = line.rstrip()

                line_split = line.split()

                func_id = line_split[1]

                if line_split[2] == args.category:

                    genome_annot[genome_id][func_id] += 1

                    all_func.add(func_id)


    all_func = sorted(list(all_func))

    outfile = open(args.output, 'w')

    headerline = ['function']

    for genome_id in genome_folders:

        headerline += [str(genome_id)]

    print('\t'.join(headerline), file=outfile)

    for func in all_func:
        
        out_row = [func]

        for genome_id in genome_folders:

            out_row += [str(genome_annot[genome_id][func])]

        print('\t'.join(out_row), file=outfile)

    outfile.close()


if __name__ == '__main__':
    main()

