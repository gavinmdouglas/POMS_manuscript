#!/usr/bin/python3

import argparse
import os
import sys
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(

            description="Read in and link two different tables (whitespace delimited) in long-format.\n\n"
                        "These are similar to mapfiles and the best way to "
                        "understand this script is through an example:\n\n"
                        "If this is FILE1:\n"
                        "AC_Bin_0   AC_0\n"
                        "AC_Bin_0   AC_4\n"
                        "AC_Bin_0   AC_17\n"
                        "AC_Bin_417 AC_4\n"
                        "\n\n"
                        "And this is FILE2:\n"
                        "AC_0   K01872\n"
                        "AC_4   K00606\n"
                        "AC_17  K02483\n"
                        "AC_17  K18344\n"
                        "\n\n"
                        "Then the output file would look like this:\n"
                        "Group  K00606  K01872  K02483  K18344\n"
                        "AC_Bin_0   1   1   1   1\n"
                        "AC_Bin_417 0   1   0   0\n\n",

epilog='''Usage example:

python link_long_tables.py -1 FILE1 -2 FILE2 -o OUTPUT

''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f1", "--file1", metavar="PATH", type=str,
                        help="Path to first input mapfile.", required=True)

    parser.add_argument("-f2", "--file2", metavar="PATH", type=str,
                        help="Path to second input mapfile.", required=True)

    parser.add_argument("-c", "--column", metavar="STR", type=str,
                        help="Name for first column in output file (\"Group\" "
                             "by default).", default="Group")

    parser.add_argument("-o", "--output", metavar="PATH", type=str,
                        help="Path to output table.", required=True)                     

    args = parser.parse_args()


    uniq_col = set()
    file2_groupings = defaultdict(set)

    with open(args.file2, "r") as FILE2:
        for line in FILE2:
            if len(line) == 0:
                continue
            line_split = line.split()
            file2_groupings[line_split[0]].add(line_split[1])
            uniq_col.add(line_split[1])


    uniq_row = set()
    file1_groupings = {}

    with open(args.file1, "r") as FILE1:
        for line in FILE1:
            if len(line) == 0:
                continue
        
            line_split = line.split()

            if line_split[0] not in file1_groupings.keys():
                file1_groupings[line_split[0]] = defaultdict(int)

            uniq_row.add(line_split[0])

            mapped_vars = file2_groupings[line_split[1]]

            for mapped_var in mapped_vars:
                file1_groupings[line_split[0]][mapped_var] += 1

    uniq_row = sorted(list(uniq_row))
    uniq_col = sorted(list(uniq_col))

    outfile = open(args.output, "w")
    print("\t".join([args.column] + uniq_col), file=outfile)
    for rowname in uniq_row:
        row_line = rowname
        
        for colname in uniq_col:
            row_line += "\t" + str(file1_groupings[rowname][colname])

        print(row_line, file=outfile)

if __name__ == '__main__':
    main()

