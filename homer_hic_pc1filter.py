#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""Description
author:     Luis AE Nagai
date:       2017-12-20
Aim:        Remove group of lines with same signal which have less than
            threshold size (recommended 1 Mb). This script remove possible
            small compartments when using HOMER scripts. Should be used in
            between runHiCpca.pl and findHiCCompartments.pl.
Input:      PC1.txt output file from HOMER runHiCpca.pl
Output:     PC1.txt file filtered

Last modification: 2017-12-22: include annotation
"""

import argparse, sys
__author__ = 'Luis AE Nagai'

parser = argparse.ArgumentParser(description="Remove group of lines with same signal which have less than threshold size. This script remove possible small compartments.")
parser.add_argument('-i','--input',help='Input file name',required=True)
parser.add_argument('-o','--output',help='Output file name',required=True)
parser.add_argument('-t','--threshold',help='Minimum size to consider a compartment',required=True)
args = parser.parse_args()

def main():

    # declare variables:
    binsize=0
    countcompsize = 0
    templist = []
    count_lines = 0

    # cleaning the output file to not append in older output
    open(args.output, 'w').close()

    with open(args.input) as fi:
        with open(args.output,"a") as fo:
            for line in fi:

                # If there is any header, shoudl start with '#'
                if line.startswith('#'):
                    fo.write(line)
                else:
                    # If first line, attribute values
                    if (count_lines == 0):
                        # Calculating the bin size in this experiment
                        binsize=int((line.strip('\n').split('\t')[3]))-int(line.strip('\n').split('\t')[2])
                        templist.append(line)
                        count_lines = 1
                        countcompsize = 1
                        signal=comp_status(float(line.strip('\n').split('\t')[-1]))
                    # From the second line
                    else:
                        # If the next line keep same signal, append
                        if comp_status(float(line.strip('\n').split('\t')[-1])) == signal:
                            templist.append(line)
                            count_lines += 1
                            countcompsize += 1
                        # If the next line shift the signal, check the size and export
                        else:
                            #print("size of compartment: {}, number of lines checked: {}".format(countcompsize,count_lines))
                            if (countcompsize*binsize) >= int(args.threshold):
                                for row in templist:
                                    fo.write(row)
                            templist=[]
                            templist.append(line)
                            signal=comp_status(float(line.strip('\n').split('\t')[-1]))
                            countcompsize = 1
                            count_lines +=1
            #including the last group
            if (countcompsize*binsize) >= int(args.threshold):
                for row in templist:
                    fo.write(row)
    print("Filtering  {}: process successfully done!".format(args.input))

def comp_status(signal):
    if signal >= 0:
        return 1
    elif signal < 0:
        return 0
    else:
        print ("Format of file is not Homer PC1.txt")
        sys.exit()


if __name__ == '__main__':
    main()
