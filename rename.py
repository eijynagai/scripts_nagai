#!/bin/env python

# Program to receive a path and collect the name

import os

def file_name(url, verbose=True):

    if verbose == True:
        print("Path: {}".format(url))

    name = url.split('/')[-1]
    name = name if '.' not in name else '.'.join(name.split('.')[:-1])

    print("File name is : {}".format(name))
    return name

path="test/path/to/somewhere.I/dont/know/bigwigwigwig.bw.BigWig"

def main():
    file_name(path)
    

if __name__=="__main__":
    main()
    
