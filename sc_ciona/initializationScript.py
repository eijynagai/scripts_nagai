#!/usr/bin/env python

import sys
import os
import re
from pathlib import Path


if len(sys.argv) != 5:
    print('The script is called with %i argument(s). Please follow the example below.' % (len(sys.argv) -1))
    print('>python {} [input_bc_matrix_directory] [output_directory] [project] [resolution]'.format(str(sys.argv[0])))
    sys.exit()

indir = sys.argv[1]
outdir = sys.argv[2]
proj = sys.argv[3]
res = float(sys.argv[4])

# check if user include '/' at end of the string and remove it
outdir = re.sub('/$', '', outdir)
outdir+= '/' + str(proj) + '_res' + str(res)

# check if indir exists
if not(os.path.isdir(indir)):
    print('Input directory does not exist!')
    sys.exit()

# create output directory
Path(outdir).mkdir(parents=True, exist_ok=True)

#------ start parameters
param = {}


# input all values from parameters
param = {
    'cr_dir'            : (indir,       'Path to data input bc matrix'),
    'dirname'           : (outdir,      'Path to output files'),
    'proj_name'         : (proj,        'Project name'),
    'qtd_min_cells'     : (3,           'Minimun number of cells per feature'),
    'qtd_min_features'  : (500,         'Minimun number of features per cell'), #1000
    'min_nFeature_RNA'  : (500,         'Minimun number of counts per feature'), #500
    'max_nFeature_RNA'  : (4500,        'Maximun number of counts per feature'),
    'max_perc_mt'       : (8,           'Maximun percentage allowed for mithocondrial genes'), #5
    'norm_scale'        : (10000,       '------'),
    'nFeatures_var'     : (3000,        'Number of variable features'),#2000
    'topx'              : (10,          'Number of top genes to represent'),
    'nPCS'              : (10,         'Number of PCs to be calculated'),
    'PCx'               : (10,         'PC threshold'),
    'RES'               : (res,         'Resolution for granularity'),
    'logFC'             : (0.25,        'log fold change'),
    'qvThres'           : (0.001,       'qv threshold'),
    'findmarkers_test'  : ('wilcox',    'test.use parameter for FindAllMarkers'),  #roc, wilcox
    'markers_path'      : ('/home/nagailae/projects/sc_ciona/seurat3_v1/script/marker_genes.txt', 'marker genes'),
    'colorplots'        : ('lightgrey', 'Color for plots')
}
#------ end parameters


# create parameters file to be used in R script
w = open(outdir + '/cfg_params.txt', 'w')
for key, val in param.items():
    if (isinstance(val[0],int) or isinstance(val[0],float)):
        w.write('## {}\n{} = {}\n\n'.format(val[1], key, val[0]))
    else:
        w.write('## {}\n{} = \'{}\'\n\n'.format(val[1], key, val[0]))

w.close()
