#!/bin/bash

# write the directory project names to be included in the search
projects=('covid19_project' \
          'signed_graphs_project' \
          'Hayashi_P10 IHEC_project' \
          'Nakajima_drosophila_velocito' \
          'nakajima_drosophila_project_v2' \
          'ChIP-seq_imputation' \
          'nakajima_drosophila_ageing')

#create a directory with the date
mkdir -p $(date +"%Y-%m-%d")

cd $(date +"%Y-%m-%d")

for project in "${projects[@]}"; do

    # Dedicated directory for the project
    mkdir -p ${project}

    # Download scrips in Jupyter, bash and python
    # -a: "archive", sync recursively, synlinks, permissions
    # -v: verbose
    rsync -av --progress \
    --include='*.ipynb' \
    --include='*.sh' \
    --include='*.py' \
    --include='*/' \
    --exclude='*' \
    nagai@euphonium:/home/nagai/${project} .
done
echo "Finished!"
