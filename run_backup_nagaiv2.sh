#!/bin/bash
# program to automatically download specific files from the server
# author: Luis Augusto Eijy Nagai (@eijynagai)

# write the directories to be searched. Double-check the names!
projects=('covid19_project' \
          'signed_graphs_project' \
          'Hayashi_P10 IHEC_project' \
          'Nakajima_drosophila_velocito' \
          'nakajima_drosophila_project_v2' \
          'ChIP-seq_imputation' \
          'nakajima_drosophila_ageing')

# create a new directory to store all files. It uses the system date
mkdir -p $(date +"%Y-%m-%d")

cd $(date +"%Y-%m-%d")

for project in "${projects[@]}"; do

    # Dedicated directory for the project
    mkdir -p ${project}

    echo "Searching on ${project} directory..."
    # Download scrips in Jupyter, bash and python
    # -a: "archive", sync recursively, synlinks, permissions
    # -m: prune empty dirs
    # -v: verbose
    rsync -amv --progress --ignore-missing-args \
    --include='*.ipynb' \
    --include='*.sh' \
    --include='*.py' \
    --include='*/' \
    --exclude='*' \
    nagai@euphonium:/home/nagai/${project} .
    echo "Done!"
done
echo "Finished!"
