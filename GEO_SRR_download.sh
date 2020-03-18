# -*- coding: utf-8 -*-

#bash GEO_SRR_download.sh input.txt

# Format of input file has to be sampleid SRR  
while read var1 var2
do
    echo Downloading $var1 $var2
        prefetch --max-size 100G $var2 --output-directory .
done < "$input"
