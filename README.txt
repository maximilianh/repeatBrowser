Repeat Browser consensus/lift building script
=============================================

requirements:
 - biopython
 - numpy
 - matplotlib (optional)
 - "muscle" aligner, in your PATH
 - bedtools
 
1) Adapt the config.py file to your environment. 

2) download Dfam.embl.gz to data/:
   wget https://www.dfam.org/releases/Dfam_3.1/families/Dfam.embl.gz

3) download hg38.fa, hg38.2bit and chrom.sizes (or your underlying genome) into data/
   from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/

4) run: python buildSeqs all
