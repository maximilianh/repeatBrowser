Repeat Browser consensus/lift building script
=============================================

requirements:
 - biopython
 - numpy
 - matplotlib (optional)
 - "muscle" aligner, in your PATH
 - bedtools
 
1) Adapt the config.py file to your environment, e.g. set the genome identifier, like hg19 or hg38, or mm10

2) download Dfam.embl.gz to data/:
   wget https://www.dfam.org/releases/Dfam_3.1/families/Dfam.embl.gz -O Dfam.embl.gz

3) download hg38.fa, hg38.2bit and chrom.sizes (or your underlying genome identifier) into data/
   from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/

4) run: python buildSeqs all

   The step "rmsk" will download the hg38 or hg19 DFAM rmsk tables from https://hgwdev.gi.ucsc.edu/~max/repBrowser/

   If you use another genome, you may need to adapt the step and use an rmsk table from hgdownload, e.g.
   https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz

