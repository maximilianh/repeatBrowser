#!/bin/bash

#############README
#This file documents the commands used to generate all RepeatBrowser tracks
#Commands are intended to be run from the data/ directory of the track hub although it moves final output to subdirectories of the appropriate genome.
#This final file was not run as a single shell script, rather each section was run individually as we added tracks.  However this file serves as documentation
#for the commands that were run to generate Repeat Browser tracks. Theoretically, it should work as a single script assuming similar directory structure.
#
#If all you want to do is mirror the Repeat Browser we recommend simply downloading the processed files using rsync.
#If you want to generate your own tracks you can follow the operations here as general reference, but will likely need to adapt them to the formatting of your
#own genomic track data but you may notice a few path inconsistencies or occasional use of non-standard scripts.
###




#########Consensus Alignment track
#	Just a simple bar spanning the length of each consensus but with info on Dfam mappings etc. Also used for searchTrix.
#
#Fields as defined by data/build.as
#table build
#"repeat brower dfam build descriptions"
#    (
#    string chrom;        "Chromosome (or contig, scaffold, etc.)"
#    uint   chromStart;   "Start position in chromosome"
#    uint   chromEnd;     "End position in chromosome"
#    string name;         "Name of item"
#    uint   score;      "Score from 0-1000"
#    char[1] strand;    "+ or -"
#    uint thickStart;   "Start of where display should be thick (start codon)"
#    uint thickEnd;     "End of where display should be thick (stop codon)"
#    uint reserved;     "Used as itemRgb as of 2004-11-22"
#    int blockCount;    "Number of blocks"
#    int[blockCount] blockSizes; "Comma separated list of block sizes"
#    int[blockCount] chromStarts; "Start positions relative to chromStart"
#    string  rbName;            "Repeat Browser Name"
#    string  dfamName;            "DFAM Record"
#    string  dfamType;            "DFAM Type"
#    string  dfamClass;            "DFAM Class"
#    string  dfamOrg;            "DFAM Organismn"
#    lstring  Comment;            "DFAM Comment"
#    string  rmName;            "RepeatMasker Names"
#    string  repCount;            "Number of copies in genome, total"
#    string  splitRepCount;            "Copies in genome, per RepeatMasker Name"
#    string  liftCount;            "Liftable copies"
#    string  seqSource;            "Sequence Source(s)"
#    string  targetLen;            "Length of 99% percentile in genome"
#    string  maxLen;            "Maximum length in genome"
#    string  alnLen;            "Length of top50 alignment"
#    string  nCount;            "Number of Ns in top50 consensus"
#    string  ratio;            "Ratio of alignment length to top50 length"
#    )
#
#	Output: RepeatConsensusAlignment.bb
###

#build.bb is an output of the buildSeqs step run by Max
#I originally took it from Max's hgwdev directory while we were dividing & conquering but a complete run of buildTracks should result in build.bb in the tracks/subdirectory
#wget https://hgwdev.gi.ucsc.edu/~max/kznf/hg38reps/tracks/build.bb 
#searchTrix need to be set as well if you want to make hub searchable: https://genome.ucsc.edu/goldenPath/help/trix.html

mv tracks/build.bb .
bigBedToBed build.bb build.bed



#########ORFs track
# 	ORFs: Using Max's summary.tsv
#
#	There are three classes of repeat browser sequences, 1) those that have exact matches to the Dfam name, 2) Those that have "int" stripped in the Dfam name 3) those that were built by buildSeqs
#
#	For exact matches, search Dfam.embl for the name using the NM record in Dfam embl and grab the whole record (ending in \\). Then grep CDS features and print tab separated field with <cons name> <start stop> <ORF name>
#	For int-strip, do the same thing.
#	For everything else (that isn't a manual consensus) do the same thing.
#	So in the end I could have just used the same command for all three cases (awk '($6 != cons), but I had to think it through and assumed cases would be different going in. This is exactly what was run.
#	This would have also been much faster splitting Dfam.embl into the appropriate records rather than searching the whole file over and over again, but not a major issue for this particular use case.
#
#	For the built consensuses I first make a fasta file of built consensuses, and make a blast database/
#	Then I subset the RepeatMasker peptide library to those species/phylogenetic groupings for TEs found in hg38
#	Then I use tblastn to find the highest scoring gag and pol for each built consensus (since built consensuses are all LINEs)
#	The tblastn bitscore must also be above 50. Because of bed numbering vs tblastn I shifted start down one (since I noticed ATG becoming TG)
#
#	Output: orfs.bb
###

#for exact matches 
awk '($6=="exact")' summary.tsv  | cut -f 1 | while read i; do sed -n -e '/^NM\s*'${i}'$/,/\/\// p' Dfam.embl | grep -A1 "FT.*CDS "| tr -s "\n" "\t" | sed "s/--/\n/g" | tr -s "." " \t" | tr -s "  " "\t" | tr -s "\"" "\t" | awk -v s=$i '{print s"\t"$3"\t"$4"\t"$7}'; done > exact_annotations.tab

#for renamed matches
awk '($6=="int-strip")' summary.tsv  | cut -f 1 | while read i; do sed -n -e '/^NM\s*'${i}'$/,/\/\// p' Dfam.embl | grep -A1 "FT.*CDS "| tr -s "\n" "\t" | sed "s/--/\n/g" | tr -s "." " \t" | tr -s "  " "\t" | tr -s "\"" "\t" | awk -v s=$i '{print s"\t"$3"\t"$4"\t"$7}'; done > intstrip_annotations.tab

#for manual matches
awk '($6!="int-strip" && $6!="exact" && $6!="cons")' summary.tsv  |  cut -f 1 | while read i; do sed -n -e '/^NM\s*'${i}'$/,/\/\// p' Dfam.embl | grep -A1 "FT.*CDS "| tr -s "\n" "\t" | sed "s/--/\n/g" | tr -s "." " \t" | tr -s "  " "\t" | tr -s "\"" "\t" | awk -v s=$i '{print s"\t"$3"\t"$4"\t"$7}'; done > manual_annotations.tab


#Dfam annotations has the annotations from Dfam.embl for sequences that are in the Repeat Browser
cat intstrip_annotations.bed manual_annotations.bed exact_annotations.bed | awk '{print $0"\t0\t+"}'| sort | uniq > Dfam_annotations.bed


#built cons
awk '($6=="cons")' ../summary.tsv  |  cut -f 1 | while read i; do sed -n -e '/^>'${i}'$/,/^>/ p' ~/hive/jferna10/RepeatBrowserHub/hg38reps/hg38reps.fa | head -n-1 >> built_cons.fa; done
cd peps/

#RepeatMasker pep library from Hubley
wget https://github.com/rmhubley/RepeatMasker/raw/master/Libraries/RepeatPeps.lib

#make a blastdb for tblastn
makeblastdb -in built_cons.fa -dbtype nucl -parse_seqids

#phylogenetic groupings for TEs found in hg38 (this list was arrived at by taking the Dfam name, reduce repeat peptide libs to peptides from those groups
cut -f 1 summary.tsv | while read i ; do sed -n -e '/^NM\s*'${i}'$/,/\/\// p' Dfam.embl | grep "^OS" | cut -f2; done >  hg38_TE_phlyo.txt
cat hg38_TE_phlyo.txt | tr -s " " "\t" | cut -f2 |  sort | uniq | tr "\n" "|"

#copy and paste output to grep,reduces peptide library significantly, prevents aligning proteins from random creatures to human
grep -E  'Afrotheria|Amniota|Boreoeutheria|Carnivora|Catarrhini|Cetartiodactyla|Euarchontoglires|Euteleostomi|Eutheria|Glires|Haplorrhini|Hominidae|Homininae|Hominoidea|Homo|Mammalia|Metatheria|Metazoa|Monotremata|Primates|Rodentia|Sauropsida|Scandentia|Simiiformes|Tetrapoda|Theria|Vertebrata|Xenarthra' RepeatPeps.lib -A1 > HumanRepeatPeps.lib 
grep ">" built_cons.fa | cut -f 2 -d">" | while read i; do grep -A1 "$i" HumanRepeatPeps.lib;done > Human_built_cons_peps.lib
sed -zi "s/--\n//g" Human_built_cons_peps.lib

#find peptides in consensuses
tblastn -db built_cons.fa -query Human_built_cons_peps.lib -outfmt 6 > built_consensus_annotations.tab

cat built_consensus_annotations.tab | sort | uniq | sort -k 2,12 -n -r > temp.tab

#take the highest scoring gag and pol, (since all built consesnsus are L1s); Also impose filter that peptide score must also be above 50
grep ">" built_cons.fa | cut -f 2 -d">" | while read i; do grep -P "\t$i\t" temp.tab | grep "gag" | sort -k 12 -n -r | head -n+1 |  awk '($12> 50) {print $2"\t"$9-1"\t"$10"\t"$1"\t"$12"\t+"}';done | sort | uniq > cons_gag.bed
grep ">" built_cons.fa | cut -f 2 -d">" | while read i; do grep -P "\t$i\t" temp.tab | grep "pol" | sort -k 12 -n -r | head -n+1 |  awk '($12> 50) {print $2"\t"$9-1"\t"$10"\t"$1"\t"$12"\t+"}';done | sort | uniq > cons_pol.bed

mv *.bed ../
cd ..

#combine and change score to 0 (avoid issues with decimals etc that make browser unhappy)
cat Dfam_annotations.bed cons_gag.bed cons_pol.bed HERV_full.bed| awk '{print $1"\t"$2"\t"$3"\t"$4"\t"0"\t"$6}' | sed "s/L1PBA1/L1PBa1/g" | bedtools sort > ORFs.bed
bedToBigBed ORFs.bed ../hg38reps/hg38reps.sizes ORFs.bb
mv ORFs.bb ../hg38reps/orfs.bb


#########knownGene32
#	liftOver of the genes to hg38reps 
#
#	Get knownGene table from UCSC
#	However liftOver struggles with the bed12 format; the only way to lift the exons of genes is split to a bed of the exons.
#	Splitting to a bed of exons -> loss of coding/non-coding, therefore first split to list of coding exons, ncRNA exons, and UTR exons
#	Then lift all three lists to hg38reps.
#	Along the way join back all the desired identifiers.
#
#	Output: 	gencode_cds.bb
#				gencode_ncRNA.bb 
#				gencode_utr.bb
###

#get known geneTable
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/knownGene32.bb
bigBedToBed knownGene32.bb knownGene32.bed
cut -f 1-12 knownGene32.bed > knownGene32_12.bed 

#sub out ENST with chr coords, make bed12
awk '$4=$1":"$2"-"$3' knownGene32.bed | tr -s " " "\t" | cut -f 1-12 > knownGene32_12.bed 

#sub out score with chr coords, build key file, this is used to join all IDs back together later
awk '$5=$1":"$2"-"$3' knownGene32.bed  | tr -s " " "\t" | cut -f 4-8,10,17  > knownGene32_keyfile.txt

#get list of exons
bedToExons knownGene32_12.bed  knownGene32_exons.bed

#if thickstart = thickstop then this gene has no coding at all (e.g. ncRNA) 
awk '$7 == $8' knownGene32_12.bed > ncRNA_knownGene32_temp.bed

#remove the ncRNA from the knownGeneFile to get exons from protein coding genes
#uniqueToFileOne is just a onliner in my bin: diff -U $(wc -l < $1) $1 $2 | sed -n 's/^-//p'
uniqueToFileOne.sh knownGene32_12.bed ncRNA_knownGene32_temp.bed  | tail -n+2 > cds_knownGene32_12.bed

#exon list for ncRNA
bedToExons ncRNA_knownGene32_temp.bed ncRNA_knownGene32_exons.bed 

#split protein coding genes to cds 
bedToExons -cdsOnly cds_knownGene32_12.bed temp.bed 
bedToExons	cds_knownGene32_12.bed knownGene32_exons.bed

#remove exons where start = stop
awk '$2 != $3' temp.bed > cds_knownGene32_exons.bed

#some exons have UTR+coding on them, split them in two 
bedtools subtract -a knownGene32_exons.bed -b cds_knownGene32_exons.bed > utr_knownGene32_exons.bed

#lift all 3 groups
liftOver ncRNA_knownGene32_exons.bed ../lift/hg38_to_hg38reps.over.chain  RepBro_ncRNA_knownGene32_12.bed RepBro_ncRNA_knownGene32_12.unmapped -multiple
liftOver cds_knownGene32_exons.bed ../lift/hg38_to_hg38reps.over.chain  RepBro_cds_knownGene32_exons.bed RepBro_cds_knownGene32_exons.unmapped -multiple
liftOver utr_knownGene32_exons.bed ../lift/hg38_to_hg38reps.over.chain  RepBro_utr_knownGene32_exons.bed RepBro_utr_knownGene32_exons.unmapped -multiple

#add thickStart and thickStop back depending on coding/non-coding
awk -v OFS="\t" '{$7=$2;$8=$2;$9=0; print $0}' RepBro_ncRNA_knownGene32_12.bed > t1.bed
awk -v OFS="\t" '{$7=$2;$8=$3;$9=0; print $0}' RepBro_cds_knownGene32_exons.bed > t2.bed
awk -v OFS="\t" '{$7=$2;$8=$2;$9=0; print $0}' RepBro_utr_knownGene32_exons.bed > t3.bed

#add all IDs back
join -1 4 -2 1 <(sort -k4 RepBro_ncRNA_knownGene32_12.bed) <(sort -k1 knownGene32_keyfile.txt | cut -f 1,2,7) | awk -v OFS="\t" '{print $2,$3,$4,$8,"0", $6,$3,$3,$1,$7}' | sort | uniq | bedtools sort > RepBro_ncRNA_knownGene32_14.bed
join -1 4 -2 1 <(sort -k4 RepBro_cds_knownGene32_exons.bed) <(sort -k1 knownGene32_keyfile.txt | cut -f 1,2,7) | awk -v OFS="\t" '{print $2,$3,$4,$8,"0", $6,$3,$4,$1,$7}'| sort | uniq | bedtools sort > RepBro_cds_knownGene32_exons_14.bed
join -1 4 -2 1 <(sort -k4 RepBro_utr_knownGene32_exons.bed) <(sort -k1 knownGene32_keyfile.txt | cut -f 1,2,7) | awk -v OFS="\t" '{print $2,$3,$4,$8,"0", $6,$3,$3,$1,$7}' | sort | uniq | bedtools sort > RepBro_utr_knownGene32_exons_14.bed

#collapse alternate transcripts from the same gene name that span the same region to one record
bedtools groupby -g 1,2,3,4,5,6,7,8 -o collapse -c 9,10 -i RepBro_cds_knownGene32_exons_14.bed > collapse_RepBro_cds_knownGene32_exons_14.bed
bedtools groupby -g 1,2,3,4,5,6,7,8 -o collapse -c 9,10 -i RepBro_ncRNA_knownGene32_14.bed | awk 'BEGIN{OFS=FS="\t"} {$10=substr($10,1,250); $9=substr($9,1,250)}1' > collapse_RepBro_ncRNA_knownGene32_14.bed
bedtools groupby -g 1,2,3,4,5,6,7,8 -o collapse -c 9,10 -i RepBro_utr_knownGene32_exons_14.bed | awk  'BEGIN{OFS=FS="\t"} {$10=substr($10,1,250); $9=substr($9,1,250)}1' > collapse_RepBro_utr_knownGene32_exons_14.bed

#convert to bigBed
bedToBigBed collapse_RepBro_cds_knownGene32_exons_14.bed ../hg38reps/hg38reps.sizes gencode_cds.bb -type=bed8+2 -as=gencode.as
bedToBigBed collapse_RepBro_ncRNA_knownGene32_14.bed ../hg38reps/hg38reps.sizes gencode_ncRNA.bb -type=bed8+2 -as=gencode.as
bedToBigBed collapse_RepBro_utr_knownGene32_exons_14.bed ../hg38reps/hg38reps.sizes gencode_utr.bb -type=bed8+2 -as=gencode.as


#########schmittges2016
#	liftOver of the Schmittges ZNF ChIP-SEQ to hg38reps 
#	
#	Copied bb from hg19, swap chromosome coordinates to name
#	Output: 160 *_hg38reps.bb files in ../hg38reps/schmittges2016
#	Output: 160 *_hg38reps.bw files in ../hg38reps/schmittges2016
###

#make bb summits
mkdir ../hg38reps/schmittges2016
ls ../hg19/schmittges2016/*.bb | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do bigBedToBed ../hg19/schmittges2016/${i}.bb temp.bed; awk '$4=$1":"$2"-"$3' temp.bed > test.bed; liftOver -multiple test.bed ../lift/hg19_to_hg38reps.over.chain ../hg38reps/schmittges2016/${i}_hg38reps.bed ../hg38reps/schmittges2016/${i}_hg38reps.unmapped; bedSort ../hg38reps/schmittges2016/${i}_hg38reps.bed ../hg38reps/schmittges2016/${i}_hg38reps.bed; bedToBigBed ../hg38reps/schmittges2016/${i}_hg38reps.bed ../hg38reps/hg38reps.sizes ../hg38reps/schmittges2016/${i}_hg38reps.bb; done

#make coverage
ls ../hg38reps/schmittges2016/*.bed | cut -f 4 -d"/" | cut -f 1 -d"." | while read i; do bedtools  genomecov -bg -split -i ../hg38reps/schmittges2016/${i}.bed -g ../hg38reps/hg38reps.sizes > temp.bg; bedGraphToBigWig temp.bg ../hg38reps/hg38reps.sizes ../hg38reps/schmittges2016/${i}.bw; done


#make trackdb
echo -e "track schmittges2016\nshortLabel Schmittges_Hughes 2016\nlongLabel Schmittges_Hughes Zinc Finger ChIP-Seq\ngroup summitCov\ntype bigWig\ncompositeTrack on\nvisibility hide\n" > trackDb.schmittges2016.txt

ls ../hg38reps/schmittges2016/*.bw | cut -f 4 -d"/" | cut -f 1 -d"."| awk '{split($1,t,"_"); print "\n\ttrack "$1"_cov\n\tparent schmittges2016\n\tshortLabel "t[1]" "t[2]"\n\tlongLabel schmittges2016 "t[1]" "t[2]"\n\tvisibility dense\n\ttype bigWig\n\tbigDataUrl schmittges2016/"$1".bw" }' >> trackDb.schmittges2016.txt 

echo -e "\ntrack schmittges2016Summits\nshortLabel Schmittges_Hughes 2016\nlongLabel Schmittges_Hughes Zinc Finger ChIP-Seq\ngroup summits\ntype bigBed\ncompositeTrack on\nvisibility hide\n" >> trackDb.schmittges2016.txt

ls ../hg38reps/schmittges2016/*.bb | cut -f 4 -d"/" | cut -f 1 -d"."| awk '{split($1,t,"_"); print "\n\ttrack "$1"_sum\n\tparent schmittges2016Summits\n\tshortLabel "t[1]" "t[2]"\n\tlongLabel schmittges2016 "t[1]" "t[2]"\n\tvisibility dense\n\ttype bigBed\n\tbigDataUrl schmittges2016/"$1".bb" }' >> trackDb.schmittges2016.txt 

mv trackDb.schmittges2016.txt ../hg38reps/


#########Imbeault_Trono 2017
#	liftOver of the Imbeault KZNF ChIP-SEQ data to hg38reps 
#	
#
#	Output: 251 *_hg38reps.bb files in ../hg38reps/imbeault2017
###


#make bb summits
mkdir ../hg38reps/imbeault2017
ls ../hg19/imbeault2017/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do awk '$4=$1":"$2"-"$3' ../hg19/imbeault2017/*${i}.bed > test.bed; liftOver -multiple test.bed ../lift/hg19_to_hg38reps.over.chain ../hg38reps/imbeault2017/${i}_hg38reps.bed ../hg38reps/imbeault2017/${i}_hg38reps.unmapped; bedSort ../hg38reps/imbeault2017/${i}_hg38reps.bed ../hg38reps/imbeault2017/${i}_hg38reps.bed; bedToBigBed ../hg38reps/imbeault2017/${i}_hg38reps.bed ../hg38reps/hg38reps.sizes ../hg38reps/imbeault2017/${i}_hg38reps.bb; done

#make coverage
ls ../hg38reps/imbeault2017/*.bed | cut -f 4 -d"/" | cut -f 1 -d"." | while read i; do bedtools  genomecov -bg -split -i ../hg38reps/imbeault2017/${i}.bed -g ../hg38reps/hg38reps.sizes > temp.bg; bedGraphToBigWig temp.bg ../hg38reps/hg38reps.sizes ../hg38reps/imbeault2017/${i}.bw; done


#make trackdb
echo -e "track imbeault2017\nshortLabel Imbeault_Trono 2017\nlongLabel Imbeault_Trono Zinc Finger ChIP-Seq\ngroup summitCov\ntype bigWig\ncompositeTrack on\nvisibility hide\n" > trackDb.imbeault2017.txt

ls ../hg38reps/imbeault2017/*.bw | cut -f 4 -d"/" | cut -f 1 -d"."| awk '{split($1,t,"_"); print "\n\ttrack "$1"_cov\n\tparent imbeault2017\n\tshortLabel "t[1]" "t[2]"\n\tlongLabel imbeault2017 "t[1]" "t[2]"\n\tvisibility dense\n\ttype bigWig\n\tbigDataUrl imbeault2017/"$1".bw" }' >> trackDb.imbeault2017.txt

echo -e "\ntrack imbeault2017Summits\nshortLabel Imbeault_Trono 2017\nlongLabel Imbeault_Trono Zinc Finger ChIP-Seq\ngroup summits\ntype bigBed\ncompositeTrack on\nvisibility hide\n" >> trackDb.imbeault2017.txt

ls ../hg38reps/imbeault2017/*.bb | cut -f 4 -d"/" | cut -f 1 -d"."| awk '{split($1,t,"_"); print "\n\ttrack "$1"_sum\n\tparent imbeault2017Summits\n\tshortLabel "t[1]" "t[2]"\n\tlongLabel imbeault2017 "t[1]" "t[2]"\n\tvisibility dense\n\ttype bigBed\n\tbigDataUrl imbeault2017/"$1".bb" }' >> trackDb.imbeault2017.txt

mv trackDb.imbeault2017.txt ../hg38reps/




#########trf
#	calculation of tandem repeats on hg38reps using trf
#	
#
#	Output: trf.bb
###

trfBig ../hg38reps/hg38reps.fa trf.out bedAt=trf.raw
cat trf.tab | gawk 'BEGIN {FS="\t"; OFS="\t"} {$4=$16; if (length($4) > 256){$4=substr($4,1,206)"_truncated"}; NR=14; print}' | cut -f-15 > trf.bed
bedSort trf.bed trf.bed
bedToBigBed trf.bed ../hg38reps/hg38reps.sizes trf.bb -type=bed4+11 #-as=



#########other_cons_aln
#	self alignments of all repeat consensuses to one another
#	
#
#	Output: other_cons_aln.bb
###

blat ../hg38reps/hg38reps.fa ../hg38reps/hg38reps.fa other_cons_aln.psl
pslToBigPsl other_cons_aln.psl stdout | sort -k1,1 -k2,2n > other_cons_aln.txt
wget https://genome.ucsc.edu/goldenPath/help/examples/bigPsl.as
bedSort other_cons_aln.txt other_cons_aln.txt
bedToBigBed as=bigPsl.as -type=bed12+13 -tab other_cons_aln.txt ../hg38reps/hg38reps.sizes other_cons_aln.bb


#########Tsankov_Meissner2014
#	liftOver of the Tsankov ChIP-SEQ to hg38reps 
#	
#
#	Output: 204 *_hg38reps.bb files in ../hg38reps/imbeault2017
###
mkdir ../hg19/tsankov2014
mkdir ../hg38reps/tsankov2014

grep -B3 "longLabel" trackDb.tsankov2014.txt | tr "\n" "\t" | tr -s "-" "\n" | cut -f 2,5 | cut -f 5,10-16 -d" " | tr -s "_" "\t" | cut -f 1,3 | sed "s/\//-/g" | tr -s " " "_" |  sed "s/\t_/\t/g" > tsankov2014.key
cat tsankov2014.key | while read i g; do mv ../hg19/tsankov2014/${i}_peaks.bb ../hg19/tsankov2014/${i}_${g}_peaks.bb; done

#make summits
ls ../hg19/tsankov2014/*.bb | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do bigBedToBed ../hg19/tsankov2014/${i}.bb temp.bed; awk '$4=$1":"$2"-"$3' temp.bed > test.bed; liftOver -multiple test.bed ../lift/hg19_to_hg38reps.over.chain ../hg38reps/tsankov2014/${i}_hg38reps.bed ../hg38reps/tsankov2014/${i}_hg38reps.unmapped; bedSort ../hg38reps/tsankov2014/${i}_hg38reps.bed ../hg38reps/tsankov2014/${i}_hg38reps.bed; bedToBigBed ../hg38reps/tsankov2014/${i}_hg38reps.bed ../hg38reps/hg38reps.sizes ../hg38reps/tsankov2014/${i}_hg38reps.bb; done

#make coverage
ls ../hg38reps/tsankov2014/*.bed | cut -f 4 -d"/" | cut -f 1 -d"." | while read i; do bedtools  genomecov -bg -split -i ../hg38reps/tsankov2014/${i}.bed -g ../hg38reps/hg38reps.sizes > temp.bg; bedGraphToBigWig temp.bg ../hg38reps/hg38reps.sizes ../hg38reps/tsankov2014/${i}.bw; done


#make trackdb
echo -e "track tsankov2014\nshortLabel Tsankov_Meissner 2014\nlongLabel Tsankov_Meissner TF ChIP-Seq in ES differentiation\ngroup summitCov\ntype bigWig\ncompositeTrack on\nvisibility hide\n" > trackDb.tsankov2014.txt

ls ../hg38reps/tsankov2014/*.bw | cut -f 4 -d"/" | cut -f 1 -d"."| awk '{split($1,t,"_"); if (t[4]!="peaks") {t[3] ="$t[3] $t[4]"}; print "\n\ttrack "$1"_cov\n\tparent tsankov2014\n\tshortLabel "t[2]" "t[3]"\n\tlongLabel tsankov2014 "t[1]" "t[2]" "t[3]"\n\tvisibility dense\n\ttype bigWig\n\tbigDataUrl tsankov2014/"$1".bw" }' >> trackDb.tsankov2014.txt 

echo -e "\ntrack tsankov2014Summits\nshortLabel Tsankov_Meissner 2014\nlongLabel Tsankov_Meissner TF ChIP-Seq in ES differentiation\ngroup summits\ntype bigBed\ncompositeTrack on\nvisibility hide\n" >> trackDb.tsankov2014.txt

ls ../hg38reps/tsankov2014/*.bb | cut -f 4 -d"/" | cut -f 1 -d"."| awk '{split($1,t,"_"); if (t[4]!="peaks") {t[3] ="$t[3] $t[4]"} ; print "\n\ttrack "$1"_sum\n\tparent tsankov2014Summits\n\tshortLabel "t[2]" "t[3]"\n\tlongLabel tsankov2014 "t[1]" "t[2]" "t[3]"\n\tvisibility dense\n\ttype bigBed\n\tbigDataUrl tsankov2014/"$1".bb" }' >> trackDb.tsankov2014.txt 

mv trackDb.tsankov2014.txt ../hg38reps/



#########Mappings to the Human Genome
#	liftOver of up to 500 repeat instance alignments to the consensus  
#	Taken directly from Max's output from the buildLiftOver
#	I just took the psls and downsamples to 500 instances and converted to bigPsl
#
#	Output: mappings_500.bb
###

wget -r –level=0 -E –ignore-length -x -k -p -erobots=off -np -N https://hgwdev.gi.ucsc.edu/~max/kznf/hg38reps/hg38/seqs/pslLift/

ls pslLift/*.psl | while read i; do shuf $i | head -n 500; done >  merged_500.psl 


awk '$10=($10":"$11"-"$12)' merged_500.psl | tr -s " " "\t" > more_500.psl


pslToBigPsl more_500.psl stdout | sort -k1,1 -k2,2n > more_500.txt
wget https://genome.ucsc.edu/goldenPath/help/examples/bigPsl.as
bedSort more_500.txt more_500.txt
bedToBigBed as=bigPsl.as -type=bed12+13 -tab more_500.txt ../hg38reps/hg38reps.sizes mappings_500.bb



#########ENCODE TFBS
#	liftOver of ENCODE TFBS collection from many cell types
#	These files required a little more parsing to keep Cell Type and factor information consistent.
#
#	Output: .bb and bw files
###
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz
gunzip wgEncodeRegDnaseClusteredV3.bed.gz

sed "s/\//-/g" encRegTfbsClusteredWithCells.hg38.bed  | awk '{split ($6,a,","); for (i in a) { print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t+") >> "../hg19/encRegTfbsClustered/encRegTfbsClusteredWithCells_"$4"_"a[i]".bed"}}' 


#make bb summits
mkdir ../hg38reps/encRegTfbsClustered
ls ../hg38/encRegTfbsClustered/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do awk '$4=$1":"$2"-"$3' ../hg38/encRegTfbsClustered/*${i}.bed > test.bed; liftOver -multiple test.bed ../lift/hg38_to_hg38reps.over.chain ../hg38reps/encRegTfbsClustered/${i}_hg38reps.bed ../hg38reps/encRegTfbsClustered/${i}_hg38reps.unmapped; bedSort ../hg38reps/encRegTfbsClustered/${i}_hg38reps.bed ../hg38reps/encRegTfbsClustered/${i}_hg38reps.bed; bedToBigBed ../hg38reps/encRegTfbsClustered/${i}_hg38reps.bed ../hg38reps/hg38reps.sizes ../hg38reps/encRegTfbsClustered/${i}_hg38reps.bb; done

#make coverage
ls ../hg38reps/encRegTfbsClustered/*.bed | cut -f 4 -d"/" | cut -f 1 -d"." | while read i; do bedtools  genomecov -bg -split -i ../hg38reps/encRegTfbsClustered/${i}.bed -g ../hg38reps/hg38reps.sizes > temp.bg; bedGraphToBigWig temp.bg ../hg38reps/hg38reps.sizes ../hg38reps/encRegTfbsClustered/${i}.bw; done


#make trackdb
echo -e "track ENCODE_TFBS\nshortLabel ENCODE_TFBS\nlongLabel ENCODE TFBS Clustered\ngroup summitCov\ntype bigWig\ncompositeTrack on\nvisibility hide\n" > trackDb.encTFBSclustered.txt

ls ../hg38reps/encRegTfbsClustered/*.bw | cut -f 4 -d"/" | cut -f 1 -d"."| awk '{split($1,t,"_"); print "\n\ttrack "$1"_cov\n\tparent ENCODE_TFBS\n\tshortLabel "t[2]" "t[3]"\n\tlongLabel ENCODE TFBS Clustered "t[2]" "t[3]"\n\tvisibility dense\n\ttype bigWig\n\tbigDataUrl encRegTfbsClustered/"$1".bw" }' >> trackDb.encTFBSclustered.txt

echo -e "\ntrack ENCODE_TFBSSummits\nENCODE_TFBS 2017\nlongLabel Imbeault_Trono Zinc Finger ChIP-Seq\ngroup summits\ntype bigBed\ncompositeTrack on\nvisibility hide\n" >> trackDb.encTFBSclustered.txt

ls ../hg38reps/encRegTfbsClustered/*.bb | cut -f 4 -d"/" | cut -f 1 -d"."| awk '{split($1,t,"_"); print "\n\ttrack "$1"_sum\n\tparent ENCODE_TFBS_Summits\n\tshortLabel "t[2]" "t[3]"\n\tlongLabel iENCODE TFBS Clustered "t[2]" "t[3]"\n\tvisibility dense\n\ttype bigBed\n\tbigDataUrl encRegTfbsClustered/"$1".bb" }' >> trackDb.encTFBSclustered.txt

mv trackDb.encTFBSclustered.txt ../hg38reps/


######encode DNAse
# Liftover of ENCODE DNAse Hypersensitivity data.
# In order to make track manipulation easier (user can check/uncheck cell types and factors) I split the file up.
#
# Output: *.bw and *.bb for ENCODE DNAse
###

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseMasterSites/wgEncodeAwgDnaseMasterSites.bed.gz
gunzip *.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseMasterSites/wgEncodeAwgDnaseMasterSources.tab

#save bed coordinate and temporarily move cell types to name field
awk -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3"*"$8,$5,"+"}' wgEncodeAwgDnaseMasterSites.bed > encode_DNAse_5.bed

#lift peaks
liftOver -multiple encode_DNAse_5.bed ../lift/hg19_to_hg38reps.over.chain encode_DNAse_5_hg38reps.bed encode_DNAse_5_hg38reps.unmapped

#split the bed coord out from the cell type 
awk -v OFS="\t" '{split($4, a, "*"); {print $1,$2,$3,a[1],$5,"+",a[2]}}' encode_DNAse_5_hg38reps.bed > temp.bed

#for every cell type make a new entry 
awk -v OFS="\t" '{split($7, a, ","); for(i in a) { if (a[i] ~/^[0-9]+$/ ){print $1,$2,$3,$4,$5,"+",a[i]}}}' temp.bed | sort -k 7b,7 > encode_DNAse_7_hg38reps.bed

#add the cell names for the numbers
sort -k1b,1  wgEncodeAwgDnaseMasterSources.tab > DNAse_key.tab

join -1 1 -2 7  DNAse_key.tab encode_DNAse_7_hg38reps.bed | awk -v OFS="\t" '{print $3,$4,$5,$6,$7,$8, $2}' >  encode_DNAse_7_hg38reps_named.bed

mv encode_DNAse_7_hg38reps_named.bed ~/hive/jferna10/RepeatBrowserHub/hg38reps/wgEncodeAwgDnaseMasterSites

awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6  >> "encode_DNAse_6_" $7 ".bed"}' encode_DNAse_7_hg38reps_named.bed
mv encode_DNAse_7_hg38reps_named.bed encode_DNAse_7_hg38reps_named.bak

ls *.bed | cut -f 1 -d"."|  while read i; do bedSort $i.bed $i.bed ; bedtools genomecov -bg -split -i $i.bed -g ../hg38reps.sizes > temp.bg; bedGraphToBigWig temp.bg ../hg38reps.sizes $i.bw; done
ls *.bed | cut -f 1 -d"."|  while read i; do bedSort $i.bed $i.bed ; bedToBigBed $i.bed ../hg38reps.sizes $i.bb -type=bed6; done



#make trackdb
echo -e "track ENCODE_DNAse\nshortLabel ENCODE DNAse Summits\nlongLabel ENCODE DNAse Hypersensitivity (125 cell types) v3 \ngroup histone\ntype bigBed\ncompositeTrack on\nvisibility hide\n" > trackDb.dnase.txt

ls *.bw | cut -f 1 -d"." |  awk '{split($1, t, "_")}{print "\n\ttrack "$1"_cov\n\tparent ENCODE_DNAse\n\tshortLabel "t[4]" "t[5]"\n\tlongLabel ENCODE_DNAse "t[4]" "t[5]"\n\tvisibility hide\n\ttype bigWig\n\tbigDataUrl wgEncodeAwgDnaseMasterSites/"$1".bw" }' >> trackDb.dnase.txt

echo -e "track ENCODE_DNAse_Coverage\nshortLabel ENCODE DNAse Coverage \nlongLabel ENCODE DNAse Hypersensitivity (125 cell types) v3 \ngroup histone\ntype bigWig\ncompositeTrack on\nvisibility hide\n" >> trackDb.dnase.txt

ls *.bb | cut -f 4 -d"/" | cut -f 1 -d"." | awk '{split($1, t, "_")}{print "\n\ttrack "$1"_sum\n\tparent ENCODE_DNAse_Coverage\n\tshortLabel "t[4]" "t[5]"\n\tlongLabel ENCODE_DNAse_Coverage "t[4]" "t[5]"\n\tvisibility hide\n\ttype bigBed\n\tbigDataUrl wgEncodeAwgDnaseMasterSites/"$1".bb" }' >> trackDb.dnase.txt



#####Theunissen
# Liftover of Primed-Naive ChIP-SEQ
#
# Output: *.bw and *.bb for Theunissen collection 
###
ls ../hg19/theunissen2016/*.bb | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do bigBedToBed ../hg19/theunissen2016/${i}.bb temp.bed; awk '$4=$1":"$2"-"$3' temp.bed > test.bed; liftOver -multiple test.bed ../lift/hg19_to_hg38reps.over.chain ../hg38reps/theunissen2016/${i}_hg38reps.bed ../hg38reps/theunissen2016/${i}_hg38reps.unmapped; bedSort ../hg38reps/theunissen2016/${i}_hg38reps.bed ../hg38reps/theunissen2016/${i}_hg38reps.bed; bedToBigBed ../hg38reps/theunissen2016/${i}_hg38reps.bed ../hg38reps/hg38reps.sizes ../hg38reps/theunissen2016/${i}_hg38reps.bb; done

#make coverage
ls ../hg38reps/theunissen2016/*.bed | cut -f 4 -d"/" | cut -f 1 -d"." | while read i; do bedtools  genomecov -bg -split -i ../hg38reps/theunissen2016/${i}.bed -g ../hg38reps/hg38reps.sizes > temp.bg; bedGraphToBigWig temp.bg ../hg38reps/hg38reps.sizes ../hg38reps/theunissen2016/${i}.bw; done


#make trackdb
echo -e "track theunissen2016\nshortLabel theunissen2016S\nlongLabel theunissen 2016 Primed Naive\ngroup summitCov\ntype bigWig\ncompositeTrack on\nvisibility hide\n" > trackDb.theunissen2016.txt

ls ../hg38reps/theunissen2016/*.bw | cut -f 4 -d"/" | cut -f 1 -d"."| awk '{split($1,t,"_"); print "\n\ttrack "$1"_cov\n\tparent theunissen2016\n\tshortLabel "t[2]" "t[3]"\n\tlongLabel theunissen2016 Primed Naive "t[2]" "t[3]"\n\tvisibility dense\n\ttype bigWig\n\tbigDataUrl theunissen2016/"$1".bw" }' >> trackDb.theunissen2016.txt

echo -e "\ntrack theunissen2016Summits\nshortLabel theunissen2016Summits\nlongLabel theunissen 2016 Primed Naive\ngroup summits\ntype bigBed\ncompositeTrack on\nvisibility hide\n" >> trackDb.theunissen2016.txt

ls ../hg38reps/theunissen2016/*.bb | cut -f 4 -d"/" | cut -f 1 -d"."| awk '{split($1,t,"_"); print "\n\ttrack "$1"_sum\n\tparent theunissen2016Summits\n\tshortLabel "t[2]" "t[3]"\n\tlongLabel theunissen2016 Primed Naive "t[2]" "t[3]"\n\tvisibility dense\n\ttype bigBed\n\tbigDataUrl theunissen2016/"$1".bb" }' >> trackDb.theunissen2016.txt

mv trackDb.theunissen2016.txt ../hg38reps/


############encode histone
# Liftover of ENCODE ChIP-SEQ for many cell types and many histone marks.
#
#	Output: *.bw  and *.bb for each combination
###
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/*broadPeak* ./
ls | while read i; do  awk -v OFS="\t" '{if ($8 > 10){print $1,$2,$3, "peak_"NR, $5, $6}}' $i > temp.bed; cat temp.bed > $i; done 


ls | while read i; do awk '$4=$1":"$2"-"$3' $i > test.bed; liftOver -multiple test.bed ../../lift/hg19_to_hg38reps.over.chain ../../hg38reps/wgEncodeUwHistone/${i}_hg38reps.bed ../../hg38reps/wgEncodeUwHistone/${i}_hg38reps.unmapped; bedSort ../../hg38reps/wgEncodeUwHistone/${i}_hg38reps.bed ../../hg38reps/wgEncodeUwHistone/${i}_hg38reps.bed; bedToBigBed ../../hg38reps/wgEncodeUwHistone/${i}_hg38reps.bed ../../hg38reps/hg38reps.sizes ../../hg38reps/wgEncodeUwHistone/${i}_hg38reps.bb; done


ls *.bed | cut -f 1 -d"." | while read i; do bedtools  genomecov -bg -split -i ${i}.bed -g ../../hg38reps/hg38reps.sizes > temp.bg; bedGraphToBigWig temp.bg ../../hg38reps/hg38reps.sizes ${i}.bw; done

#make trackdb
echo -e "track ENCODE_UW_Histone\nshortLabel ENCODE UW Histone\nlongLabel ENCODE UW Histone\ngroup histone\ntype bigWig\ncompositeTrack on\nvisibility hide\n" > trackDb.uw_histone.txt

ls ../hg38reps/wgEncodeUwHistone/*.bb | cut -f 4 -d"/" | cut -f 1 -d"."| sed 's/\([^[:blank:]]\)\([[:upper:]]\)/\1*\2/g' | awk '{split($1,t,"*"); print "\n\ttrack "$1"_cov\n\tparent ENCODE_UW_Histone\n\tshortLabel "t[5]" "t[6]"\n\tlongLabel ENCODE_UW_Histone "t[5]" "t[6]" "t[7]" "t[8]" "t[9]"\n\tvisibility dense\n\ttype bigWig\n\tbigDataUrl wgEncodeUwHistone/"$1".bw" }' >> trackDb.uw_histone.txt

echo -e "\ntrack ENCODE_UW_HistoneSummits\nshortLabel ENCODE_UW_HistoneSummits\nlongLabel ENCODE_UW_HistoneSummits\ngroup histone\ntype bigBed\ncompositeTrack on\nvisibility hide\n" >> trackDb.uw_histone.txt

ls ../hg38reps/wgEncodeUwHistone/*.bb | cut -f 4 -d"/" | cut -f 1 -d"." | sed 's/\([^[:blank:]]\)\([[:upper:]]\)/\1*\2/g'| awk '{split($1,t,"_"); print "\n\ttrack "$1"_sum\n\tparent ENCODE_UW_HistoneSummits\n\tshortLabel "t[5]" "t[6]"\n\tlongLabel ENCODE_UW_HistoneSummits "t[5]" "t[6]" "t[7]" "t[8]" "t[9]"\n\tvisibility dense\n\ttype bigBed\n\tbigDataUrl wgEncodeUwHistone/"$1".bb" }' >> trackDb.uw_histone.txt


sed 's/\*//g' trackDb.uw_histone.txt > temp.txt

cat temp.txt > trackDb.uw_histone.txt

mv trackDb.uw_histone.txt ../hg38reps/

####Meta summits
# Meta summits are peaks called on the genomic summit mappings to the Repeat Browser Consensus.
# I imposed a minimum score of 100 which seems to be rather stringent.  One can lower this ($5 >100) statement
# Meta summits calculated for all ChIP-SEQ collections.
#
###

ls ../hg38reps/imbeault2017/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do macs2 callpeak -t ../hg38reps/imbeault2017/${i}.bed --keep-dup all --outdir macs2 --nomodel --call-summits -n ${i}; done
cd macs2
cat *_summits* | awk '$5 > 100'| sed "s/_Trono_hg38reps_peak_.*\t/\t/g" | sort -k5 -n > imbeault2017_summit_summary.bed
bedSort imbeault2017_summit_summary.bed imbeault2017_summit_summary.bed 
bedtools slop -i imbeault2017_summit_summary.bed -g ../../hg38reps/hg38reps.sizes -b 10 > temp.bed
mv temp.bed imbeault2017_summit_summary.bed

ls ../hg38reps/imbeault2017/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do macs2 callpeak -t ../hg38reps/imbeault2017/${i}.bed --keep-dup all --outdir macs2 --nomodel --call-summits -n ${i}; done

cat *_summits* | awk '$5 > 100' | sort -k5 -n > imbeault2017_summit_summary.bed
bedSort imbeault2017_summit_summary.bed imbeault2017_summit_summary.bed 
bedtools slop -i imbeault2017_summit_summary.bed -g ../../../hg38reps/hg38reps.sizes -b 10 > temp.bed

awk '{if($5>1000){$5=1000}}{printf ("%s\t%.0f\t%.0f\t%s\t%.0f\t%s\n", $1, $2, $3, $4, $5, "+")}' imbeault2017_summit_summary.bed > temp.bed

awk -v OFS="\t" '{{split($4,a,"_")}{$4=a[1]; print $0}}' temp.bed > imbeault2017_summit_summary.bed

bedToBigBed imbeault2017_summit_summary.bed ../../../hg38reps/hg38reps.sizes imbeault2017_summit_summary.bb -type=bed6


###schmitges
ls ../hg38reps/schmittges2016/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do macs2 callpeak -t ../hg38reps/schmittges2016/${i}.bed --keep-dup all --outdir macs2/schmittges2016 --nomodel --call-summits --extsize 30 -n ${i}; done
cd macs2/schmittges2016
cat *_summits* | awk '$5 > 100' | sort -k5 -n > schmittges2016_summit_summary.bed
bedSort schmittges2016_summit_summary.bed schmittges2016_summit_summary.bed 
bedtools slop -i schmittges2016_summit_summary.bed -g ../../../hg38reps/hg38reps.sizes -b 10 > temp.bed
awk '{if($5>1000){$5=1000}}{printf ("%s\t%.0f\t%.0f\t%s\t%.0f\t%s\n", $1, $2, $3, $4, $5, "+")}' temp.bed > schmittges2016_summit_summary.bed
bedSort schmittges2016_summit_summary.bed temp.bed

awk -v OFS="\t" '{{split($4,a,"_")}{$4=a[2]"_"a[3]; print $0}}' temp.bed > schmittges2016_summit_summary.bed
bedToBigBed schmittges2016_summit_summary.bed ../../../hg38reps/hg38reps.sizes schmittges2016_summit_summary.bb -type=bed6


####tsankov
mkdir macs2/tsankov2014
ls ../hg38reps/tsankov2014/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do awk '{print $0 "\t+"}' ../hg38reps/tsankov2014/$i.bed > tmp.bed; cat tmp.bed > ../hg38reps/tsankov2014/$i.bed; done
ls ../hg38reps/tsankov2014/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do macs2 callpeak -f BED -t ../hg38reps/tsankov2014/${i}.bed --keep-dup all --outdir macs2/tsankov2014 --nomodel --call-summits --extsize 30 -n ${i}; done

cd macs2/tsankov2014
cat *_summits* | awk '$5 > 100' | sort -k5 -n > tsankov2014_summit_summary.bed
bedSort tsankov2014_summit_summary.bed tsankov2014_summit_summary.bed 
bedtools slop -i tsankov2014_summit_summary.bed -g ../../../hg38reps/hg38reps.sizes -b 10 > temp.bed
awk '{if($5>1000){$5=1000}}{printf ("%s\t%.0f\t%.0f\t%s\t%.0f\t%s\n", $1, $2, $3, $4, $5, "+")}' temp.bed > tsankov2014_summit_summary.bed
bedSort tsankov2014_summit_summary.bed temp.bed

awk -v OFS="\t" '{{split($4,a,"_")}{$4=a[2]"_"a[3]; print $0}}' temp.bed > tsankov2014_summit_summary.bed
bedToBigBed tsankov2014_summit_summary.bed ../../../hg38reps/hg38reps.sizes tsankov2014_summit_summary.bb -type=bed6


##encode histone
mkdir macs2/wgEncodeUwHistone

ls ../hg38reps/wgEncodeUwHistone/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do awk -v OFS="\t" '{$6="+"}{print $0}' ../hg38reps/wgEncodeUwHistone/$i.bed > tmp.bed; cat tmp.bed > ../hg38reps/wgEncodeUwHistone/$i.bed; done
ls ../hg38reps/wgEncodeUwHistone/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do macs2 callpeak -f BED -t ../hg38reps/wgEncodeUwHistone/${i}.bed --keep-dup all --outdir macs2/wgEncodeUwHistone --nomodel --call-summits --extsize 30 -n ${i}; done

cd macs2/wgEncodeUwHistone
cat *_summits* | awk '$5 > 100' | sort -k5 -n > wgEncodeUwHistone_summit_summary.bed
bedSort wgEncodeUwHistone_summit_summary.bed wgEncodeUwHistone_summit_summary.bed 
bedtools slop -i wgEncodeUwHistone_summit_summary.bed -g ../../../hg38reps/hg38reps.sizes -b 10 > temp.bed
awk '{if($5>1000){$5=1000}}{printf ("%s\t%.0f\t%.0f\t%s\t%.0f\t%s\n", $1, $2, $3, $4, $5, "+")}' temp.bed > wgEncodeUwHistone_summit_summary.bed
bedSort wgEncodeUwHistone_summit_summary.bed temp.bed


cut -f 4 temp.bed | sed 's/\([^[:blank:]]\)\([[:upper:]]\)/\1*\2/g' | tr "_" "*"| awk '{split($1,a,"\*")}{$1=a[5]"_"a[6]"_"a[9]"_"a[10]}{print $0}' | sed s/"Hotspots_"//g | sed s/"_hg38reps"//g > temp.txt
cut -f 1-3 wgEncodeUwHistone_summit_summary.bed  > temp1.txt 
cut -f 5-6 wgEncodeUwHistone_summit_summary.bed > temp3.txt
paste temp1.txt temp.txt temp3.txt > wgEncodeUwHistone_summit_summary.bed

bedToBigBed wgEncodeUwHistone_summit_summary.bed ../../../hg38reps/hg38reps.sizes wgEncodeUwHistone_summit_summary.bb -type=bed6

##encode TFBS
mkdir macs2/encRegTfbsClustered


ls ../hg38reps/encRegTfbsClustered/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do macs2 callpeak -f BED -t ../hg38reps/encRegTfbsClustered/${i}.bed --keep-dup all --outdir macs2/encRegTfbsClustered --nomodel --call-summits --extsize 30 -n ${i}; done

cd macs2/encRegTfbsClustered
cat *_summits* | awk '$5 > 100' | sort -k5 -n > encRegTfbsClustered_summit_summary.bed
bedSort encRegTfbsClustered_summit_summary.bed encRegTfbsClustered_summit_summary.bed 
bedtools slop -i encRegTfbsClustered_summit_summary.bed -g ../../../hg38reps/hg38reps.sizes -b 10 > temp.bed
awk '{if($5>1000){$5=1000}}{printf ("%s\t%.0f\t%.0f\t%s\t%.0f\t%s\n", $1, $2, $3, $4, $5, "+")}' temp.bed > encRegTfbsClustered_summit_summary.bed
bedSort encRegTfbsClustered_summit_summary.bed temp.bed

awk -v OFS="\t" '{{split($4,a,"_")}{$4=a[2]"_"a[3]; print $0}}' temp.bed > encRegTfbsClustered_summit_summary.bed

bedToBigBed encRegTfbsClustered_summit_summary.bed ../../../hg38reps/hg38reps.sizes encRegTfbsClustered_summit_summary.bb -type=bed6


#dnase
mkdir macs2/wgEncodeAwgDnaseMasterSites


ls ../hg38reps/wgEncodeAwgDnaseMasterSites/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do macs2 callpeak -f BED -t ../hg38reps/wgEncodeAwgDnaseMasterSites/${i}.bed --keep-dup all --outdir macs2/wgEncodeAwgDnaseMasterSites --nomodel --call-summits --extsize 30 -n ${i}; done

cd macs2/wgEncodeAwgDnaseMasterSites
cat *_summits* | awk '$5 > 100'| sort -k5 -n > wgEncodeAwgDnaseMasterSites_summit_summary.bed
bedSort wgEncodeAwgDnaseMasterSites_summit_summary.bed wgEncodeAwgDnaseMasterSites_summit_summary.bed 
bedtools slop -i wgEncodeAwgDnaseMasterSites_summit_summary.bed -g ../../../hg38reps/hg38reps.sizes -b 10 > temp.bed
awk '{if($5>1000){$5=1000}}{printf ("%s\t%.0f\t%.0f\t%s\t%.0f\t%s\n", $1, $2, $3, $4, $5, "+")}' temp.bed > wgEncodeAwgDnaseMasterSites_summit_summary.bed
bedSort wgEncodeAwgDnaseMasterSites_summit_summary.bed temp.bed

awk -v OFS="\t" '{{split($4,a,"_")}{$4=a[4]"_"a[5]; print $0}}' temp.bed > wgEncodeAwgDnaseMasterSites_summit_summary.bed

bedToBigBed wgEncodeAwgDnaseMasterSites_summit_summary.bed ../../../hg38reps/hg38reps.sizes wgEncodeAwgDnaseMasterSites_summit_summary.bb -type=bed6

##theunissen2016 
mkdir macs2/theunissen2016 


ls ../hg38reps/theunissen2016/*.bed | cut -f 4 -d "/" | cut -f 1 -d"." | while read i; do macs2 callpeak -f BED -t ../hg38reps/theunissen2016/${i}.bed --keep-dup all --outdir macs2/theunissen2016 --nomodel --call-summits --extsize 30 -n ${i}; done

cd macs2/theunissen2016
cat *_summits* | awk '$5 > 100' | sort -k5 -n > theunissen2016_summit_summary.bed
bedSort theunissen2016_summit_summary.bed theunissen2016_summit_summary.bed 
bedtools slop -i theunissen2016_summit_summary.bed -g ../../../hg38reps/hg38reps.sizes -b 10 > temp.bed
awk '{if($5>1000){$5=1000}}{printf ("%s\t%.0f\t%.0f\t%s\t%.0f\t%s\n", $1, $2, $3, $4, $5, "+")}' temp.bed > theunissen2016_summit_summary.bed
bedSort theunissen2016_summit_summary.bed temp.bed

awk -v OFS="\t" '{{split($4,a,"_")}{$4=a[1]"_"a[2]"_"a[3]; print $0}}' temp.bed > theunissen2016_summit_summary.bed

bedToBigBed theunissen2016_summit_summary.bed ../../../hg38reps/hg38reps.sizes theunissen2016_summit_summary.bb -type=bed6
#################