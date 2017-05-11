# This script processes the reference transcripts from the most recent Human genome assembly to pull out all GIMAP sequence transcripts, 
# find their genomic location (chromosome) using blast and identifying top hits, and uses biomaRt in R 

#data files loaded in home directory /Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/GENOME_TRANSCRIPTOME_SEQUENCES/human\ genome/GRCh38/Big_Data_Final_project

#Step 1: Put all sequences with a GIMAP header in one file.
	#use awk command to look at lines starting with >, and that have GIMAP somewhere in the line, and then print that until the next > symbol.

#!/bin/bash
awk '/^>/ { p = ($0 ~ />*GIMAP*/)} p'  interim_GRCh38.p10.RefSeqs.fasta > H_sapiens_GIMAP_IDs.fa

#Step 2: identify any GIMAP-like sequences that have yet to be identified in the genome by using the GIMAP 5 sequence as a BLAST (this is a well studied GIMAP member, see reference 1)
# query and the ref_seq_transcripts file as the database. 
awk '/^>/ { p = ($0 ~ />*GIMAP5*/)} p' H_sapiens_GIMAP_IDs.fa > H_sapiens_GIMAP5_seq.fa #find sequence for GIMAP5
/usr/local/ncbi/blast/bin/makeblastdb -in interim_GRCh38.p10.RefSeqs.fasta -title H_sapiens_ref_seq.db -dbtype nucl -out H_sapiens_ref_seq.db #make a blast db with the Ref_seqs file
/usr/local/ncbi/blast/bin/blastn -query H_sapiens_GIMAP5_seq.fa -db H_sapiens_ref_seq.db -out H_sapiens_GIMAP5_blastn_to_ref_seq.txt -task blastn-short #need to use blastn for nucleotide to nucleotide blast, then the blastn-short is needed to get good matches for short sequences

#Step 3: Add any sequences with an E value of < E-10, and a bitscore > 1000 to the file of collected GIMAP sequences
	#manual file inspection shows two entries with E value = 0, and all others with E value >10. None added to H_sapiens_GIMAP_IDs.fa file

#Step 4: Find the genomic location of each Seq ID in the H_sapiens_GIMAP_IDs.fa file
	#blast the transcript sequences to the reference genome sequences
	/usr/local/ncbi/blast/bin/makeblastdb -in GRCh38_latest_genomic.fna -title GRCh38_latest_genomic.db -dbtype nucl -out GRCh38_latest_genomic.db
	/usr/local/ncbi/blast/bin/blastn -query H_sapiens_GIMAP_IDs.fa -db GRCh38_latest_genomic.db -out H_sapiens_GIMAP_blastn_to_GRCh38_latest_genomic_eval0.txt -evalue 1e-150 -outfmt 7 #get as close to zero as possible and use the tabular format with comments in blastn
	
	#isolate each top hit
	cat H_sapiens_GIMAP_blastn_to_GRCh38_latest_genomic_eval0.txt |awk '/hits found/{getline;print}' | grep -v "#" > top_hits.txt
	echo Top hits found
	
	#make file with just the locations and Ids of each of the top hits (query acc, subject acc, start, end). Get multiple fields by separating them with commas
	cut -f 1,2,7,8 top_hits.txt > top_hit_acc_location.txt
	
	#Compile sequences for all unique subject acccessions
	 cut -f 2 top_hits.txt > top_hit_subject_acc.txt
	 awk '{for (i=1;i<=NF;i++) a[$i]++} END{for (i in a) printf i" ";print ""}' top_hit_subject_acc.txt >top_hit_subject_acc_unique.txt #take out repeat sequences in list to find unique
	 cat top_hit_subject_acc_unique.txt
	 awk '/^>/ { p = ($0 ~ />NC_000007.14/)} p' GRCh38_latest_genomic.fna  > top_hit_subject_acc_seq.txt #find first unique sequence
	 awk '/^>/ { p = ($0 ~ />NC_000004.12/)} p' GRCh38_latest_genomic.fna  >> top_hit_subject_acc_seq.txt #append second unique sequence to first 
	 
#Step 5: #Use biomaRt in R for that accession and find all the ensembl proteins flanking the GIMAP region, and compiles the collected data to be used for further downstream analysis
	#write R script to look at the annotated human genome for chromosome 4 and 7 
	Rscript session-info.R
	Rscript print-args.R
	#use position data from top_hit_acc_location.txt, print out all results to the command line to be saved into text file and further processed. 
	Rscript analysis.R 
	

#reference 1: Schulteis RD, Chu H, Dai X, Chen Y, Edwards B, Haribhai D, Williams CB, et al. Impaired survival of peripheral T cells, disrupted NK/NKT cell development, and liver failure in mice lacking Gimap5. Blood. 2008;112:4905–4914.
#Reference 2: For help with awk commands, http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/linux.html
#Reference 3: For help with step #5 http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/subsetFASTAFASTAQ.html
#reference 4: http://pythonforbiologists.com/index.php/manipulating-blast-output-with-command-line-tools/
#reference 5: http://www.bioconductor.org/packages//2.7/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
#reference 5: Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), "A greedy algorithm for aligning DNA sequences", J Comput Biol 2000; 7(1-2):203-14. 
#reference 6: http://genomicsclass.github.io/book/pages/bioc1_annoOverview.html
#reference 7: https://rachelss.github.io/BigDataAnalysis/12-cmdline.html