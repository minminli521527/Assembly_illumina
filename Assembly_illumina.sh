### 1. Download genome data, build index ###
# For clarity of expression and convenience of operation, use bgzip to decompress this sequence file and rename it to E.coli_K12_MG1655.fa.
gzip -dc GCF_000005845.2_ASM584v2_genomic.fna.gz > E.coli_K12_MG1655.fa
# Use samtools to create an index for it, which is to prepare for other data analysis tools (such as GATK) to quickly obtain any sequence on fasta.
/Tools/common/bin/samtools faidx E.coli_K12_MG1655.fa
#At this time, an E.coli_K12_MG1655.fa.fai file will be generated. In addition to convenient other tools, you can get the sequence of any position in the fasta file or any complete chromosome sequence through such an index. It is easy to complete the extraction of the specific region sequence of the reference sequence (or any fasta file). for example:
samtools faidx E.coli_K12_MG1655.fa NC_000913.3:1000000-1000200
# Obtained this sequence of E. coli K12 reference sequence: >NC_000913.3:1000000-1000200.



### 2. Download sequencing data ###
/Tools/sratoolkit/2.8.2/bin/fastq-dump --split-files SRR1770413
### SRA converted to fastq
/Tools/sratoolkit/2.8.2/bin/fastq-dump --split-files SRR1770413.sra
# Batch conversion
for id in SRR*
do
echo $id
/home/bio-soft/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump --split-3 $id
done
### Use bgzip (gzip is not recommended) to compress it into a .gz file, which can save space and will not affect the subsequent data analysis.
/Tools/common/bin/bgzip -f SRR1770413_1.fastq
/Tools/common/bin/bgzip -f SRR1770413_1.fastq
# Batch compression
ls *.fastq|while read id;do bgzip $id;done
### Linux environment SRA toolkit installation
# download
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
#Unzip
tar -xvzf sratoolkit.2.9.2-ubuntu64.tar.gz
# Check (you can see all kinds of software)
ls sratoolkit.2.9.2-ubuntu64/bin
# Configure environment variables (depending on individual circumstances, slightly different)
# (1)vi ~/.bashrc
# (2)export PATH=$PATH:/public/home/zffang/software/sratoolkit.2.9.2-ubuntu64/bin
# (3)source ~/.bashrc (make the configuration effective)



### 3. Data Filtering ###
###3.1 SGA filtering
#(1) Install SGA
conda install -c bioconda sga
#(2) Remove reads containing N in the PE library data. Here, only the PE data that needs to be assembled by contig is processed.
sga preprocess -o pe180.pp.fastq --pe-mode 1 pe180.1.fastq pe180.2.fastq
sga preprocess -o pe500.pp.fastq --pe-mode 1 pe500.1.fastq pe500.2.fastq
#(3) Use the index command to build an index on the reads data. When using multiple sets of data for contig assembly, use merge to merge the indexes of multiple sets of data. Example:
sga index -a ropebwt -t 8 --no-reverse pe180.pp.fastq
sga index -a ropebwt -t 8 --no-reverse pe500.pp.fastq
sga merge -t 8 -p pe -r pe180.pp.fastq pe500.pp.fastq
#(4) Use the correct command to correct reads
sga correct -t 8 -e 0.04 -m 45 -c 5 -o pe.ec.fastq pe.fastq
#(5) The filter of reads removes identical duplicates in the data; removes reads with low frequency kmer.
#First reconstruct the index of the corrected data
sga index -a ropebwt -t 8 pe.ec.fastq
sga filter -x 2 -t 8 -o pe.ec.f.fastq --homopolymer-check --low-complexity-check pe.ec.fastq



###3.2 Trim-galore filter
#(1) Install Trim-galore
conda create -n trim-galore python=2.7
conda activate trim-galore
conda install -c bioconda trim-galore -y
#(2) Use trim_galore to control the batch quality of double-ended *gz files and output to the current folder
cd /home/minminli/Bioinformatics/assembly/illumina/data
ls ./*_1.fastq.gz> ./1
ls ./*_2.fastq.gz> ./2
paste ./1 ./2> ./config
#New loop script qc.sh
cat> ./qc.sh
#Write the loop in the qc.sh script:
source trim-galore
dir='/home/minminli/Bioinformatics/assembly/illumina/data'
cat /home/minminli/Bioinformatics/assembly/illumina/data/config |while read id
do
           arr=($id)
           fq1=${arr[0]}
           fq2=${arr[1]}
trim_galore -q 15 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2
done
source deactivate trim-galore
#Run script
bash /home/minminli/Bioinformatics/assembly/illumina/data/qc.sh config

###3.3 Use fastqc+multiqc software to check data quality
find ./ -name "*.gz" |xargs fastqc -t 3 -o ./result
multiqc *fastqc.zip --pdf
#The results are saved in /home/minminli/Bioinformatics/assembly/illumina/data.



### 4. SPAdes genome assembly ###
#(1) Install SPAdes
conda install -c bioconda spades
# Run the test
spades.py --test
#(2)Assembly
spades.py -t 2 -k 21,33,55,77 \
--pe1-1 SRR3156160_pe_1_val_1.fq.gz --pe1-2 SRR3156160_pe_2_val_2.fq.gz \
--pe2-1 SRR3166543_pe_1_val_1.fq.gz --pe2-2 SRR3166543_pe_2_val_2.fq.gz \
--mp1-1 SRR3156163_jump_MP8KB_1_val_1.fq.gz --mp1-2 SRR3156163_jump_MP8KB_2_val_2.fq.gz \
--mp2-1 SRR3156596_MP20KB_1_val_1.fq.gz --mp2-2 SRR3156596_MP20KB_2_val_2.fq.gz \
--mp3-1 SRR3157034_MPUNK_1_val_1.fq.gz --mp3-2 SRR3157034_MPUNK_2_val_2.fq.gz \
-o ./spades_output



### 5. SSPACE assembly contig-scaffold ###
#(1) Install SSPACE
#Download software http://www.baseclear.com/lab-products/bioinformatics-tools/sspace-standard/
tar zxf SSPACE-STANDARD-3.0_linux-x86_64.tar.gz
#After decompressing the software package, run the perl program in the software folder to run SSPACE. The software main catalog contains some software instructions and examples, among which the README file describes in detail.
./SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl
#If you run the software v3.0 on Ubuntu 14.04, an error will occur. After the configure source code is revised, it can run successfully. The two amendments are as follows:
#################################################################################
[cc lang=”perl”] line124
#~ require “getopts.pl”;
use Getopt::Std;
[/cc]

[cc lang=”perl”] line136
#~ &Getopts(‘m:o:v:p:k:a:z:s:b:n:l:x:T:g:r:S:’);
getopt(‘m:o:v:p:k:a:z:s:b:n:l:x:T:g:r:S:’);
[/cc]
#################################################################################
#(2) Create libraries.txt, the content is as follows.
#################################################################################
lib1 bowtie SRR3156160_pe_1_val_1.fq.gz SRR3156160_pe_2_val_2.fq.gz 150 0.25 FR
lib1 bowtie SRR3166543_pe_1_val_1.fq.gz SRR3166543_pe_2_val_2.fq.gz 200 0.25 FR
lib2 SRR3157034_MPUNK_1_val_1.fq.gz SRR3157034_MPUNK_2_val_2.fq.gz 3000 0.5 RF
lib2 SRR3156163_jump_MP8KB_1_val_1.fq.gz SRR3156163_jump_MP8KB_2_val_2.fq.gz 8000 0.5 RF
lib2 SRR3156596_MP20KB_1_val_1.fq.gz RR3156596_MP20KB_2_val_2.fq.gz 20000 0.5 RF
unpaired bowtie unpaired_reads1.fasta
#################################################################################
#(3)Assemble contig-scaffold
perl /home/minminli/software/SSPACE/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl -l libraries.txt -s contigs_abyss.fa -k 3 -a 0.3 -x 0 -m 32 -o 20 -n 15 -p 1 -v 0 -b SSPACE_scaffolds 



### 6. Statistics contig or scaffolds information ###
#6.1 QUAST software
#(1) Install QUAST
conda create -n quast python=3.5
conda install -c bioconda quast
#(2)Run
quast.py scaffold.fasta -t 8 -o ./quast_outputfile

#6.2 BUSCO software
#(1) Install BUSCO
conda create -n busco python=3.7
conda install -c bioconda -c conda-forge busco
#Check whether the python file in the miniconda/envs/busco/bin directory is busco (BUSCO 4.0.5) or run_BUSCO.py (BUSCO 3.0.2)
#Run the following code to display the version (BUSCO 3.0.2)
run_BUSCO.py -h
#Or run the following code to display the version (BUSCO 4.0.5)
busco -h
#(2)Download the corresponding database file
# Browser open the following URL, https://busco.ezlab.org/, download the corresponding file according to BUSCO 4.0.5(synechococcales_odb10)/BUSCO 3.0.2(bateria_odb9) version, tar -xzvf decompress
#(3)Run
#Note: In the fasta file, the contig name generated by some assembly tools is in this form >contig/1/12345, and BUSCO will report an error when this fasta file is running
#BUSCO version 3.0.2 running
run_BUSCO.py -i scaffoldsK21.fasta -l bateria_odb9 -o busco -m genome
#Result interpretation: busco folder: Because busco is set in the -o option above, the folder name is suffixed. In this folder, one file is the most important. Is short_summary_busco.txt
#Or BUSCO version 4.0.5 running
busco -i scaffoldsK21.fasta -l synechococcales_odb10 -o busco -m genome
#Result interpretation: busco folder: Because busco is set in the -o option above, the folder name is suffixed. In this folder, one file is the most important. Is short_summary.specific.synechococcales_odb10.busco.txt
#Drawing: After the execution is completed, use generate_plot.py to draw, and gather all the short_summary_busco.txt results detected by BUSCO into a folder my_summaries
#BUSCO version 3.0.2 running
generate_plot.py -h
generate_plot.py –wd ./my_summaries
#Or BUSCO version 4.0.5 running
generate_plot.py -h
python3 /home/minminli/miniconda3/envs/busco/bin/generate_plot.py -wd ./my_summaries



### 7. Variation analysis ###



### 8. Genome annotation ###
#8.1 Repeat sequence annotation
#8.1.1 Based on the ab initio annotation method
#(1) RepeatModeler software installation
conda create -n annotation
conda activate annotation
conda install -c bioconda repeatmasker -y
RepeatMasker -h
conda install -c bioconda repeatmodeler -y
RepeatModeler -h
#(2)Run
#(2.1) Build a genome index database
#-name Arabidopsis_assembly: prefix, indicates the name of the database; -engine ncbi: indicates the use of rmblast
BuildDatabase -name Arabidopsis_assembly -engine ncbi Arabidopsis_assembly.fasta
#(2.2) Build library
#Result: This step is very time-consuming. The results of the operation are stored in the RM_.xxx folder, and the two main files Arabidopsis_assembly-families.fa/consensi.fa and Arabidopsis_assemblyh-families.stk are obtained. The former is the consensus repeat sequence found, and the sequence id will indicate which repeat it belongs to. Sequence family, if it cannot be classified, use the label "Unkown"; the latter is a seed alignment file in Stockholm format, which can be uploaded to the Dfam_consensus database using util/dfamConsensusTool.pl.
#Parameter: -database: the same as the previous step; -engine: the same as the previous step; -pa: the number of threads.
RepeatModeler -pa 2 -database Arabidopsis_assembly -engine ncbi
#(2.3) Combine RepeatMasker to find repetitive sequences in the genome
#Use the library consensus sequence obtained in the previous step as the reference repeat sequence, and use RepeatMasker to find the repeat sequence in the target genome.
#Result: The output path specified by -dir must be established in advance, otherwise there will be no result; -poly to identify and annotate microsatellites; -nolow to remove low-complexity repetitive sequences; -norna to remove small RNA sequences; -gff to organize into standard gff files structure
mkdir repeat_denovo
RepeatMasker -nolow -no_is -norna -engine ncbi -pa 2 -html -gff -poly -lib Arabidopsis_assembly-families.fa -dir ./repeat_denovo Arabidopsis_assembly.fasta
#result:
#1. *.out.gff: Repetitive sequence genome annotation file, similar to gene annotation, the most important result. Through this *.out/*.out.gff file, you can locate the position of the special type of repetitive sequence element that you want to pay attention to in the genome, and then you can write a script to extract this sequence based on the position information, or Further study their functions and so on.
#2. *.tbl: Repetitive sequence annotation result report, "bases masked" is the total length of the repetitive sequence and the proportion in the genome
#3. *.out.html: Detailed results of web version, same as RepeatMasker online annotation result report
#4. *.masked: Replace the large items annotated as repetitive sequence regions with the genome of N
#5. *.out: RepeatMasker default input result format, the information is basically related to gff. SW score The score based on the Smith-Waterman algorithm comparison; the substitution rate of the Div% comparison interval compared with the consensus sequence; the percentage of base deletions (deleted bases) in the query sequence for Del%; Ins% in the repeat library sequence Percentage of missing bases in the middle (bases inserted); the repetitive sequence to be masked in the Query sequence input; Position begin; Position end; Query left The number of bases in the query sequence that exceeds the upper region, += compared to the library The sense strand of the repeated sequence, if it is complementary, use "c" to indicate;
Matching repeat The name of the repeating sequence; Repeat family(class) The type of the repeating sequence; Position begin;
Position end; The number of bases from the left end of the Query left alignment region to the left end of the repeated sequence, the sequence ID of the alignment
#6. *.cat.gz: Sequence and repeat sequence alignment file
#7. *.poly.out: The microsatellite annotations in the prediction result *.out are identified and organized into a table separately, with the same file structure as *.out
#(2.4) Optional, the previous step has been generated through the -gff parameter
# Organize the content in *.out into the structure type of the standard gff file (-gff parameter generation). It mainly contains information such as the position and structure of the repeated sequence.
less -S *.out
perl repeat_to_gff.pl *.out
less -S *.out.gff



#8.1.2 Based on the same origin annotation method
#(1) RepeatMasker software installation
conda activate annotation
conda install -c bioconda repeatmasker -y
RepeatMasker -h
#(2)Run
#(2.1) Create a result output directory
mkdir repeat_homology
#(2.2) Run the program
#Parameter: parallel is the number of threads to choose; species is the name of the species, see help for common species, write lowercase Latin genus names or full names in quotation marks if there are none; html and gff are the output results in html and gff format for easy viewing and downstream analysis; dir output Result directory; genome fa file must be placed at the end of all parameters; Species specified by -species, otherwise the human repetitive sequence database is the default comparison. For conventional model species or taxa, you can directly enter the general name. A safer way is Add the Latin name of the species after the "-species" parameter. Because there is a space between the genus name and the species name, it needs to be enclosed in quotation marks. If RepeatMasker does not prompt an exception, that is, there is a target species in the library, it will run normally, wait The result is fine; -engine can also not specify ncbi.
RepeatMasker -nolow -no_is -norna -engine ncbi -pa 2 -species arabidopsis -html -gff -dir ./repeat_homology Arabidopsis_assembly.fasta
#(2.3)Result
#1. *.out.gff: Repetitive sequence genome annotation file, similar to gene annotation, the most important result
#2. *.tbl: Repetitive sequence annotation result report information summary table overview
#3. *.out.html: Detailed results of web version, same as RepeatMasker online annotation result report
#4. *.masked: Replace the large items annotated as repetitive sequence regions with the genome of N
#5. *.out: RepeatMasker default input result format, the information is basically related to gff
#6. *.cat.gz: Sequence and repeat sequence alignment file

#8.1.3 Search for tandem repeats in DNA sequences
#(1) TRF software installation
conda activate annotation
conda install -c bioconda trf -y
trf -h
#(2)Run
#(2.1) Run
mkdir repeat_trg
trf Arabidopsis_assembly.fasta 2 7 7 80 10 50 500 -f -d -m
#(2.2)View results
#1-2 columns; 3 columns of repeat unit length; 4 columns of repeat unit copy number; 5 columns of repeat unit length; 6 columns of matching percentage; 7 columns of indel percentage; 8 columns of scores; 9-12 columns of ACGT base numbers ; 13 columns of moisture value; 14 columns of unit repeat sequence; 15 columns of entire repeat sequence
less -S Arabidopsis_assembly.fasta.2.7.7.80.10.50.500.dat
#(2.3)dat-gff conversion
perl repeat_to_gff.pl Arabidopsis_assembly.fasta.2.7.7.80.10.50.500.dat
less -S Arabidopsis_assembly.fasta.2.7.7.80.10.50.500.dat.gff


# 8.2 Gene structure prediction
##8.2.1 denovo based on ab initio prediction
#First construct the repeat-mask genome, on this basis, use August, Genescan, GlimmerHMM, Geneid and SNAP to predict the coding region
#Predict the gene structure through the existing probability model, and the accuracy in predicting the splicing site and UTR region is low

## 8.2.1.1 AUGUSTUS
#Generate gff file and protein sequence
#(1) AUGUSTUS software installation
conda activate busco
augustus
#(2) Run
#View all species information
augustus --species=help
#Run the following series of commands to get the standard gff format file Arabidopsis_assembly.fasta.masked.temp.joined.augustsus.gff
augustus --species=arabidopsis Arabidopsis_assembly.fasta.masked> Arabidopsis_assembly.fasta.masked.gff
join_aug_pred.pl <Arabidopsis_assembly.fasta.masked.gff | grep -v'^#'> Arabidopsis_assembly.fasta.masked.temp.joined.gff
bedtools sort -i Arabidopsis_assembly.fasta.masked.temp.joined.gff> Arabidopsis_assembly.fasta.masked.temp.joined.augustsus.gff
less Arabidopsis_assembly.fasta.masked.temp.joined.augustsus.gff
###!!! Note the conversion format: the gff format generated in the previous step is the gtf format:
/home/minminli/miniconda3/envs/maker/bin/genemark_gtf2gff3 cut_Arabidopsis.fasta.masked.gff> cut_Arabidopsis.fasta.finaly.gff
#or:
/home/minminli/miniconda3/envs/maker/bin/genemark_gtf2gff3 Arabidopsis_assembly.fasta.masked.temp.joined.augustsus.gff> cut_Arabidopsis.fasta.masked.finaly.gff
# Extract the amino acid fasta sequence from the initial Arabidopsis_assembly.fasta.masked.gff file:
sed -n'/^#/p' Arabidopsis_assembly.fasta.masked.gff | sed -n'/start/,/\]/p' | sed's/# start gene />/g;s/protein sequence \ = \[//g;s/#//g;s/\]//g;s/^\s//g'> Arabidopsis_assembly.fasta.masked.pep.fa
#(3) Or without using RNA-seq data, use the following method for multiple chromosomes in parallel Augustus
mkdir 01-augustsus && cd 01-augustsus
ln ../Arabidopsis_assembly.fasta.masked
#Result file is in genome.fa.split
seqkit split -i Arabidopsis_assembly.fasta.masked
#Parallel processing
find Arabidopsis_assembly.fasta.masked.split/ -type f -name "*.masked" | parallel -j 30 augustus --species=arabidopsis --gff3=on >> temp.gff
join_aug_pred.pl < temp.gff  | grep -v '^#' > temp.joined.gff
bedtools sort -i temp.joined.gff > augustsus.gff
less augustsus.gff

## 8.2.1.2 GeneMark-ES
#Generate gff file
#If after the program is executed, no *.gtf file is generated or the file is empty, it means that an error occurred during the program execution. The GeneMark-ES program uses the motif search method to derive the branch point model. GeneMark-ES needs a lot of memory and hard disk space in the process of execution.
#(1) GeneMark-ES software installation
#Software link: http://exon.gatech.edu/GeneMark/, register to download, decompression is available
tar xzvf gmes_linux_64.tar.gz
cd gmes_linux_64
#gm_key is the pass of the software and needs to be copied to the home directory
gunzip gm_key_64.gz
cp gm_key_64 ~/.gm_key_64
#(2) An error was reported when running ./gmes_petap.pl of GeneMark-ET, install the following perl modules
#Install cpan
yum -y install perl
yum -y install perl-CPAN
#Run cpan, then add force before installation
cpan
force install YAML
force install Hash::Merge
force install Logger::Simple
force install Parallel::ForkManager
#If force cannot be installed, run the following command in the terminal to install:
curl -L http://cpanmin.us | perl---sudo YAML
curl -L http://cpanmin.us | perl---sudo Hash::Merge
curl -L http://cpanmin.us | perl---sudo Logger::Simple
curl -L http://cpanmin.us | perl---sudo Parallel::ForkManager
#(3) Run the software
#Use the gmes_petap.pl script, --ES is the algorithm, --cores 2 is the thread, and *.masked is the genome with repetitive sequences removed fa
/home/minminli/Bioinformatics/software/genemark-es/gmes_linux_64/gmes_petap.pl --ES --cores 2 --sequence cut_Arabidopsis.fasta.masked
#Generate genemark_es.gtf file
less genemark_es.gtf
#Use genemark_gtf2gff3 to convert gtf to gff
/home/minminli/miniconda3/envs/maker/bin/genemark_gtf2gff3 genemark_es.gtf> genemark_es.gff
#(4) In addition, you can use the --ET algorithm to run transcriptome data
#Where seq.fasta.masked is the genome sequence of your species, introns.gff is the location information of recording introns, and --cores is the multithreading parameter.
gmes_petap.pl --ET introns.gff --sequence seq.fasta.masked --cores 60
#How to get introns.gff?
#Method 1: Bet_to_gff.pl is provided under the TopHat2 software catalog of the RNA-seq comparison tool to convert the junctions.bed in the output result of TopHat2 to introns.gff
bed_to_gff.pl --bed tophat_out/junctions.bed --gff introns.gff --label tophat2
#Method 2: The star_to_gff.pl script under the STAR software converts the STAR output file "SJ.out.tab" that records intron information to gff format
star_to_gff.pl --star SJ.out.tab --gff SJ.gff --label STAR

## 8.2.1.3 FragGenescan
#CDSRegion Forecast
#(1) FragGenescan software installation (Hidden Markov Model model)
conda activate gene-structure-prediction
conda install -c bioconda fraggenescan
FragGeneScan
#(2)Run
run_FragGeneScan.pl -genome=./Arabidopsis_assembly.fasta.masked -out=./Arabidopsis_assembly.fasta.masked  -complete=0  -train=illumina_5


##8.2.2 Homology-based prediction based on homologous sequence
#Progeins of relative species(Ensembl,NCBI,BLAST,Genewise)/Known proteins(SwissPort database)/Expression mapping(EST,cDNA,RNA-Seq:BLAT software)
#Download the complete protein set of several other representative animals and plants, use TblastN to align the protein sequence to the sequence of the preliminary assembly result, the threshold of E-value is 1e-5. Use Solar software for BLAST hits of different proteins To merge. GeneWise predicts the complete gene structure based on the corresponding gene region of each BLAST hit
#Auxiliary annotation through RNA-seq data of species can more accurately determine exon regions and splicing sites.
#Some gene proteins are highly conserved among similar species, so you can use the existing high-quality related species annotation information to determine exon boundaries and splicing sites by sequence alignment

## 8.2.2.1 Genewise(111 ???)
#(1) Genewise software installation
conda activate gene-structure-prediction
conda install -c bioconda wise2
#If an error is reported, add the environment variable of wisecfg
echo'export WISECONFIGDIR=/opt/biosoft/wise2.4.1/wisecfg/' >> ~/.bashrc
source ~/.bashrc
#(2)Run
genewise -help
genewise Arabidopsis_thaliana.TAIR10.pep.all.fa Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -pseudo -genes -both -gff> Arabidopsis_assembly.genewise.out.gff


##8.2.3 Prediction based on transcriptome (111 ???)
#Transcriptome prediction: Use Tophat to align the RNA-seq data to the sequence of the preliminary assembly result, and then use cufflinks to assemble the transcript to form a gene model.


##8.2.4 Evidencemodeler(EVM)(111 ???)
#Each method has its own advantages and disadvantages, so in the end, it needs to be integrated with EvidenceModeler (EVM) and GLEAN tools to merge into a complete gene structure. Based on reliable gene structure, functional annotations, protein functional domain annotations, gene ontology annotations, pathway annotations, etc. can be followed.



# 8.3 Gene function annotation
# 8.3.1 InterProScan
#(1) InterProScan 5.0 software installation
# (1.1) First download InterProScan 5.0, download address: ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/, select the latest version interproscan-5.41-78.0-64-bit.tar. gz and interproscan-5.41-78.0-64-bit.tar.gz.md5 can be, then download
#Use md5 value to verify whether the download is complete
md5sum -c interproscan-5.41-78.0-64-bit.tar.gz.md5
#Use md5 value to verify whether the download is complete
md5sum -c panther-data-14.1.tar.gz.md5

# (1.2) Then download Panther Models, download address ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/, also choose the latest version panther-data-14.1.tar.gz And panther-data-14.1.tar.gz.md5, then download
# (1.3) Before installing the above, first check whether the server meets the installation requirements:
        64-bit Linux
        Perl (default on most Linux distributions)
        Python 2.7.x only
        Oracle’s Java JDK/JRE version 8 (required by InterProScan 5.17-56.0 onwards). Earlier InterProScan release versions required Java 6 (version 6u4 and above) or Java 7.
        Environment variables set:
            $JAVA_HOME should point to the location of the JVM
            $JAVA_HOME/bin should be added to the $PATH
        For details, please check https://github.com/ebi-pf-team/interproscan/wiki/InstallationRequirements
# (1.4) After the configuration is complete, you can install it, in fact, just unzip it:
tar -zxvf interproscan-5.41-78.0-64-bit.tar.gz
cd interproscan-5.41-78.0-64-bit/interproscan-5.41-78.0/data/
#Remember to move the panther-data-14.1.tar.gz compressed package to the above directory first
tar -zxvf panther-data-14.1.tar.gz
# (1.5) Next is an optional option, depending on whether you need Match Lookup Service or not, because I am localized and do not want to operate on the Internet, so this operation will be prohibited
vim interproscan-5.41-78.0-64-bit/interproscan-5.41-78.0/interproscan.properties
# Then comment# Remove the following line of code
precalculated.match.lookup.service.url=http://www.ebi.ac.uk/interpro/match-lookup
# (1.6) After the above steps, the installation of InterProScan is basically completed, enter ./interproscan.sh in the terminal, and you will see the usage information
./interproscan.sh
#(2) Use of InterProScan 5.0
#parameter settings:
# -appl,–applications are used to specify which databases in Interpro are used, all databases are by default
# -b,–output-file-base is used to specify the path or folder of the output file, the default is the path of the input file
# -f,–formats is used to specify the suffix of the output file, the protein sequence will output TSV, XML and GFF3
# -i,–input The input file should generally be in fasta format without other special symbols. Both nucleic acid and protein can be used, but protein is recommended. After all, the protein file is relatively small.
#database：
# CDD
# COILS
# Gene3D
# HAMAP
# MOBIDB
# PANTHER
# Pfam
# PIRSF
# PRINTS
# ProDom
# PROSITE (Profiles and Patterns)
# SFLD
# SMART (unlicensed components only by default - this analysis has simplified post-processing that includes an E-value filter, however you should not expect it to give the same match output as the fully licensed version of SMART)
# SUPERFAMILY
# TIGRFAMs
# The following databases are available in interproscan 5, but require permission:
# Phobius (licensed software)
# SignalP
# SMART (licensed components)
# TMHMM
# After the installation, we can take the test file in the Interproscan folder to test, if there is no error, it means that InterProScan can run normally
./interproscan.sh -i test_proteins.fasta -f tsv
#Finally, just check the result. If there is no special requirement, the use of InterProScan 5.0 is just like that.
less test_proteins.fasta.tsv
# Analyze the protein result seed_proteins.faa produced by Genemak-ES:
cp seed_proteins.faa seed_proteins.fa
#Remove the * symbol in the sequence
sed -i "s/*//g" seed_proteins.fa
# ./interproscan.sh --goterms --pathways -i seed_proteins.fa
'/home/minminli/Bioinformatics/software/interproscan/interproscan-5.41-78.0-64-bit/interproscan-5.41-78.0/interproscan.sh' --goterms --pathways -i seed_proteins.fa
#View Results
less seed_proteins.fa.gff3
less seed_proteins.fa.tsv
less seed_proteins.fa.xml



# 8.3.2 KEGG
#KASS :http://www.genome.jp/tools/kaas/
#It is recommended to upload the amino acid sequence of the encoded protein for functional annotation. It is not recommended to use the nucleic acid sequence of the gene. In eukaryotes, due to variable splicing, one gene corresponds to multiple encoded proteins. Therefore, it is best to associate the corresponding gene id with the KO annotation result according to the corresponding relationship between the gene id and the protein id (which can be obtained from the genome annotation gff file).
#html”>“[BRITE hierarchies]”>“KEGG Orthology (KO)”, download the comment form.



# 8.4 ncRNA prediction
# 8.4.1 tRNAscan-SE predict tRNAs
#tRNAscan-SE is a very good method to detect tRNA pseudogenes.
# (1) tRNAscan-SE software installation
conda activate gene-structure-prediction
conda install trnascan-se
tRNAscan -h
# (2) Run
# View parameters
tRNAscan-SE
#Default uses eukaryotes; -O is suitable for mitochondria and chloroplasts. If this parameter is selected, only Cove will be used for analysis, the search speed will be very slow, and pseudogenes detection will not be given; -C only use Cove for tRNA analysis, although The accuracy is improved to a certain extent, but it will be extremely slow; -o saves the results to a file; -f saves the tRNA secondary structure results to a file; -m saves the statistical results to a file; cut_Arabidopsis.fasta.masked means remove Genome sequence after repetition.
tRNAscan-SE -o tRNA.cover.out -f rRNA.cover.ss -m tRNA.cover.stats cut_Arabidopsis.fasta.masked
#or
tRNAscan-SE -C -o tRNA.cover.out -f rRNA.cover.ss -m tRNA.cover.stats cut_Arabidopsis.fasta.masked
# (3) View results
# The main result is in *.out
less tRNA.cover.out
less rRNA.cover.ss
less tRNA.cover.stats



#8.4.2 barrnap predicts rRNA
# Compared with RNAmmer, in addition to basic bacteria, archaea, and eukaryotes, new mitochondrial rRNA predictions are added.
# (1) barrnap software installation
conda activate gene-structure-prediction
conda install barrnap
# (2) Run
# View parameters
barrnap
barrnap -h
#run
#--kingdom parameter specifies the species type, bac represents bacteria, arc represents archaea, euk represents eukaryotes, mito represents metazoan mitochondria; --threads specifies the number of parallel threads
barrnap --kingdom euk --threads 8 --quiet small.fna> rRNA.gff3
#The prediction result is saved in GFF3 format
less rRNA.gff3



# 8.4.3 infernal predict ncRNA
# (1) infernal software installation
conda activate gene-structure-prediction
conda install -c bioconda infernal
cmscan -h
cmpress -h
# (2) Download Rfam_ CM data file
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.1/Rfam.cm.gz
gunzip Rfam.cm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.1/Rfam.clanin
# (3) Use Infernal's program cmpress to compress Rfam.cm and create an index # It takes about one minute
cmpress Rfam.cm
# Output as
#Working...    done.
#Pressed and indexed 3016 CMs and p7 HMM filters (3016 names and 3016 accessions).
#Covariance models and p7 filters pressed into binary file:  Rfam14.1.cm.i1m
#SSI index for binary covariance model file:                 Rfam14.1.cm.i1i
#Optimized p7 filter profiles (MSV part)  pressed into:      Rfam14.1.cm.i1f
#Optimized p7 filter profiles (remainder) pressed into:      Rfam14.1.cm.i1p
# (4) Run
# Determine the nucleotide sequence to be queried or the size of the genome cut_Arabidopsis.fasta.masked, as a parameter of subsequent commands
# There is a line in the output result, Total # of residues: 1630948 is what we need.
esl-seqstat cut_Arabidopsis.fasta.masked
# Consider that the genome is double-stranded and the unit of the parameters used in the next step is Million. The calculation formula is: Z = total * 2 * CMmumber/1000000, so the number of models in the CM database must be calculated, which is used in Rfam14.0 version
#Calculate the number of models in the CM database and get a value of 6032. Z = total * 2 * CMmumber/1000000 = 1630948*2*6032/1000000 = 19675
less Rfam.cm | grep'NAME'|sort|wc -l
#run
# Rfam.clanin is the downloaded claninfo file, and the path must be provided
# Rfam.cm cm file downloaded
# cut_Arabidopsis.fasta.masked Sequence to be queried
# cut_Arabidopsis.fasta.masked.cmscan output result
# cut_Arabidopsis.fasta.masked.tblout Output results in tabular format
# -Z: The size of the query sequence, in M. Calculate by esl-seqstat or write a program yourself, remember to multiply by 2, multiply by the number of CM models, and divide by 1000000
#--cut_ga: Output results that are not less than the Rfam GA threshold. This is the threshold for Rfam to authenticate the RNA family, and sequence scores not lower than this threshold are considered true homologous sequences.
#--tblout: table output.
#--fmt 2: A format of table output.
cmscan -Z 19675 --cut_ga --rfam --nohmmonly --tblout cut_Arabidopsis.fasta.masked.tblout --fmt 2 --cpu 24 --clanin Rfam.clanin Rfam.cm cut_Arabidopsis.fasta.masked> cut_Arabidopsis.fasta. masked.cmscan
# For an input sequence of 500M size, 48 threads, it takes 7 hours, preferably in the background
nohup cmscan -Z 19675 --cut_ga --rfam --nohmmonly --tblout cut_Arabidopsis.fasta.masked.tblout --fmt 2 --cpu 24 --clanin Rfam14.1.clanin Rfam.cm cut_Arabidopsis.fasta.masked > cut_Arabidopsis.fasta.masked.cmscan &
# (5) Output result
In # filename/cut_Arabidopsis.fasta.masked.tblout, each row has 27 columns. The key ones are `target name`, `accession`, `query name`, `seq from`, `seq to`, `strand`, `E-value`, `score`.
The #`olp` column indicates the overlapping information of the query sequence, and `*` indicates that on the same chain, there is no sequence that overlaps the query sequence and there is a match in the Rfam database. This is the query sequence that needs to be retained. `^` means that there is no sequence on the same chain that matches the Rfam database better than this query sequence, and it needs to be retained. `=` means that there is a sequence that matches the Rfam database better than this query sequence on the same chain and should be ignored.
# The final result can be obtained by the following command.
grep -v'=' cut_Arabidopsis.fasta.masked.tblout> cut_Arabidopsis.fasta.deoverlapped.tblout
# (6) Format conversion
# Process the file into excel form, only keeping the information we need
#Get tblout.final, convert the result format, and extract the necessary columns and non-overlapping areas or areas with high scores in the overlapping areas.
#= indicates that there is a sequence that matches the Rfam database better than this query sequence on the same chain and should be ignored.
# Remove the comment line at the beginning of # Extract the required columns. Initially set the separator to a tab character.
#If it is the first line, print the name of the line starting from line 3, if the 20 column is not an equal sign, and does not start with a hash sign, print 2, 3, 4 and other columns to output
awk'BEGIN{OFS="\t";}{if(FNR==1) print "target_name\taccession\tquery_name\tquery_start\tquery_end\tstrand\tscore\tEvalue"; if(FNR>2 && $20!="= "&& $0!~/^#/) print $2,$3,$4,$10,$11,$12,$17,$18; }'cut_Arabidopsis.fasta.tblout> cut_Arabidopsis.fasta.tblout.final.xls
# (7) Download Rfam notes
#First download the notes of the Rfam family, click on the Rfam: Home page website, select SEARCH--Entry type--select all check boxes, Submit, and copy the obtained form and organize it into the format rfamannotation.txt divided by the TAB key. Extract the second semicolon information in the third column, this is each ncRNA class
cat rfamannotation.txt | awk'BEGIN {FS=OFS="\t"}{split($3,x,";");class=x[2];print $1,$2,$3,class}'> rfam_anno_class. TXT
# (8) Count the number of ncRNA
#Read file 1 first, and store column 2 and column 4, namely name and class, in dictionary a according to the relationship. Read file 2 again, each row, use its own column 1, which is the name, to query the dictionary, get the class, store it in the type variable, and count
awk'BEGIN{OFS=FS="\t"}ARGIND==1{a[$2]=$4;}ARGIND==2{type=a[$1]; if(type=="") type="Others "; count[type]+=1;}END{for(type in count) print type, count[type];}' rfam_anno_class.txt cut_Arabidopsis.fasta.deoverlapped.tblout.final.xls
# (9) Count the total length of each type of RNA
awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$2]=$4;}ARGIND==2{type=a[$1]; if(type=="") type="Others"; if($6=="-")len[type]+=$4-$5+1;if($6=="+")len[type]+=$5-$4+1}END{for(type in len) print type, len[type];}' rfam_anno_class.txt cut_Arabidopsis.fasta.deoverlapped.tblout.final.xls



### 9. Collinear drawing ###
#9.1 Sibelia
#(1) Install Sibelia
conda activate denovo
conda install -c bioconda sibelia
#(2)Run
Sibelia -h
Sibelia -s loose.fasta
circos -conf circos.conf

#9.2

MCScanX


