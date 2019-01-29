# Genomic Origin Through Taxonomic CHAllenge (GOTTCHA)

GOTTCHA is an application of a novel, gene-independent and signature-based metagenomic
taxonomic profiling method with significantly smaller false discovery rates (FDR) that is 
laptop deployable. Our algorithm was tested and validated on twenty synthetic and mock 
datasets ranging in community composition and complexity, was applied successfully to data
generated from spiked environmental and clinical samples, and robustly demonstrates 
superior performance compared with other available tools.

-------------------------------------------------------------------
## SYSTEM REQUIREMENT

Linux (2.6 kernel or later) or Mac (OSX 10.6 Snow Leopard or later) operating system
with minimal 8 GB of RAM is recommended. Perl v5.8 or above is required. The C/C++
compiling enviroment might be required for installing dependencies. Systems may vary.
Please assure that your system has the essential software building packages (e.g. build-essential
for Ubuntu, XCODE for Mac...etc) installed properly before running the installing 
script.

GOTTCHA was tested successfully on our Linux servers (Ubuntu 12.10 w/ Perl v5.14.2; 
Ubuntu 10.04 w/ Perl v5.10.1) and Macbook Pro laptops (MAC OSX 10.8 w/ XCODE v5.1).

-------------------------------------------------------------------
## QUICK START 

This is an example of profiling a "test.fastq" file using GOTTCHA with a species-level
pre-computed bacterial database. The testing FASTQ file comes along with the GOTTCHA package
in the "test" directory. More details are stated in the INSTRUCTION section.

1. Obtaining GOTTCHA package:

        $ git clone https://github.com/LANL-Bioinformatics/GOTTCHA.git gottcha

2. Installing GOTTCHA:

        $ cd gottcha
        $ ./INSTALL.sh

3. Downloading lookup table and species-level database from our web server:
 
        $ wget https://edge-dl.lanl.gov/gottcha/GOTTCHA_database_v20150825/GOTTCHA_lookup.tar.gz
        $ wget https://edge-dl.lanl.gov/gottcha/GOTTCHA_database_v20150825/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species.tar.gz
   
   If you have any difficulty obtaining the databases, please contact Po-E Li <po-e@lanl.gov>.

4. Unpacking and decompressing archives:

        $ tar -zxvf GOTTCHA_lookup.tar.gz
        $ tar -zxvf GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species.tar.gz
	
5. Running gottcha.pl: 

        $ bin/gottcha.pl             \
             --threads 8             \
             --outdir ./             \
             --input test/test.fastq \
             --database database/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species

6. Enjoying the result at './test.gottcha.tsv'.

-------------------------------------------------------------------
## DETAIL INSTRUCTIONS

The detail of steps in the above section will be descrbed in this section. Note that all 
instructions in this document use pre-computed databases downloaded from our web site.

If you are looking for instructions to build a CUSTOM database and/or running GOTTCHA
step-by-step, please read README_FULL.md.

-------------------------------------------------------------------
### Obtaining GOTTCHA

The source codes can be downloading from [here](https://bitbucket.org/poeli/gottcha).
The pre-computed databases need to be downloaded separately from our web server.
Please see below in the [Obtaining Pre-computed Databases] section.
       
You can use "git" to obtain the package:

        $ git clone https://github.com/poeli/gottcha

or download the compressed archive in
 [zip](https://github.com/poeli/gottcha/get/master.zip),
 [gz](https://github.com/poeli/gottcha/get/master.tar.gz) or 
 [bz2](https://github.com/poeli/gottcha/get/master.tar.bz2).

-------------------------------------------------------------------
### Installation

The GOTTCHA profiling and database-generating scripts are primarily Perl-based, and require at
least Perl v5.8 with dependencies installed properly (listed in README_FULL.md).
The splitrim tool is written in [D](http://www.dlang.org) that requires an appropriate D 
compiler to complie it. GOTTCHA utilizes [BWA](https://github.com/lh3/bwa) with the BWA-MEM algorithm
for read mapping. You can either keep "dmd" and "bwa" in your system path or simply 
run the installation script - INSTALL.sh. This script will check and try to install missing 
tools and dependencies:

    	$ ./INSTALL.sh

After running INSTALL.sh successfully, the binaries and related scripts will be stored
in the ./bin directory.

-------------------------------------------------------------------
### Obtaining Pre-computed Databases

Databases of unique genome segments at multiple taxonomic levels (e.g. family, species, 
genus, strain-level, etc.) are used for taxonomic classification of reads. Variants of 
these databases, in which all human 24-mers were removed were also generated and used in 
this study. These 24-mers were derived from the GRCh37.p10 (Genome Reference Consortium),
HuRef (J. Craig Venter Institute), and CHM1_1.0 (Washington U. School of Medicine) 
assemblies and include unplaced scaffolds. For example, 
GOTTCHA_BACTERIA_c3514_k24_u24_xHUMAN3x.species.tar.gz is a GOTTCHA bacterial 
species-level signature database that was produced by eliminating shared 24-mer (k24) sequences from
3514 bacterial replicons (c3514; includes both chromosomes and plasmids) and 3 human genomes (xHuman3X), 
while retaining a minimum of 24bp of unique fragments (u24).

The compressed database archives are available for users to download from our web server:
 
 > https://edge-dl.lanl.gov/gottcha/GOTTCHA_database_v20150825/

GOTTCHA requires a taxanomic lookup table (GOTTCHA_lookup.tar.gz) and a pre-computed
database (e.g: GOTTCHA_BACTERIA_c3514_k24_u24_xHUMAN3x.species.tar.gz) to classify reads.
These signature databases could be huge. We highly recommend that users also download the
corresponding *.md5 file for verification.

You can use the 'wget' command to download both archives, one at a time:

        $ wget https://edge-dl.lanl.gov/gottcha/GOTTCHA_database_v20150825/GOTTCHA_lookup.tar.gz
        $ wget https://edge-dl.lanl.gov/gottcha/GOTTCHA_database_v20150825/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species.tar.gz

Then use 'tar' to unpack and decompress both archives:

        $ tar -zxvf GOTTCHA_lookup.tar.gz
        $ tar -zxvf GOTTCHA_BACTERIA_c3514_k24_u24_xHUMAN3x.species.tar.gz

Files will be expanded to ./database directory by default.

 > Note: The plasmid related results and the option need the new parsed database to work properly.
 > For users who downloaded the databases before 30th March 2015, we encourage you to download the corresponding
 > new parsed database (*.parsedGOTTCHA.dmp) to have a taste of new plasmid related feature. GOTTCHA v1.0 still supports old 
 > databases but plasmid relative results will always be shown as zero due to absence of plasmid information.
 > The new .dmp files are about 130MB size in gzipped format. Please unzipped the file and use it to
 > replace the old one. The general location of file is "gottcha/parsed_database_only/<DATABASE_NAME>.parsedGOTTCHA.dmp.gz". For
 > example, to download the new dmp file for species level database:
 >
 >      $ wget ftp://ftp.lanl.gov/public/genome/gottcha/parsed_database_only/GOTTCHA_BACTERIA_c3514_k24_u24_xHUMAN3x.species.parsedGOTTCHA.dmp.gz
 > 

Here is a list of the available pre-computed databases. Note that these databases are
also available in FASTA format at the FASTA/ directory:

  * GOTTCHA_BACTERIA_c3514_k24_u24.class.tar.gz (4.6GB)
  * GOTTCHA_BACTERIA_c3514_k24_u24.family.tar.gz (4.6GB)
  * GOTTCHA_BACTERIA_c3514_k24_u24.genus.tar.gz (4.5GB)
  * GOTTCHA_BACTERIA_c3514_k24_u24.order.tar.gz (4.6GB)
  * GOTTCHA_BACTERIA_c3514_k24_u24.phylum.tar.gz (4.6GB)
  * GOTTCHA_BACTERIA_c3514_k24_u24.species.tar.gz (4.3GB)
  * GOTTCHA_BACTERIA_c3514_k24_u24.strain.tar.gz (3.9GB)
  * GOTTCHA_BACTERIA_c3514_k24_u24_xHUMAN3x.genus.tar.gz (4.5GB)
  * GOTTCHA_BACTERIA_c3514_k24_u24_xHUMAN3x.species.tar.gz (4.3GB)
  * GOTTCHA_BACTERIA_c3514_k24_u24_xHUMAN3x.strain.tar.gz (3.8GB)
  * GOTTCHA_VIRUSES_c3498_k85_u24.genus.tar.gz (71MB)
  * GOTTCHA_VIRUSES_c3498_k85_u24.species.tar.gz (68MB)
  * GOTTCHA_VIRUSES_c3498_k85_u24.strain.tar.gz (68MB)
  * GOTTCHA_VIRUSES_c3498_k85_u24_xHUMAN3x.genus.tar.gz (71MB)
  * GOTTCHA_VIRUSES_c3498_k85_u24_xHUMAN3x.species.tar.gz (68MB)
  * GOTTCHA_VIRUSES_c3498_k85_u24_xHUMAN3x.strain.tar.gz (68MB)

Note: If you have any difficulty obtaining the databases, please contact Po-E Li <po-e@lanl.gov>.

-------------------------------------------------------------------
### Running GOTTCHA

The procedure includes 3 major steps: (1) split-trimming the input data, (2) mapping reads 
to a GOTTCHA database using BWA, and (3) profiling/filtering the results. These steps have
been wrapped into a sigle script called 'gottcha.pl'. User will need to provide a FASTQ file
as input and specify the location and name of the database.

Here is the general usage to run GOTTCHA:

 > $ bin/gottcha.pl -i \<FASTQ\> -d \<PATH/DATABASE_PREFIX\>

We provided a testing FASTQ file and example output in ./test. The following command
is an example that runs "test.fastq" through GOTTCHA using a species-level database with
8 threads:

        $ bin/gottcha.pl             \
             --threads 8             \
             --mode all              \
             --input test/test.fastq \
             --database database/GOTTCHA_BACTERIA_c3514_k24_u24_xHUMAN3x.species

In this case, we specify the output mode to "all" using "--mode all" option that
gives us two output files and all intermediate ouptuts stored in "test_temp"
directory. Both outputs are plain text files in tab-separated values format:
a summary table "test.gottcha.tsv" and a full information table 
"test.gottcha_full.tsv".

-------------------------------------------------------------------
### Interpreting Results

GOTTCHA reports profiling results in a neat summary table (*.gottcha.tsv) by default.
The tsv file will list the organism(s) at all taxonomic levels from STRAIN to PHYLUM,
their linear length, total bases mapped, linear depth of coverage, and the normalized 
linear depth of coverage. The linear depth of coverage (LINEAR_DOC) is used to calculate 
relative abundance of each organism or taxonomic name in the sample.

Summary table:

 Column            | Description
 ----------------- | -------------------------------------------------------------------
 LEVEL             | taxonomic rank
 NAME              | taxonomic name
 REL_ABUNDANCE     | relative abundance (equivalent to NORM_COV by default)
 LINEAR_LENGTH     | number of non-overlapping bases covering the signatures
 TOTAL_BP_MAPPED   | sum total of all hit lengths recruited to signatures
 HIT_COUNT         | number of hits recruited to signatures
 HIT_COUNT_PLASMID | number of hits recruited to signatures
 READ_COUNT        | number of reads recruited to signatures
 LINEAR_DOC        | linear depth-of-coverage (TOTAL_BP_MAPPED / LINEAR_LENGTH)
 NORM_COV          | normalized linear depth-of-coverage (LINEAR_DOC / SUM(LINEAR_DOC in certain level))

There are two report modes available. Other than a summary table, "full" report 
mode will report a table with more detail information from unfiltered results. 
The explanation of each column in the full report can be found in README_REPORT.md
The "all" report mode will keep all output files that were generated by each 
profiling step.

-------------------------------------------------------------------
### Visulizing Results using Krona

[Krona](http://sourceforge.net/p/krona/home/krona/) is an interactive browser that allows 
the exploration of hierarchical data with pie charts. Assuming you have Krona installed properly,
we are going to create a Krona chart from a text file listing abundance and lineages. 
You will find <PREFIX>.lineage.tsv file when you run gottcha.pl in "all" output mode.

Use 'ktImportText' and save the chart to "test.krona.html":

        $ ktImportText test_temp/test.lineage.tsv -o test.krona.html

------------------------------------------------------------------
## CITATION

Tracey Allen K. Freitas, Po-E Li, Matthew B. Scholz and Patrick S. G. Chain (2015) **Accurate read-based metagenome characterization using a hierarchical suite of unique signatures,** Nucleic Acids Research (DOI: 10.1093/nar/gkv180)

-------------------------------------------------------------------
## AUTHORS

Tracey Allen K. Freitas, Po-E Li, Matthew B. Scholz, Patrick S. G. Chain
Bioscience Division, Los Alamos National Laboratory, Los Alamos, NM 87545

-------------------------------------------------------------------
## ACKNOWLEDGEMENTS

We would like to thank Jason Gans for critical discussions on classification and machine 
learning techniques, and Shihai Feng for the generation of synthetic datasets.

This project is funded by U.S. Defense Threat Reduction Agency [R-00059-12-0 and R-00332-13-0 to P.S.G.C.].

-------------------------------------------------------------------
## CHANGE LOG
Version 1.0b:
> Bug fix for database inconsistency.

Version 1.0a:
> 1. Support multiple input files.
> 2. Add "--dumpSam" option to dump the mapping result.
> 3. Fix minor bugs.

Version 1.0:
> 1. Report the number of reads that hit to plasmids and provide an option to ignore them.
> 2. Report the number of READ_COUNT.

Version 0.9e:
> 1. Fix minor bugs
> 2. Amend the display of runtime information

Version 0.9d:
> 1. Splitrim allows lower-case bases in fastq file
> 2. Fix bugs that fail to specify output directory
> 3. Minor bug fix

Version 0.9c:
> 1. Fix FASTQ header compatibility
> 2. Provide more information while running GOTTCHA

Version 0.9b (05/12/2014):
> 1. Add '--stDir' option for pre-splitrimmed input file
> 1. Performance improvement

Version 0.9a (04/30/2014):
> 1. Provide bwaOpt option for user to use their own parameters to run BWA-MEM.
> 2. Use absolute path to run system calls
> 3. Provide 'relAbui' option to choose column to calculate relative abundance.

Version 0.9:
> 1. Initial release
