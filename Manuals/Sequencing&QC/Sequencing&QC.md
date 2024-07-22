<img src="https://coursesandconferences.wellcomeconnectingscience.org/wp-content/themes/wcc_courses_and_conferences/dist/assets/svg/logo.svg" width="200" height="200">

[<<< Go back to Manual Contents Page](https://github.com/WCSCourses/GenEpiLAC2024/blob/main/Manuals/Manual_main.md)

<br>

# Module 1 - Sequencing & QC - Costa Rica 2024 <!-- omit in toc -->

## Table of contents <!-- omit in toc -->
- [Introduction](#introduction)
  - [Brief recap](#brief-recap)
  - [Commonly used file formats for NGS data](#commonly-used-file-formats-for-NGS-data)
    - [FASTQ](#fastq)
    - [FASTA](#fasta)
    - [SAM and BAM](#sam-and-bam)
  - [Quality control for FastQ data](#quality-control-for-FastQ-data)
- [Practical Exercise](#practical-exercise)
  - [QC of short data](#QC-of-short-data)
  - [QC of long data](#QC-of-long-data)
  - [Bonus](#bonus)
  - [Short quiz on Linux](#short-quiz-on-Linux)

## Module Overview and Aims

In this module, we are going to get familiar with several common file formats used for sequence data and then we are going to perform some quality control (QC) on some FASTQ-formatted sequence data obtained from both short and long reads sequencing technologies.

# [Introduction](#introduction)

## [Brief recap](#brief-recap)

Before we start, during the whole week, you will be using most of the concepts you´ve learnt in the pre-course "Introduction to Linux for biologists". Below you can find a list of frequently used commands for you to remember them: 

| Command | What it does |
| ------ | ------------ |
| ls | Lists the contents of the current directory  |
| mkdir  | Makes a new directory  |
| mv | Moves or renames a file  |
| cp  | Copies a file  |
| rm | Removes a file  |
| mkdir  | Makes a new directory  |
| cat | Concatenates files  |
| more  | Displays the contents of a file one page at a time  |
| head | Displays the first ten lines of a file  |
| tail  | Displays the last ten lines of a file  |
| cd | Changes current working directory  |
| pwd	 | Prints working directory  |
| find  | Finds files matching an expression  |
| grep | Searches a file for patterns  |
| wc  | Counts the lines, words, characters, and bytes in a file  |
| kill | Stops a process  |
| jobs | Lists the processes that are running  |

```
Command line shortcuts 
●	Up/Down arrows: Previous commands 
●	!!: Reruns previous command 
●	Tab: Auto complete 
●	Tab+Tab: All available options 
●	Ctrl+a: Move cursor to start of line 
●	Ctrl+e: Move cursor to end of line 
●	Alt+: Alternates between terminals 
●	Ctrl+l: Clear screen (or Command+k on Mac) 
●	Ctrl+c: Terminates the running program 
●	Ctrl+z: Suspends the running program 
●	Ctrl+w: Removes a previous word 
●	Ctrl+d: Logout 
●	Ctrl+d(in a command): Removes a character 
●	Ctrl+u: Removes till the beginning 
```
## [Commonly used file formats for NGS data](#commonly-used-file-formats-for-NGS-data)

### FASTQ

The widely used FASTA file format has the great advantage of simplicity. However, this simplicity can be restrictive if we want to include additional data/metadata in addition to the sequence.

Given the non-negligible error rates of NGS technologies, often we need to accompany our sequence data with quality scores that estimate our confidence in the accuracy of the sequence data. As we will see later, this allows us to perform quality control checks and filter-out poor-quality data before performing analyses.

**FASTQ** is a simple text-based format that allows us to include quality scores. A single sequence is represented by four lines of text:

    @ERR8261968.1 1 length=97
    ACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTA
    +ERR8261968.1 1 length=97
    CCCCCFDDFFFFGGGGGGGGGGHHHHHHHHHHHGGGGHHHHHHHHHHHHHHHGHHGHHIIHHGGGGGGHHHHHHHHHHHHHHHHHHHGGGHHHHHHH

- The first line is a 'header' containing a unique identifier for the sequence and, optionally, further description.
- The second line contains the actual nucleotide sequence.
- The third line is redundant  and can be safely ignored. Sometimes it simply repeats the first line. Sometimes it is blank or just contains a '+' character.
- The fourth line contains a string of characters that encode quality scores for each nucleotide in the sequence on [ASCII code](https://en.wikipedia.org/wiki/ASCII). 

Each single character encodes a score, typically a number between 0 and 40; this score is encoded by a single character.

| Character | ASCII | FASTQ quality score (ASCII – 33) 
| --|--|--
| ! | 33 | 0
| “ | 34 | 1
| # | 35 | 2
| $ | 36 | 3
| % | 37 | 4
| ... | ... | ...
| C | 67 | 34
| D | 68 | 35
| E | 69 | 36
| F | 70 | 37
| G | 71 | 38
| H | 72 | 39
|40 | 73 | 40

So, in the example above, we can see that most of the positions within the 97-nucleotide sequence have scores in the high 30s, which indicates a high degree of confidence in their accuracy.
- A score of 30 denotes a 1 in 1000 chance of an error, i.e. 99.9% accuracy.
- A score of 40 denotes a 1 in 10,000 chance of an error, i.e. 99.99% accuracy.

You can read more about the FastQ file format and quality scores here:

Cock, P. J., Fields, C. J., Goto, N., Heuer, M. L., & Rice, P. M. (2010). The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. *Nucleic Acids Research*, **38**, 1767–1771. https://doi.org/10.1093/nar/gkp1137.

### FASTA

Among the most common and simplest file formats for representing nucleotide sequences is FASTA.  Essentially, each sequence is represented by a 'header' line that begins with a '>', followed by lines containing the actual nucleotide sequence. By convention, the first 'word' in the header line is a unique identifier, which is usually the accession number. Consider this portion of the whole genome of *Staphylococcus aureus* MRSA252 as an example of a FASTA-formatted nucleotide sequence:

```
  >NC_002952.2 Staphylococcus aureus subsp. aureus MRSA252, complete sequence
  CGATTAAAGATAGAAATACACGATGCGAGCAATCAAATTTCATAACATCACCATGAGTTTGATCCAAAGCATGAGTGTTT
  ACAATGTTTGAATACCTTATACAGTTCTTATACATACTTTATAAATTATTTCCCAAGCTGTTTTGATACACACACTAACA
  GATACTCTATAGAAGGAAAAGTTATCCACTTATGCACACTTATACTTTTTAGAATTGTGGATAATTAGAAATTACACACA
  AAGTTATACTATTTTTAGCAACATATTCACAGGTATTTGACATATAGAGAACTGAAAAAGTATAATTGTGTGGATAAGTC
  GTCCAACTCATGATTTTATAAGGATTTATTTATTGATATTTACATAAAAATACTGTGCATAACTAATAAGCAGGATAAAG
  TTATCCACCGATTGTTATTAACTTGTGGATAATTATTAACATGGTGTGTTTAGAAGTTATCCACGGTTGTTATTTTTGTG
  TATAACTTAAAAATTTAAGAAAGATGGAGTAAATTTATGTCGGAAAAAGAAATTTGGGAAAAAGTGCTTGAAATTGCTCA
  AGAAAAATTATCAGCTGTAAGTTACTCAACTTTCCTAAAAGATACTGAGCTTTACACGATCAAAGATGGTGAAGCTATCG
  TATTATCGAGTATTCCTTTTAATGCAAATTGGTTAAATCAACAATATGCTGAAATTATCCAAGCAATCTTATTTGATGTT
  GTAGGCTATGAAGTAAAACCTCACTTTATTACTACTGAAGAATTAGCAAATTATAGTAATAATGAAACTGCTACTCCAAA
  AGAAGCAACAAAACCTTCTACTGAAACAACTGAGGATAATCATGTGCTTGGTAGAGAGCAATTCAATGCCCATAACACAT
  TTGACACTTTTGTAATCGGACCTGGTAACCGCTTCCCACATGCAGCGAGTTTAGCTGTGGCCGAAGCACCAGCCAAAGCG
  TACAATCCATTATTTATCTATGGAGGTGTTGGTTTAGGAAAAACCCATTTAATGCATGCCATTGGTCATCATGTTTTAGA
  TAATAATCCAGATGCCAAAGTGATTTACACATCAAGTGAAAAATTCACAAATGAATTTATTAAATCAATTCGTGATAACG
```

- The first line begins with '>' indicating that it is the header line.
- This is immediately followed by 'NC_002952.2', which is the accession number for [this sequence in the NCBI Assembly database](https://www.ncbi.nlm.nih.gov/assembly/GCF_000011505.1).
- Then follows the actual nucleotide sequence, split over several lines.

It is very common to combine multiple sequences into a single multi-FASTA file like this:

```
>ctxB1
ATGATTAAATTAAAATTTGGTGTTTTTTTTACAGTTTTACTATCTTCAGCATATGCACAT
GGAACACCTCAAAATATTACTGATTTGTGTGCAGAATACCACAACACACAAATACATACG
CTAAATGATAAGATATTTTCGTATACAGAATCTCTAGCTGGAAAAAGAGAGATGGCTATC
ATTACTTTTAAGAATGGTGCAACTTTTCAAGTAGAAGTACCAGGTAGTCAACATATAGAT
TCACAAAAAAAAGCGATTGAAAGGATGAAGGATACCCTGAGGATTGCATATCTTACTGAA
GCTAAAGTCGAAAAGTTATGTGTATGGAATAATAAAACGCCTCATGCGATTGCCGCAATT
AGTATGGCAAATTAA
>ctxB10
ATGATTAAATTAAAATTTGGTGTTTTTTTTACAGTTTTACTATCTTCAGCATATGCACAT
GGAACACCTCAAAATATTACTGATTTGTGTGCAGAATACCCCAACACACAAATATATACG
CTAAATGATAAGATATTTTCGTATACAGAATCTCTAGCTGGAAAAAGAGAGATGGCTATC
ATTACTTTTAAGAATGGTGCAATTTTTCAAGTAGAAGTACCAGGTAGTCAACATATAGAT
TCACAAAAAAAAGCGATTGAAAGGATGAAGGATACCCTGAGGATTGCATATCTTACTGAA
GCTAAAGTCGAAAAGTTATGTGTATGGAATAATAAAACGCCTCATGCGATTGCCGCAATT
AGTATGGCAAATTAA
>ctxB11
ATGATTAAATTAAAATTTGGTGTTTTTTTTACAGTTTTACTATCTTCAGCATATGCACAT
GGAACACCTCAAAATATTACTGATTTGTGTGCAGAATACCCCAACACACAAATACATACG
CTAAATGATAAGATATTTTCGTATACAGAATCTCTAGCTGGAAAAAGAGAGATGGCTATC
ATTACTTTTAAGAATGGTGCAACTTTTCAAGTAGAAGTACCAGGTAGTCAACATATAGAT
TCACAAAAAAAAGCGATTGAAAGGATGAAGGATACCCTGAGGATTGCATATCTTACTGAA
GCTAAAGTCGAAAAGTTATGTGTATGGAATAATAAAACGCCTCATGCGATTGCCGCAATT
AGTATGGCAAATTAA
>ctxB12
ATGATTAAATTAAAATTTGGTGTTTTTTTTACAGTTTTACTATCTTCAGCATATGCACAT
GGAACACCTCAAAATATTACTGATTTGTGTGCAGAATACCACAACGCACAAATATATACG
CTAAATGATAAGATATTGTCGTATACAGAATCTCTAGCTGGAAACAGAGAGATGGCTATC
ATTACTTTTAAGAATGGTGCAACTTTTCAAGTAGAAGTACCAGGTAGTCAACATATAGAT
TCACAAAAAAAAGCGATTGAAAGGATGAAGGATACCCTGAGGATTGCATATCTTACTGAA
GCTAAAGTCGAAAAGTTATGTGTATGGAATAATAAAACGCCTCATGCGATTGCCGCAATT
AGTATGGCAAATTAA
>ctxB2
ATGATTAAATTAAAATTTGGTGTTTTTTTTACAGTTTTACTATCTTCAGCATATGCACAT
GGAACACCTCAAAATATTACTGATTTGTGTGCAGAATACCACAACACACAAATACATACG
CTAAATGATAAGATATTGTCGTATACAGAATCTCTAGCTGGAAAAAGAGAGATGGCTATC
ATTACTTTTAAGAATGGTGCAACTTTTCAAGTAGAAGTACCAGGTAGTCAACATATAGAT
TCACAAAAAAAAGCGATTGAAAGGATGAAGGATACCCTGAGGATTGCATATCTTACTGAA
GCTAAAGTCGAAAAGTTATGTGTATGGAATAATAAAACGCCTCATGCGATTGCCGCAATT
AGTATGGCAAATTAA
```

If you want a more detailed history of the FASTA file format, then you could take a look at the Wikipedia page here: https://en.wikipedia.org/wiki/FASTA_format.

### SAM and BAM

A **SAM file** (usually named *.sam) is used to represent aligned sequences. It is particularly useful for storing the results of genomic or transcriptomic sequence reads aligned against a reference genome sequence. The **BAM file** format is a compressed form of SAM. This has the disadvantage that it is not readable by a human but has the advantage of being smaller than the corresponding SAM file and thus easier to share and copy between locations.

You can read about SAM and BAM formats here:
 - Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, **25**, 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
-  [https://samtools.github.io/hts-specs/SAMv1.pdf](https://samtools.github.io/hts-specs/SAMv1.pdf).

We can view BAM files graphically using specialised genome browsers software such as (we will be working with one of them in the Genome Visualization Tools module):
- [IGV](https://igv.org/)
- [Tablet](https://ics.hutton.ac.uk/tablet/)
- [Artemis / BAMview](http://sanger-pathogens.github.io/Artemis/Artemis/) 

An example **SAM file** is shown in the image below, along with the 11 mandatory fields of this standard:

![SAMformat](https://user-images.githubusercontent.com/65819144/230166077-264a6161-ec8c-47d7-87b7-8e073080a548.jpg)

(_Image from [Data Wrangling and Processing for Genomics](https://datacarpentry.org/wrangling-genomics)_)

The meaning of each of these mandatory fields are detailed below:

![SAMformat_fields](https://user-images.githubusercontent.com/65819144/230166180-fb086380-413a-4eef-8781-d0ba2b7e215b.jpg)

## [Quality control for FastQ data](#quality-control-for-FastQ-data)

Modern high throughput sequencers can generate hundreds of millions of sequences in a single run. Before analysing this sequence to draw biological conclusions you should always perform some simple quality control checks to ensure that the raw data looks good and there are no problems or biases in your data which may affect how you can usefully use it.

Most sequencers will generate a QC report as part of their analysis pipeline, but this is usually only focused on identifying problems which were generated by the sequencer itself. To check the quality of the reads generated, we will be using the software [**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) which will provide a QC report that can spot problems which originate either in the sequencer or in the starting library material. 

FastQC can be run in one of two modes. It can either run as a stand alone interactive application for the immediate analysis of small numbers of FastQ files, or it can be run in a non-interactive mode where it would be suitable for integrating into a larger analysis pipeline for the systematic processing of large numbers of files.

![fastqc](https://user-images.githubusercontent.com/65819144/230931477-4ed32628-11d5-4b1a-bfe6-70633185eeec.png)

FastQC will highlight any areas where this library looks unusual and where you should take a closer look. The program is not tied to any specific type of sequencing technique and can be used to look at libraries coming from a large number of different experiment types (Genomic Sequencing, ChIP-Seq, RNA-Seq, BS-Seq etc etc).

Different softwares for assesing the quality of reads have been developed. Another useful software for checking the quality of long-reads is **NanoPlot**. NanoPlot assess characteristics relevant to long-read genome sequencing and produces different number of plots depending on the data. A detailed view can be seen on [here](https://github.com/wdecoster/NanoPlot#plots-generated).

It is very common to have some quality metrics fail after running FastQC, and this may or may not be a problem for your downstream application. But don't worry there are softwares developed to filter poor quality reads and trim poor quality bases or adapters from our samples. In this module we will be working with [**TrimGalore**](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for short reads and with [**Filtlong**](https://github.com/rrwick/Filtlong) for Oxford Nanopore reads.

>**Note**: **Porechop** is a tool for finding and removing adapters from Oxford Nanopore reads. Adapters on the ends of reads are trimmed off, and when a read has an adapter in its middle, it is treated as chimeric and chopped into separate reads. We won't be using this software in this module but you could try it with your data! More information [here](https://github.com/rrwick/Porechop).

# [Practical Exercise](#practical-exercise)

For this module, we will be working with files which are inside the folder "Module_1_Sequencing_QC". Do you remember how to change to this directory? Give it a try!

## [QC of short data](#QC-of-short-data)

Inside the directory "Module_1_Sequencing_QC", we have 10 compressed fastq files as fastq.gz. You can check that the files are there. Do you remember how to do that?

We will work with files:
```
ARIMSS995-11_1.fastq.gz
ARIMSS995-11_2.fastq.gz
```

As we saw in the introduction, this is the format we will get from an Illumina paired-end run. To have a better look at their structure, we will uncompress them:
```
gzip -d ARIMSS995-11_1.fastq.gz
```
Wait until the command line is shown again on the screen and then type:
```
gzip -d ARIMSS995-11_2.fastq.gz
```

The files should be decompressed now and ready to work. Let’s have a look at one of the fastq files:

```
head ARIMSS995-11_1.fastq
```

You should see something like this:

```
@M03346:10:000000000-AL260:1:1101:28195:9448 1:N:0:3
CCTCATCAGCCGCATTGCGACCAACTTCTGGATTAGCGCCAGCGCCCAGTCCTTTG
GTGATACCGCTACCGATTTGAATCGTCTGTCCAACCGCTGTTTTACGCAGCGCTTGT
GCATCGGTATTTACCGCGAAGAATTCAACACCTTCAATACGCTCGCGCACCATGTGT
TCAACAGCATTACCGCCGCCGCCGCCGACGCCGATGACTTTAATCACCGCGTCATT
GGTAAGTTCCATTGGTTCAAACATA
+
ABCCCFFFFFFCGGGGGGGGGGGHGHHHHHGHHHHHHGGGGGGGGGGGHHHHH
HHHHHGHHHHHGGGGGHGGGGGHHHHHHHHHGHHHHHHGGGGGFHHHHHGGG
GGGGGGGHHHHHHGGFHGHHHHHGGGGGGGHHHHHHHHGGHHHHHHGHEGGG
GGGGGGGGGGGGGGGGGGGGGGGFFFFFFFFFFFFFFFFFCFFFFFFAFFFFFFFFF
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFB
@M03346:10:000000000-AL260:1:1101:5008:22354 1:N:0:3
TCTGTGTGGGGCGCTGATTATAGACGTTACGGCTGCAATCTTCATCTGTTTGTATGG
ATATCGCTTCTCGGGTATTGTTCTCACCCGGTGCTATCGCCGTTTTTCCTGCAAGGT
TTTATTTTCATGTCACAGTCTACATCCGTGCTTCGTCGTAATGGATTTACTTTTAAACA
GTTTTTTGTTGCTCACGATCGCTGTGCGATGAAAGTGGGAACGGATGGCATTTTGC
TGGGCGCATGGGCACCGGTGGG
+
AAABAFDFAAB>EAEGEGGFGGHHHFFFHAFGGGGGGHHHHHHGHGHHHHHHHHH
HHBFFHGGFHEEGGGEG?/GHHHHGHFFHH1EEEEGHFHFGE?EBDGFGHDHFFHH
HGGGHGGHHHHGDHHH?<<GHFHFHHDGHFDGGGGHHDGDGGHFEHHEHHHDDD
HGHGFCCFCGHHHGGGGGGFGFFFE?E.ADFGG/0.CEFFGGFFF/EEF.AD.DBFFFFF
FFF/;/FF.;@-ADFE99..;B--;-
```

**Remembering the fastq file format, how many reads do you have represented in the previous box?**

Now, let’s check how many reads we have in ARIMSS995-11_1.fastq.gz and ARIMSS995-11_2.fastq.gz, and double check we have the same number in the R1 and R2 files. We can use the `wc` command with the -l (letter elle) option to count the number of lines in the file:

```
wc -l ARIMSS995-11_1.fastq
```
**How many paired reads are in the sample?**
>Hint: Think about how many lines a single read takes up in FASTQ format.

**Does the number of reads in R1 and R2 files match?**

Let's examine the quality of these sequence data using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). We will only use the two fastq files we have been working with and two additional ones located in the same directory named 'Module_1_Sequencing_QC'. The file names are: **untrimmed_1.fastq.gz** and **untrimmed_2.fastq.gz**. 

First, for the sake of simplicity, let’s compress again the two fastq files belonging to sample ARIMSS995-11:

```
gzip ARIMSS995-11_*.fastq
```

Going back to FastqQC, we can launch the graphical interface by simply executing ``fastqc`` on the Terminal command line. However, it is often more convenient to use the software in the command-line mode. Execute the following command in the Terminal:

    fastqc *.fastq.gz

> With this command we will be running fastqc in all the files that have the `.fastq.gz` extension. Nonetheless we will be focusing only in samples ARIMSS995-11 and untrimmed.

> **Note**: instead of compressing the files again with gzip, you could also use the command `fastqc *.fastq*`. By adding a star symbol after `.fastq`, fastqc will run in all files that include `.fastq` as part of the extension (zipped and unzipped)
         

You will see some messages like this on your screen:

    Started analysis of ARIMSS995-11_1.fastq.gz
    Approx 5% complete for ARIMSS995-11_1.fastq.gz
    Approx 10% complete for ARIMSS995-11_1.fastq.gz
    ...
    Approx 90% complete for ARIMSS995-11_1.fastq.gz
    Approx 95% complete for ARIMSS995-11_1.fastq.gz
    Analysis complete for ARIMSS995-11_1.fastq.gz
    Started analysis of ARIMSS995-11_2.fastq.gz
    Approx 5% complete for ARIMSS995-11_2.fastq.gz
    Approx 10% complete for ARIMSS995-11_2.fastq.gz
    Approx 15% complete for ARIMSS995-11_2.fastq.gz
    ...
    Approx 90% complete for ARIMSS995-11_2.fastq.gz
    Approx 95% complete for ARIMSS995-11_2.fastq.gz
    Analysis complete for ARIMSS995-11_2.fastq.gz


Now, execute the command ``ls -lh`` and you should see some new files have appeared. You should see something like this:

    
    total 909M
    -rw-rw-r-- 1 manager manager 773K Apr 10 16:45 ARIMSS995-11_1_fastqc.html
    -rw-rw-r-- 1 manager manager 513K Apr 10 16:45 ARIMSS995-11_1_fastqc.zip
    -rwxrwx--- 1 manager manager  72K Mar 23 10:09 ARIMSS995-11_1.fastq.gz
    -rw-rw-r-- 1 manager manager 789K Apr 10 16:45 ARIMSS995-11_2_fastqc.html
    -rw-rw-r-- 1 manager manager 525K Apr 10 16:45 ARIMSS995-11_2_fastqc.zip
    -rwxrwx--- 1 manager manager  67K Mar 23 10:09 ARIMSS995-11_2.fastq.gz
    -rw-rw-r-- 1 manager manager 602K Apr 10 16:45 ERR1009125_1_fastqc.html
    -rw-rw-r-- 1 manager manager 289K Apr 10 16:45 ERR1009125_1_fastqc.zip
    -rwxrwx--- 1 manager manager 101M Mar 23 10:10 ERR1009125_1.fastq.gz
    -rw-rw-r-- 1 manager manager 598K Apr 10 16:45 ERR1009125_2_fastqc.html
    -rw-rw-r-- 1 manager manager 291K Apr 10 16:45 ERR1009125_2_fastqc.zip
    -rwxrwx--- 1 manager manager 106M Mar 23 10:09 ERR1009125_2.fastq.gz
    -rw-rw-r-- 1 manager manager 605K Apr 10 16:46 ERR1009142_1_fastqc.html
    -rw-rw-r-- 1 manager manager 297K Apr 10 16:46 ERR1009142_1_fastqc.zip
    -rwxrwx--- 1 manager manager  79M Mar 23 10:09 ERR1009142_1.fastq.gz
    -rw-rw-r-- 1 manager manager 600K Apr 10 16:46 ERR1009142_2_fastqc.html
    -rw-rw-r-- 1 manager manager 295K Apr 10 16:46 ERR1009142_2_fastqc.zip
    -rwxrwx--- 1 manager manager  84M Mar 23 10:09 ERR1009142_2.fastq.gz
    -rw-rw-r-- 1 manager manager 563K Apr 10 16:46 ERR200457_1_fastqc.html
    -rw-rw-r-- 1 manager manager 327K Apr 10 16:46 ERR200457_1_fastqc.zip
    -rwxrwx--- 1 manager manager 144M Mar 23 10:10 ERR200457_1.fastq.gz
    -rw-rw-r-- 1 manager manager 564K Apr 10 16:46 ERR200457_2_fastqc.html
    -rw-rw-r-- 1 manager manager 322K Apr 10 16:46 ERR200457_2_fastqc.zip
    -rwxrwx--- 1 manager manager 148M Mar 23 10:09 ERR200457_2.fastq.gz
    -rwxrwx--- 1 manager manager  26K Mar 23 10:10 Metadata.csv
    -rwxrwx--- 1 manager manager 4.7M Mar 23 10:09 reference_Ss046.fasta
    -rwxrwx--- 1 manager manager  33K Mar 23 10:09 Ssonei.txt
    -rw-rw-r-- 1 manager manager 656K Apr 10 16:46 untrimmed_1_fastqc.html
    -rw-rw-r-- 1 manager manager 378K Apr 10 16:46 untrimmed_1_fastqc.zip
    -rwxrwx--- 1 manager manager 109M Mar 23 10:09 untrimmed_1.fastq.gz
    -rw-rw-r-- 1 manager manager 664K Apr 10 16:46 untrimmed_2_fastqc.html
    -rw-rw-r-- 1 manager manager 391K Apr 10 16:46 untrimmed_2_fastqc.zip
    -rwxrwx--- 1 manager manager 126M Mar 23 10:09 untrimmed_2.fastq.gz


We are most interested in the HTML files, which contain the FastQC reports for our fastq files. Let's open the HTML files generated for the two samples we are working with.

Use the following command as an example:

    firefox ARIMSS995-11_*.html &

> **Note**: Don't worry if error messages appear on the command line, just hit enter to get back the prompt.

You should then see something like this:

![fastqc1](https://user-images.githubusercontent.com/65819144/230941522-3129125c-9423-4475-b617-1cc868f93786.jpg)


There is a lot of QC information in these reports. Feel free to explore these in your own time and take a look at the [FastQC homepage](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) where you can find the [explanation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) to each report FastQC is generating and we recommend you to watch the tutorial video at http://www.youtube.com/watch?v=bz93ReOv87Y.

For now, we are just going to look at
- Basic statistics
- Per-base sequence quality 
- Adapter content

**Analysing FastQC reports, how many reads are there? Does this answer match your previous answer for sample ARIMSS995-11 (based on ``wc`` command)?**

**With respect to quality scores, which of the two files has better-quality data: *ARIMSS995-11_1.fastq.gz* or *ARIMSS995-11_2.fastq.gz*?**

**Besides sequence differences, are there any other differences related to run quality between the two samples we analysed (ARIMSS995-11 and untrimmed)?**

**Are these datasets contaminated with any Illumina sequencing adapter oligonucleotides?**

Now, we are going to look at how we can remove poor data and adapter contamination by trimming and filtering. We will use [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) by executing the following command on the Terminal:

We will need to do some minor trimming (quality 25, length 50) as well as checking/removing Illumina adapter sequences:
```
trim_galore -q 25 --length 50 --paired --illumina --fastqc ARIMSS995-11_1.fastq.gz ARIMSS995-11_2.fastq.gz
```

-q 25 = trim the 3’ end of the reads – remove nucleotides less than Phred Quality 25

--length 50 = after adapter and quality trimming, remove reads less than length 50bp

--paired = the names of the paired FASTQ files to analyses in order

--illumina = adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter AGATCGGAAGAGC instead of the default auto-detection of adapter sequence

--fastqc = run FastQC in the default mode on the FastQ file once trimming is complete

Once trim_galore has finished, check the outputs. You should see that two new FASTQ (.fq) files have been created by trim_galore:

ARIMSS995-11_1_val_1.fq.gz

ARIMSS995-11_2_val_2.fq.gz

**How many paired reads are left in the sample after trimming? Compare them with the fastq files obtained from the sequencer**

If you have extra time, you could try trimming the 'untrimmed' sample!

## [QC of long data](#QC-of-long-data)

The file we will working with for this part of the exercise is in a Google Drive link, so let's download it to our VM:

[https://drive.google.com/file/d/1QqNDYqGfQbkcJrP4nchXvdL94BCvAUQf/view?usp=sharing](https://drive.google.com/file/d/1QqNDYqGfQbkcJrP4nchXvdL94BCvAUQf/view?usp=sharing)


Now we want to move that file into our current working directory (Module_1_Sequencing_QC):

```
mv /home/manager/Downloads/CTMA_1441_longds.fastq.gz .
```

Remember to check that you have moved correctly this file!

Let's check the quality of the reads we have downloaded:

```
NanoPlot -o nanoplot/CTMA_1441_longds --tsv_stats --info_in_report --N50 --no_static --prefix CTMA_1441_longds_ --fastq CTMA_1441_longds.fastq.gz
```

Let's take a look at the outputs generated by Nanoplot:

```
ls nanoplot/CTMA_1441_longds
```

You'll see different .html and .txt files with different plots and information which are summarized in the report: 

```
firefox nanoplot/CTMA_1441_longds/CTMA_1441_longds_NanoPlot-report.html &
```

**How many reads are there in CTMA_1441_longds.fastq.gz?**
**What is the median read length, N50 and median quality score?**

Now, let's improve the quality of our long reads using Filtlong. We will remove the worst reads (short length and low accuracy) with --keep_percent 90 and just keep those reads which are longer than 1000 kb (--min-length 1000):

```
filtlong --min_length 1000 --keep_percent 90 CTMA_1441_longds.fastq.gz | gzip > CTMA_1441_longds_filt.fastq.gz
```

Adapt the previous code to run NanoPlot on the newly created long-read filtered file (CTMA_1441_longds_filt.fastq.gz)

**Examine the results for the filtered file and try to answer the questions you answered for the unfiltered file. Do you find any difference among these results?**

## [Bonus](#bonus)

#### What if you had a large number of FastQC reports to analyze?

Multiqc (https://multiqc.info/) is a tool that summarizes different types of NGS reports (not just FastQC).

First, install the tool:

    pip install multiqc

In the folder with reports, run:

    multiqc .

See summarized report in a browser:

    firefox multiqc_report.html &

## [Short quiz on Linux](#short-quiz-on-Linux)

If after this module you have time and you feel you need to practice more some of the commands you've learnt in the pre-course "Introduction to Linux for biologists", here's a short quiz to do so: 

1. Extract the first 15 lines from file “reference_Ss046.fasta” and save the output into “output.fa” 
2. How many files are there in the Module_1_Sequencing_QC directory? 
3. Create a new file called “Module_1_Sequencing_QC.txt” which contains a list of the file names present in the Module_1_Sequencing_QC directory.
4. Get the list of countries that contributed *S. sonnei* genome sequences for the publication previously mentioned and save it to “countries.txt” 
5. Extract the Assembly accessions of the sequences in “Ssonei.txt” and save it to “assemblies.txt”. How many are there?
6. Given a file with different sequences of Latin America, how do you count the ones submitted by Chile? (we know they should have the word “CHI” in the line).
7.  .............is the command used to create a new directory.
8. Command used to create an empty file.
9. "cat" is the command used to ...................
10. ............. command is used to count the total number of lines, words and characters in a file. 
11. Which command would you use to know the location of your current working directory?

[<<< Go back to Manual Contents Page](https://github.com/WCSCourses/GenEpiLAC2024/blob/main/Manuals/Manual_main.md)

<br>
