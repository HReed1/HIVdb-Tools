# HIVdb Tools

The main purpose of these programs is to organize and extract relevant information generated from HIV-1 PR, RT, and IN amplicons. These programs were designed to work specifically with an HIV1-Resistance study for Dr. Jeanne Jordan's Lab, working under the DC Cohort sponsored by NIAID-NIH(U01 AI69503-03S2).

### Prerequisites

Perl-5.22.1 or later (https://www.perl.org/get.html)
Python 3.5.2 or later (https://www.python.org/downloads/release/python-352/)
Linux OS: Written on Ubuntu 16.04 LTS (https://www.ubuntu.com/download/desktop)

## hivdb_organize.py
This program separates PR, RT, and IN gene regions into their own FASTA files for use in further analysis.

### Getting Started:
First, you will need to have a multi-fasta for the consensus sequeneces of your samples.

```
>Sample1
GCATGCATGCT
>Sample2
ACGTAGCTAGC
>Sample3
GTAGCTAGTGT
...
```

Next, you will have to make a "sampleList.txt" file that holds all of your sequence headers without the ">". To do this you will need to use this Perl one-liner to quickly make the file:

```
perl -ne 'if(/^>(.+)/){print"$1\n";}' my_cons.fasta > sampleList.txt

```
### Usage:

```
python3 ./hivdb_organize.py
```

### Output:

You will find multiple multi-fasta files created from this program if the respective gene regions are found:
```
AllPR.fasta
AllRT.fasta
AllIN.fasta
AllPR_RT.fasta
AllRT_IN.fasta
AllReigons.fasta
```

##hivdb_extract.py
Using a [HyDRA](https://hydra.canada.ca/)  Mutation-Report file, this program extracts AA mutations and separates them into individual files for input into the [Stanford HIVdb Program](https://hivdb.stanford.edu/hivdb/by-mutations/).

### Getting Started:
To begin, you will need to have your Mutation-Report file (.csv) ready. Make sure that the spreadsheet is saved in the .csv format:

```
PL-16,RT,NNRTI,No,V,108,I,15.03,1071
PL-93,PR,Accessory,No,L,10,F,15.41,1382
PL-73,RT,NRTI,Yes,M,184,V,15.42,10459
PL-77,RT,NRTI,Yes,M,184,V,15.59,5473
PL-16,RT,NNRTI,No,L,100,V,16.13,998

```

### Usage
```
python3 ./hivdb_extract.py
```

### Output
After running, there will be a Variant file (.txt) made for each sample containing a relevant mutation. To evaluate in the Stanford HIVdb, you will have to copy and paste each section into their respective text box:

Variant_1.txt
```
RT
V90I
S68G

PR
A71T

IN
L74I
```

## Authors

**Harrison Reed** - *Initial work* - [HReed1](https://github.com/HReed1)

