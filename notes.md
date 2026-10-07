XX = completed



using claude for help on starting this project


to open use ubuntu and type "code ." to get to wsl visual studio


notes on data:
508_1: looks normal except there is an overepresented sequence "ACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTT"
508_2: all looks normal
509_1: overepresented sequence "ACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTT"
509_2: overepresented sequence "GTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAA"

should see if the seq in 508_1 and 509_1 are the same and what that could mean
use stuff i learned from bioclass to see if they match, can do this in a different project folder



XX Look up the 8 samples' SRA run IDs from GSE52778 and save them in a samples.csv with donor and treatment columns



SAMPLES.CSV: got this from claude so be a bit wary, double check at a later point

using only the untreated and the dexamethasone pair to see a basic comparison between treated and untreated patients, less complicated than comparing treatments at the moment, want to make sure the treatment does anything in the first place


when done see if i can find published papers using this dataset to compare to see if i got it right


commands to keep in mind:

echo ""name of folder" >> .gitignore adds it to the gitignore
wget "link" pulls the data from the website
salmon index builds the index
salmon quant quantifies the index (both are self explanatory)



references:
GENCODE home page
GENCODE human release history
Ensembl blog, February 2026


mapping rates/top transcripts

ENST00000331825.11      871     712.552 23859.362387    10393.093
ENST00000387405.1       66      7.805   16721.377101    79.779
ENST00000678508.1       1778    1619.472        12433.137092    12309.029
ENST00000620041.5       830     671.571 11973.388913    4915.628
ENST00000361624.2       1542    1383.472        11485.308326    9713.655
ENST00000387409.1       66      7.805   9140.352774     43.609
ENST00000361851.1       207     71.985  8721.136983     383.781
ENST00000387400.1       73      7.356   8498.057201     38.215
ENST00000368719.9       434     277.420 6663.663340     1130.111
ENST00000362079.2       784     625.598 6446.937889     2465.575


snakemake notes:

an addition/variation of python

runs rules that are written in a Snakefile  (no extension)

Example rule:

rule bwa_map:
    input:  //inputs these files
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output: //creates these file
        "mapped_reads/{sample}.bam"
    shell: //shell commands that are in the rule set
        "bwa mem {input} | samtools view -Sb - > {output}"




runs via command line

activate enviro:  conda activate "enviro name"

dry run: snakemake -n
wet run: snake make --cores 4  ##the 4 is number of cores can change this number to increase or decrease



can use python scripts to enhance 

example:

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysam import VariantFile

quals = [record.qual for record in VariantFile(snakemake.input[0])]
plt.hist(quals)

plt.savefig(snakemake.output[0])



to make a histogram,

