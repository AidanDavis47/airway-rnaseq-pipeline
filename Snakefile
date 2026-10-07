##Snake file with rules

#not sure if the comments will mess up the file, they dont
#visual code see this as a python file which is a little annoying since it highlights syntax errors

import csv #import csv

configfile: "config.yaml" #gets the setting from the config file, see file for more info

with open(config["samples"]) as f:
    SAMPLES = [row["sample"] for row in csv.DictReader(f)] #reads the samples from the samples csv (duh)
    
#want to add constraints to wild cards to make sure nothing weird happens
wildcard_constraints:
    s = r"SRR\d+", ##SRR because samples start with that
    r = "1|2"


#now implement thr rules for snake make (keep in mind snake make work in reverse kinda)


#tells snakemake the files wanted
rule all:
    input:
        expand("results/salmon/{s}/quant.sf", s=SAMPLES),
        "results/multiqc/multiqc_report.html"
        
#downloads the reads       
rule download:
    output: ##determines what the output files should be, located in the data file with {s} being the sample name
        "data/{s}_1.fastq.gz",
        "data/{s}_2.fastq.gz"
    params:
        n = config["reads_per_sample"] ##checks config file for how many reads to do per sample
    log: "logs/downloads/{s}.log" ##create a log for the run
    shell:
        "fastq-dump -X {params.n} --split-files --gzip -O data {wildcards.s} &> {log}" ##the command that will be ran 
  
#checks to make sure that there is no troublsome data, a qualtiy control step        
rule fastqc:
    input: "data/{s}_{r}.fastq.gz" ##gets the sample data
    output: #will output files to results, one html and one zip
        "results/fastqc/{s}_{r}_fastqc.html",
        "results/fastqc/{s}_{r}_fastqc.zip"
    log: "logs/fastqc/{s}_{r}.log" #create a log 
    shell:
        "fastqc {input} -o results/fastqc &> {log}" #command that will br ran, taking the input an outputting to results


#this builds that actual index      
rule salmon_index:
    input: config["transcripts"] #gets the transcripts to index
    output: directory(config["salmon_index"]) #determin the output location? looks like directory tells snake make we are looking for a folder and not just a single file, good to know
    threads: config["threads"] #determines how many threads will be used
    log: "logs/salmon_index.log" ##create a log
    shell:
        "salmon index -t {input} -i {output} --gencode -p {threads} &> {log}" ##shell command that will be ran

#determines the expression of each transcript
rule salmon_quant:
    input: ##gets the inputs to quantify, we have the r1 and r2 because the original files were split in half i beleive
        r1 = "data/{s}_1.fastq.gz",
        r2 = "data/{s}_2.fastq.gz",
        index = config["salmon_index"]
    output: "results/salmon/{s}/quant.sf" #determines what the output will be in results folder
    threads: config["threads"] #sets how many threads will be used
    log: "logs/salmon/{s}.log" #create log
    shell: #shell commands that will be ran
        "salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} "
        "-p {threads} -o results/salmon/{wildcards.s} &> {log}"

#this creats a summary as an html page       
rule multiqc:
    input: ##inputs, still a little iffy on this whole expand command, will take more of a look at it later
        expand("results/fastqc/{s}_{r}_fastqc.zip", s=SAMPLES, r=["1", "2"]),
        expand("results/salmon/{s}/quant.sf", s=SAMPLES)
    output: "results/multiqc/multiqc_report.html" #output which will be in html
    log: "logs/multiqc.log" #create log
    shell: #run the shell command
        "multiqc results/fastqc results/salmon -o results/multiqc --force &> {log}"