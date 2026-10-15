##this will be used to create a dataset that will be used to do quick automated test

##has 200 genes and their transcripts, also 20000 read pairs per sample that map to the genes

#imports
import csv, gzip, random, shutil, subprocess
from pathlib import Path

sig_genes = 100 # the most significant genes
ctrl_genes = 100 # genes that dont change much that are used as controls
max_reads = 20000 # max amount of reads per sample
must_genes = {"FKBP5", "TSC22D3", "PER1", "CRISPLD2", "DUSP1"} # these are the genes that we need, know steroids
samples = { #samples that we will be pulling the genes from
    "SRR1039508": ("N61311", "untreated"),
    "SRR1039509": ("N61311", "dexamethasone"),
    "SRR1039512": ("N052611", "untreated"),
    "SRR1039513": ("N052611", "dexamethasone")
}

out = Path("tests/data") #output to a tests folder
temp = Path("tmp_testdata") #temporary test data
(out / "data").mkdir(parents=True, exist_ok=True) #make a new folder for data
(out / "ref").mkdir(parents=True, exist_ok=True) #make a new folder for ref
temp.mkdir(exist_ok=True) #make a temp folder
random.seed(42) #psudo random

#choose genes: looking for significant ones and control genes and the must genes which are known steroid genes
with open("results/deseq2/deseq2_results.csv") as f:
    rows = [r for r in csv.DictReader(f) if r["padj"] != "NA"]
sig = [r["gene"] for r in rows[:sig_genes]] ##get significant genes
steady = [r["gene"] for r in rows if float(r["padj"]) > 0.5 and float(r["baseMean"]) > 200] #find steady genes that dont fluctuate will be used for the control genes
controls = random.sample(steady, min(ctrl_genes, len(steady))) ##get the control genes
must = [r["gene"] for r in rows if r["symbol"] in must_genes] #get the known genes which are needed
genes = set(sig) | set(controls) | set(must) ##create our collection of genes
print(f"Keeping {len(genes)} genes") ##little helpful debug statement could remove at a later point probably

#filter out the genes whose transcripts are not in the gencode file
with open("ref/tx2gene.csv") as f:
    txs = {r["tx"] for r in csv.DictReader(f) if r["gene"] in genes}

fasta_out = out / "ref" / "transcripts.fa.gz" ##this will be the outputted fasta file
keep = False
with gzip.open("ref/gencode.v48.transcripts.fa.gz", "rt") as fin, gzip.open(fasta_out, "wt") as fout: ##open the gencode and fasta files
    for line in fin: ## the actual loop that filters
        if line.startswith(">"):
            keep = line[1:].split("|")[0] in txs
        if keep:
            fout.write(line)
            
# map reads to transcripts and keep them
subprocess.run(["salmon", "index", "-t", str(fasta_out), "-i", str(temp / "idx"), "--gencode"], check=True) ##runs a command that builds a salmon index 
for s in samples:
    sam = temp / f"{s}.sam"
    subprocess.run(["salmon", "quant", "-i", str(temp / "idx"), "-l", "A", ##runs a commnad that does a salmon quant
                    "-1", f"data/{s}_1.fastq.gz", "-2", f"data/{s}_2.fastq.gz",
                    "-p", "4", "-o", str(temp / s), f"--writeMappings={sam}"], check=True) 
    with open(sam) as f:
        names = sorted({line.split("\t", 1)[0] for line in f if not line.startswith("@")}) ##reads the sam file that was created above
    names = random.sample(names, min(max_reads, len(names))) #get a random set
    (temp / f"{s}.txt").write_text("\n".join(names) + "\n") #save the random sample
    for r in ("1", "2"):
        reads = subprocess.run(["seqtk", "subseq", f"data/{s}_{r}.fastq.gz", str(temp /f"{s}.txt")], check=True, capture_output=True).stdout ##another command that will run seqtk and set reads to it
        with gzip.open(out / "data" / f"{s}_{r}.fastq.gz", "wb") as fo: ## save the reads in a fastq file
            fo.write(reads)
    print(f"{s}: kept {len(names)} read pairs")
    

## now just configure for the test run

with open(out / "samples.csv", "w") as f:
    f.write("sample,donor,treatment\n")
    for s, (donor, treatment) in samples.items():
        f.write(f"{s},{donor},{treatment}\n")

(out / "config.yaml").write_text("""samples: samples.csv
reads_per_sample: 1000000
transcripts: ref/transcripts.fa.gz
salmon_index: ref/salmon_index
threads: 2
min_counts: 10
min_samples: 2
""")

shutil.rmtree(temp) ## clean up unimportant files

print("Test data should now be in tests/data")

