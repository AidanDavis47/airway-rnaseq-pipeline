##these are the actual tests that are being done

##imports
import csv, shutil, subprocess
from pathlib import Path
import pytest ##needed to run actual tests

REPO = Path(__file__).resolve().parents[1] #find the project folder 

@pytest.fixture(scope="session")

##this will get the test data and copy it to a temporary folder and then run the pipeline there
def run_pipeline(tmp_path_factory):
    d = tmp_path_factory.mktemp("pipeline")
    shutil.copytree(REPO / "tests" / "data", d, dirs_exist_ok=True)
    subprocess.run(["snakemake", "-s", str(REPO/"Snakefile"), "--configfile", "config.yaml", "--cores", "2"], cwd=d, check=True) ## this is what actually tells the pipeline to run
    return d

##simple helper functions

##function to get a list of the ID's of samples that the tests will use
def get_samples(d):
    with open(d / "samples.csv") as f:
        return [r["sample"] for r in csv.DictReader(f)]
    
##function that will get the results that the tests will use
def get_results(d):
    with open(d / "results/deseq2/deseq2_results.csv") as f:
        return list(csv.DictReader(f))
    
## TESTS

##Test 1: Samples are quantified
## this is a test to make sure that the samples are actually getting quantified, we are checking to make sure the quant files exists and contain information
def test_SamplesQuant(run_pipeline):
    for s in get_samples(run_pipeline):
        q = run_pipeline / "results/salmon" / s / "quant.sf"
        assert q.exists() and q.stat().st_size > 0, f"missing {q}"
        
##Test 2: Make sure that the multiqc report is made
def test_Report(run_pipeline):
    assert (run_pipeline / "results/multiqc/multiqc_report.html").exists()
    
##Test 3: Data coloums are correct
##runs get_results and check to see if the coloumns that are expected exists, expected are (baseMean, log2FoldChange, lfcSE, stat, pvalue, pdaj, gene, symbol)
def test_Columns(run_pipeline): ##REMINDER THIS IS THE CORRECT WAY TO SPELL COLUMNS MAY HAVE ERRORS APPEAR IF I SPELT IT WRONG IN OTHER AREAS
    expected = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene", "symbol"}
    assert expected <= set(get_results(run_pipeline)[0].keys())

##Test 4: Check to see if FKBP5 goes up with dexamethasone, basically checking to see if the pipeline is actually getting the right information
## this will also use the get_results function since we are checking the results basically
def test_FKBP5Increase(run_pipeline):
    row = next(r for r in get_results(run_pipeline) if r["symbol"] == "FKBP5")
    assert float(row["log2FoldChange"]) > 1 ##checking to see if the change has increased, 1 would mean there was not increase i think, would that mean less than 1 would mean a decrease?

    

