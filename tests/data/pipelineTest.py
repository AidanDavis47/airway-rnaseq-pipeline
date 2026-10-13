##these are the actual tests that are being done

##imports
import csv, shutil, subprocess
from pathlib import Path
import pytest ##needed to run actual tests

REPO = Path(__file__).resolve().parents[1] #find the project folder 

@pytest.fixture(scope="session")

##this will get the test data and copy it to a temporary folder and then run the pipeline there
def run_pipeline(temp_path):
    d = temp_path.mktemp("pipeline")
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
def test_SamplesQuant(run_dir):
    for s in get_samples(run_dir):
        q = run_dir / "results/salmon" / s / "quant.sf"
        assert q.exists() and q.stat().st_size > 0, f"missing {q}"
