"""minion_validator.py runs the functional tests of the minion pipeline, for both nanopolish and medaka workflows.

For each test dataset, it will make the following checks for each workflow:

    * check for existence of the final consensus sequence
    * check the final consensus matches the expected consensus
    * check for existence of the final VCF file (generated by artic_vcf_filter)
    * check all expected variants are reported
    * check no unexpected variants are reported
    * if there is an expected deletion, the consensus sequence will be checked for it's absence

NOTE:

    * only a basic grep is used for checking the deletion in the consensus - test will fail if deletion sequence is present multiple times (I should improve this part of the test...)

"""
from Bio import SeqIO
from tqdm.auto import tqdm
import argparse
import errno
import glob
import os
import pytest
import requests
import sys
import tarfile
import vcf


from . import pipeline

# help pytest resolve where test data is kept
TEST_DIR = os.path.dirname(os.path.abspath(__file__))

# dataDir is where to store the test data locally
dataDir = TEST_DIR + "/../test-data/"

# testData is a lookup of sampleIDs to download urls
testData = {
    "MT007544": "https://raw.githubusercontent.com/artic-network/fieldbioinformatics/master/test-data/MT007544/MT007544.fastq",
    "CVR1": "https://artic.s3.climb.ac.uk/validation-sets/CVR1.tgz",
    "NRW01": "http://artic.s3.climb.ac.uk/validation-teams/NRW01.tgz",
    "SP1": "http://artic.s3.climb.ac.uk/validation-teams/SP1.tgz"
}

# refConsensuses is a nested dict of sample IDs and their expected consensus sequences from the nanopolish workflow
refConsensuses = {
    "CVR1": dataDir + "consensus-sequences/CVR1.consensus.fasta",
    "NRW01": dataDir + "consensus-sequences/NRW01.consensus.fasta",
    "SP1": dataDir + "consensus-sequences/SP1.consensus.fasta"
}

# refMedakaConsensuses is a nested dict of sample IDs and their expected consensus sequences from the medaka workflow
refMedakaConsensuses = {
    "MT007544": dataDir + "consensus-sequences/MT007544.consensus.medaka.fasta",
    "CVR1": dataDir + "consensus-sequences/CVR1.consensus.medaka.fasta",
    "NRW01": dataDir + "consensus-sequences/NRW01.consensus.medaka.fasta",
    "SP1": dataDir + "consensus-sequences/SP1.consensus.medaka.fasta"
}

# nanopolishTestVariants is a nested dict of sample IDs and their expected variants when using the nanopolish workflow
nanopolishTestVariants = {
    "CVR1":
        {
            # pos: (ref, alt, type, count)
            241: ['C', 'T', "snp", 1],
            3037: ['C', 'T', "snp", 1],
            12733: ['C', 'T', "snp", 2], # count 2 indicates this is reported twice in the VCF (due to pool groups)
            14408: ['C', 'T', "snp", 1],
            23403: ['A', 'G', "snp", 1],
            27752: ['C', 'T', "snp", 1],
            28881: ['G', 'A', "snp", 1],
            28882: ['G', 'A', "snp", 1],
            28883: ['G', 'C', "snp", 1]
        },
    "NRW01":
        {
            1440: ['G', 'A', "snp", 1],
            2891: ['G', 'A', "snp", 1],
            4655: ['C', 'T', "snp", 1],
            8422: ['G', 'A', "snp", 1],
            22323: ['C', 'T', "snp", 2],
            29546: ['C', 'A', "snp", 2]
        },
    "SP1":
        {
            241: ['C', 'T', "snp", 1],
            3037: ['C', 'T', "snp", 1],
            14408: ['C', 'T', "snp", 1],
            23403: ['A', 'G', "snp", 1]
        },
}

# medakaTestVariants is a nested dict of sample IDs and their expected variants when using the medaka workflow
medakaTestVariants = {
    "MT007544":
        {
            # pos: (ref, alt, type, count)
            29749: ["ACGATCGAGTG", 'A', "del", 1]
        },
    "CVR1":
        {
            241: ['C', 'T', "snp", 1],
            3037: ['C', 'T', "snp", 1],
            12733: ['C', 'T', "snp", 1],
            12733: ['C', 'T', "snp", 2],
            14408: ['C', 'T', "snp", 1],
            23403: ['A', 'G', "snp", 1],
            27752: ['C', 'T', "snp", 1],
            28881: ['GGG', 'AAC', "indel", 1],
        },
    "NRW01":
        {
            1440: ['G', 'A', "snp", 1],
            2891: ['G', 'A', "snp", 1],
            4655: ['C', 'T', "snp", 1],
            8422: ['G', 'A', "snp", 1],
            22323: ['C', 'T', "snp", 2],
            29546: ['C', 'A', "snp", 2]
        },
    "SP1":
        {
            241: ['C', 'T', "snp", 1],
            3037: ['C', 'T', "snp", 1],
            14408: ['C', 'T', "snp", 1],
            23403: ['A', 'G', "snp", 1]
        },
}

# dataChecker will run before the tests and download all the test datasets if not present
@pytest.fixture(scope="session", autouse=True)
def dataChecker(numValidations):
    print("checking for validation datasets...")
    for sampleID, url in testData.items():
        targetPath = dataDir + sampleID
        if os.path.exists(targetPath) == False:
            print("\tno data for {}" .format(sampleID))
            print("\tmaking dir at {}" .format(targetPath))
            os.mkdir(targetPath)
            print("\tdownloading from {}" .format(url))
            try:
                download(url, dataDir, sampleID)
            except Exception as e:
                print("download failed: ", e)
                sys.exit(1)
        else:
            print("\tfound data dir for {}" .format(sampleID))
    print("validation datasets ready\n")

# download will download and untar a test dataset
def download(url, dataDir, sampleID):
    filename = url.rsplit('/', 1)[1] 
    with open(f'{dataDir}/{filename}', 'wb+') as f:
        response = requests.get(url, stream=True)
        total = int(response.headers.get('content-length'))
        if total is None:
            f.write(response.content)
        else:
            with tqdm(total=total, unit='B', unit_scale=True, desc=filename) as pbar:
                for data in tqdm(response.iter_content(chunk_size=1024)):
                    f.write(data)
                    pbar.update(1024)

    tar = tarfile.open(dataDir + "/" + filename, "r:gz")
    tar.extractall(dataDir)
    tar.close()
    os.remove(dataDir + "/" + filename)

# genCommand will create the minion command for the requested workflow (nanopolish or medaka)
def genCommand(sampleID, workflow):
    if workflow=="nanopolish":
        return [
            "minion",
            "--read-file",
            TEST_DIR + "/../test-data/" + sampleID + "/*.fast[aq]",
            "--scheme-directory",
            TEST_DIR + "/../test-data/primer-schemes",
            "--fast5-directory",
            TEST_DIR + "/../test-data/" + sampleID + "/fast5",
            "--sequencing-summary",
            TEST_DIR + "/../test-data/" + sampleID + "/" + sampleID + "_sequencing_summary.txt",
            "nCoV-2019/V1",
            sampleID
        ]
    if workflow=="medaka":    
        return [
            "minion",
            "--medaka",
            "--read-file",
            TEST_DIR + "/../test-data/" + sampleID + "/*.fast[aq]",
            "--scheme-directory",
            TEST_DIR + "/../test-data/primer-schemes",
            "nCoV-2019/V1",
            sampleID
    ]

# cleanUp will remove all files generated by the pipeline for a given sampleID
def cleanUp(sampleID):
    fileList = glob.glob("{}*" .format(sampleID))
    for filePath in fileList:
        try:
            os.remove(filePath)
        except:
            print("Error while deleting file : ", filePath)

# checkConsensus will return 1 if a subsequence is present in a consensus file, or 0 if absent
def checkConsensus(consensusFile, subSeq):
    for record in SeqIO.parse(open(consensusFile, 'r'), 'fasta'):
        if subSeq in record.seq:
            return 1
    return 0

# runner is the test runner
def runner(workflow, numValidations):

    # get the workflow data
    if workflow == "nanopolish":
        data = nanopolishTestVariants
    elif workflow == "medaka":
        data = medakaTestVariants
    else:
        print("error setting up test runner")
        assert False

    # check the number of validation datasets requested
    counter = numValidations
    if (numValidations < 0) or (numValidations > len(testData)):
        counter = len(testData)

    # run the data through the requested workflow, then validate the output
    print("running the {} workflow with {} of the validation datasets...\n" .format(workflow, counter))
    for sampleID, expVariants in data.items():

        # break when the requested number of datasets have been run
        if counter == 0:
            break

        # generate the command
        cmd = genCommand(sampleID, workflow)

        # setup and run the minion parser
        parser = pipeline.init_pipeline_parser()

        # parse the arguments
        try:
            args = parser.parse_args(cmd)
        except SystemExit:
            print("failed to parse valid command for `artic minion`")
            assert False

        # run the minion pipeline
        try:
            args.func(parser, args)
        except SystemExit:
            print("artic minion exited early with an error")
            assert False

        # check the ARTIC consensus was created
        consensusFile = "%s.consensus.fasta" % sampleID
        assert os.path.exists(consensusFile) == True, "no consensus produced for {}" .format(sampleID)
        testSeqs = SeqIO.parse(open(consensusFile, 'r'), 'fasta')
        testConsensus = next(testSeqs)

        # check the ARTIC consensus sequence matches the one on record
        if workflow == "medaka":
            refSeqs = SeqIO.parse(open(refMedakaConsensuses[sampleID], 'r'), 'fasta')
        else:
            refSeqs = SeqIO.parse(open(refConsensuses[sampleID], 'r'), 'fasta')
        refConsensus = next(refSeqs)
        assert testConsensus.seq == refConsensus.seq, "produced ARTIC consensus does not match expected consensus for {}" .format(sampleID)

        # check the ARTIC VCF was created
        vcfFile = "%s.pass.vcf.gz" % sampleID
        assert os.path.exists(
            vcfFile) == True, "no VCF produced for {}" .format(sampleID)

        # open the VCF and check the reported variants match the expected
        for record in vcf.Reader(filename=vcfFile):
            if record.POS in expVariants:
                assert record.REF == expVariants[record.POS][0], "incorrect REF reported in VCF for {} at position {}" .format(sampleID, record.POS)
                assert str(record.ALT[0]) == expVariants[record.POS][1], "incorrect ALT reported in VCF for {} at position {}" .format(sampleID, record.POS)

                # if this is an expected deletion, check the consensus sequence for it's absence
                if expVariants[record.POS][2] == "del":
                    assert checkConsensus(consensusFile, record.REF) == 0, "expected deletion for {} was reported but was left in consensus" .format(sampleID)

                    # also check that the VCF record is correctly labelled as DEL
                    assert record.is_deletion, "deletion for {} not formatted correctly in VCF" .format(sampleID)

                # if this is an expected indel, check that the VCF record is correctly labelled as INDEL
                if expVariants[record.POS][2] == "indel":
                    assert record.is_indel, "indel for {} not formatted correctly in VCF" .format(sampleID)

                # else, check that the VCF record is correctly labelled as SNP
                if expVariants[record.POS][2] == "snp":
                    assert record.is_snp, "snp for {} not formatted correctly in VCF" .format(sampleID)

                # decrement/remove the variant from the expected list, so we can keep track of checked variants
                expVariants[record.POS][3] -= 1
                if (expVariants[record.POS][3] == 0):
                    del expVariants[record.POS]
            else:
                print("unexpected variant found for {}: {} at {}" .format(sampleID, str(record.ALT[0]), record.POS))
                assert False  

        # check we've confirmed all the expected variants         
        if len(expVariants) != 0:
            print("variants missed during test for {}" .format(sampleID))
            for key, val in expVariants.items():
                print("\t{}: {}" .format(key, val))
            assert False

        # clean up pipeline files
        cleanUp(sampleID)
        counter -= 1

# test_NanopolishMinion is the unit test runner to test the minion pipeline with the nanopolish workflow
@pytest.mark.env("nanopolish")
def test_NanopolishMinion(numValidations):
    runner("nanopolish", numValidations)

# test_MedakaMinion is the unit test runner to test the minion pipeline with the medaka workflow
@pytest.mark.env("medaka")
def test_MedakaMinion(numValidations):
    runner("medaka", numValidations)
