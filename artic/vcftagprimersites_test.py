# test_vcftagprimersites.py is the unit tests for the vcf primer site tagging
import pytest

from artic import vcftagprimersites


def test_read_bed_file():

    # process the nCoV-2019 primer scheme
    bedfile = vcftagprimersites.read_bed_file(
        "../test-data/nCoV-2019.scheme.bed")

    # check the the alts have been collapsed into a canonical primer site
    assert len(bedfile) == 196, "alts where not collapsed"

    # TODO: this is a starting point for the unit tests - more to come...
