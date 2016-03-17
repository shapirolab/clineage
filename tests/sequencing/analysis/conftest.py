import pytest
import os
import uuid
from django.conf import settings

from sequencing.analysis.models import AdamMSVariations, AdamHistogram
from sequencing.analysis.models import LEFT, RIGHT
from misc.utils import get_unique_path

from tests.sequencing.runs.conftest import *
from tests.sequencing.analysis.pu_28727_adam_ms_variations import VARS_28727
from tests.lib_prep.workflows.conftest import *
from tests.targeted_enrichment.amplicons.conftest import *
from tests.sequencing.analysis.fixture_generators import \
    generate_samplereads, generate_mergedreads, generate_readsindex, \
    generate_adammarginassignment, generate_adamampliconreads, generate_amsv, \
    generate_adamhistogram
from tests.utils import to_fixture
from tests.sequencing.analysis.reads_dict_tools import flatten_nested_list_dict
from tests.sequencing.analysis.reads_dict import READS_DICT_ADAM


@pytest.fixture(scope="session")
def reads_dict_adam():
    return flatten_nested_list_dict(READS_DICT_ADAM)


samplereads_bc1 = to_fixture(generate_samplereads, "bc1")
mergedreads_bc1 = to_fixture(generate_mergedreads, "bc1")
readsindex_bc1_M = to_fixture(generate_readsindex, "bc1", "M")
readsindex_bc1_F = to_fixture(generate_readsindex, "bc1", "F")
adammarginassignment_bc1_F = to_fixture(generate_adammarginassignment, "bc1", "F")
adamampliconreads_bc1_F_28734 = to_fixture(generate_adamampliconreads, "bc1", "F", "28734")
amsv_28727 = to_fixture(generate_amsv, "28727")
amsv_28734 = to_fixture(generate_amsv, "28734")
adamhistogram = to_fixture(generate_adamhistogram, "bc1", "F", "28734")

samplereads_bc2 = to_fixture(generate_samplereads, "bc2")
readsindex_bc2_M = to_fixture(generate_readsindex, "bc2", "M")
adammarginassignment_bc2_M = to_fixture(generate_adammarginassignment, "bc2", "M")
adamampliconreads_bc2_M_28734 = to_fixture(generate_adamampliconreads, "bc2", "M", "28734")
adamhistogram_plus3 = to_fixture(generate_adamhistogram, "bc2", "M", "28734")

@pytest.fixture()
def merged_reads_stripped_fasta():
    return {
        (
            'M00393:167:000000000-ABF3N:1:1101:20583:16769',
            'AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1101:7043:8470',
            'AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1102:6086:12380',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1106:14231:8476',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1109:4497:5075',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1112:13016:19696',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:1112:20425:16124',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:2101:12593:17954',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:2104:19595:12515',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:2109:21650:4662',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:2115:26810:9430',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ), (
            'M00393:167:000000000-ABF3N:1:2117:19711:1998',
            'AAAGGCTTCTCCCCATTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC',
        ),
    }


@pytest.fixture()
def reads_matches(pu_28734):
    return {
        'M00393:167:000000000-ABF3N:1:1101:20583:16769':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:1101:7043:8470':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:1102:6086:12380':
            {(pu_28734,  LEFT), (pu_28734,  RIGHT)},
        'M00393:167:000000000-ABF3N:1:1106:14231:8476':
            {(pu_28734,  LEFT), (pu_28734,  RIGHT)},
        'M00393:167:000000000-ABF3N:1:1109:4497:5075':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:1112:13016:19696':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:1112:20425:16124':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:2101:12593:17954':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:2104:19595:12515':
            {(pu_28734, LEFT), (pu_28734, RIGHT)},
        'M00393:167:000000000-ABF3N:1:2109:21650:4662':
            {(pu_28734,  LEFT), (pu_28734,  RIGHT)},
        'M00393:167:000000000-ABF3N:1:2115:26810:9430':
            {(pu_28734,  LEFT), (pu_28734,  RIGHT)},
        'M00393:167:000000000-ABF3N:1:2117:19711:1998':
            {(pu_28734,  LEFT), (pu_28734,  RIGHT)}
    }


@pytest.fixture()
def reads_amplicons(pu_28734):
    return {
        ('M00393:167:000000000-ABF3N:1:1101:20583:16769',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1101:7043:8470',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1102:6086:12380',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1106:14231:8476',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1109:4497:5075',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1112:13016:19696',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:1112:20425:16124',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:2101:12593:17954',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:2104:19595:12515',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:2109:21650:4662',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:2115:26810:9430',
            pu_28734),
        ('M00393:167:000000000-ABF3N:1:2117:19711:1998',
            pu_28734),
    }


@pytest.fixture()
def reads_by_amplicons(pu_28734):
    return {
        pu_28734: {
            'M00393:167:000000000-ABF3N:1:1101:20583:16769',
            'M00393:167:000000000-ABF3N:1:1101:7043:8470',
            'M00393:167:000000000-ABF3N:1:1102:6086:12380',
            'M00393:167:000000000-ABF3N:1:1106:14231:8476',
            'M00393:167:000000000-ABF3N:1:1109:4497:5075',
            'M00393:167:000000000-ABF3N:1:1112:13016:19696',
            'M00393:167:000000000-ABF3N:1:1112:20425:16124',
            'M00393:167:000000000-ABF3N:1:2101:12593:17954',
            'M00393:167:000000000-ABF3N:1:2104:19595:12515',
            'M00393:167:000000000-ABF3N:1:2109:21650:4662',
            'M00393:167:000000000-ABF3N:1:2115:26810:9430',
            'M00393:167:000000000-ABF3N:1:2117:19711:1998',
        }
    }


@pytest.fixture()
def pu_28727_adam_ms_variations():
    return VARS_28727


@pytest.fixture()
def require_adammsvariations(amsv_28727, amsv_28734):
    pass
