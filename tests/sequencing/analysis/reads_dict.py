
from tests.sequencing.analysis.reads_dict_tools import \
    get_fastq_record_triplet, ReadIdGen


ID = ReadIdGen()


def bc1(sr1, sr2, srm):
    return get_fastq_record_triplet(
        read_id=str(ID),
        sr1=sr1,
        sr2=sr2,
        srm=srm,
        barcode1="ACACAC",  # FIXME: match the barcode to BC1.
        barcode2="AGAGAG",  # FIXME: match the barcode to BC1.
    )


READS_DICT1 = {
"bc1": {
    "M": {
        "28734": {
            ((3,7),): [
                bc1(
                    sr1="AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC",
                    sr2="GGAAGTAGTGTGTGCTTGGACTTTTCTTCTTCTTCTTCTTCTTCTCTGAGGGGAGAGGGTAGAGCAGCTGAGAGAGAGCGAGAGCGCTTCGGAACACACAGGAGATGTCGCGCGTGCTTCCTCTAAGAGATTTTCTCTTTGGAAGGGGGAGAAGCCTTTAGATCGGAAGAGCGTCGTGGGGGGAAAGGGTGGAGGGGATGGGGTAGATATCGGGGGGCGG",
                    srm="AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC",
                ),bc1(
                    sr1="AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC",
                    sr2="GGAAGTAGTGTGTGCTTGGACTTTTCTTCTTCTTCTTCTTCTTCTCTGAGGGGAGAGGGTAGAGCAGCTGAGAGAGAGCGAGAGCGCTTCGGAACACACAGGAGATGTCGCGCGTGCTTCCTCTAAGAGATTTTCTCTTTGGAATGGGGAGAAGCCTTT",
                    srm="AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC",
                ),
            ],
        },
    },
    "F": {
        "28734": {
            ((3,7),): [
            ],
        },
    },
},
}


