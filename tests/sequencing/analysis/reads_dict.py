
from tests.sequencing.analysis.reads_dict_tools import \
    get_fastq_record_triplet, ReadIdGen


ID = ReadIdGen()


ASSEMBLED = "ASSEMBLED"
UNASSEMBLED = "UNASSEMBLED"


def bc1(sr1, sr2, srm=None):
    return get_fastq_record_triplet(
        read_id=str(ID),
        sr1=sr1,
        sr2=sr2,
        srm=srm,
        barcode1="TCCGCGAA",
        barcode2="GTACTGAC",
    )


READS_DICT_ADAM = {
"bc1": {
    ASSEMBLED: {
        28734: {
            frozenset(((3,7),)): [
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
            frozenset(((3,6),)): [
                bc1(
                    sr1="AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC",
                    sr2="GGAAGTAGTGTGTGCTTGGACTTTTCTTCTTCTTCTTCTTCTCTGAGGGGAGAGGGTAGAGCAGCTGAGAGAGAGCGAGAGCGCTTCGGAACACACAGGAGATGTCGCGCGTGCTTCCTCTAAGAGATTTTCTCTTTGGAAGGGGGAGAAGCCTTTAGATCGGAAGAGCGTCGTGGGGGGAAAGGGTGGAGGGGATGGGGTAGATATCGGGGGGCGG",
                    srm="AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC",
                ),bc1(
                    sr1="AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC",
                    sr2="GGAAGTAGTGTGTGCTTGGACTTTTCTTCTTCTTCTTCTTCTCTGAGGGGAGAGGGTAGAGCAGCTGAGAGAGAGCGAGAGCGCTTCGGAACACACAGGAGATGTCGCGCGTGCTTCCTCTAAGAGATTTTCTCTTTGGAATGGGGAGAAGCCTTT",
                    srm="AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC",
                ),
            ],
        },
    },
    UNASSEMBLED: {
        28734: {
            frozenset(((3,7),)): [
                bc1(
                    sr1="AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC",
                    sr2="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                ),
            ],
            frozenset(((3,5),)): [
                bc1(
                    sr1="AAAGGCTTCTCCCCACTCCAAAGAGAAAATCTCTTAGAGGAAGCACGCGCGACATCTCCTGTGTGTTCCGAAGCGCTCTCGCTCTCTCTCAGCTGCTCTACCCTCTCCCCTCAGAGAAGAAGAAGAAGAAAAGTCCAAGCACACACTACTTCC",
                    sr2="AAAAAAAAAAAAAAAAAAAAAAAAGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                ),
            ],
        },
    },
},
}


