from Bio.Seq import Seq, Alphabet
from Bio.SeqRecord import SeqRecord


_DESC_FMT = "{} {}:N:0:{}+{}"
_FWD_FMT = _DESC_FMT.format("{}", "1", "{}", "{}")
_REV_FMT = _DESC_FMT.format("{}", "2", "{}", "{}")


R1 = "R1"
R2 = "R2"
RM = "RM"
READS = "READS"


def strip_fasta_records(fasta_records):
    for rec in fasta_records:
        yield rec.id, str(rec.seq)


def sr_to_tup(sr):
    return (str(sr.seq),) + tuple(getattr(sr, k) for k in [
        "id",
        "name",
        "description",
    ]) + (tuple(sr.letter_annotations["phred_quality"]),)


def rc_sr_to_tup(sr):
    return (str(sr.seq.reverse_complement()),) + tuple(getattr(sr, k) for k in [
        "id",
        "name",
        "description",
    ]) + (tuple(sr.letter_annotations["phred_quality"]),)


def srs_to_tups(srs):
    for sr in srs:
        yield sr_to_tup(sr)


def rc_srs_to_tups(srs):
    for sr in srs:
        yield rc_sr_to_tup(sr)


def get_fastq_record_triplet(read_id, sr1, sr2, srm, barcode1, barcode2):
    seq_r1 = Seq(sr1, Alphabet.SingleLetterAlphabet())
    seq_r2 = Seq(sr2, Alphabet.SingleLetterAlphabet())
    if srm is not None:
        seq_rm = Seq(srm, Alphabet.SingleLetterAlphabet())
    fwd_desc = _FWD_FMT.format(read_id, barcode1, barcode2)
    rev_desc = _REV_FMT.format(read_id, barcode1, barcode2)
    rec1 = SeqRecord(
        seq_r1,
        id=read_id,
        name=read_id,
        description=fwd_desc,
        letter_annotations={"phred_quality": [40] * len(seq_r1)},
    )
    rec2 = SeqRecord(
        seq_r2,
        id=read_id,
        name=read_id,
        description=rev_desc,
        letter_annotations={"phred_quality": [30] * len(seq_r2)},
    )
    if srm is not None:
        recm = SeqRecord(
            seq_rm,
            id=read_id,
            name=read_id,
            description=fwd_desc,
            letter_annotations={"phred_quality": [40] * len(seq_rm)},
        )
    if srm is not None:
        return {R1: rec1, R2: rec2, RM: recm}
    else:
        return {R1: rec1, R2: rec2}


_CUSTOM_READ_ID_FMT = "M00321:123:000000000-ABCDE:1:1111:22222:{:05}"


def _read_id_gen(start=1, FMT=_CUSTOM_READ_ID_FMT):
    for i in range(start,100000):
        yield FMT.format(i)


class ReadIdGen(object):
    def __init__(self, start=1):
        self._gen = _read_id_gen(start)

    def __str__(self):
        return next(self._gen)
