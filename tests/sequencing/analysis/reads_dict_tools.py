from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


_DESC_FMT = "{} {}:N:0:{}+{}"
_FWD_FMT = _DESC_FMT.format("{}", "1", "{}", "{}")
_REV_FMT = _DESC_FMT.format("{}", "2", "{}", "{}")


R1 = "R1"
R2 = "R2"
RM = "RM"


def get_fastq_record_triplet(read_id, sr1, sr2, srm, barcode1, barcode2):
    seq_r1 = Seq(sr1)
    seq_r2 = Seq(sr2)
    seq_rm = Seq(srm)
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
        letter_annotations={"phred_quality": [40] * len(seq_r2)},
    )
    recm = SeqRecord(
        seq_rm,
        id=read_id,
        name=read_id,
        description=fwd_desc,
        letter_annotations={"phred_quality": [40] * len(seq_rm)},
    )
    return {R1: rec1, R2: rec2, RM: recm}


_CUSTOM_READ_ID_FMT = "M00321:123:000000000-ABCDE:1:1111:22222:{:05}"


def _read_id_gen(start=1, FMT=_CUSTOM_READ_ID_FMT):
    for i in xrange(start,100000):
        yield FMT.format(i)


class ReadIdGen(object):
    def __init__(self, start=1):
        self._gen = _read_id_gen(start)

    def __str__(self):
        return next(self._gen)


def _iterate_nested_keys(nested_dict):
    if isinstance(nested_dict, dict):
        for k, inner_d in nested_dict.iteritems():
            for inner_k in _iterate_nested_keys(inner_d):
                yield (k,) + inner_k
    else:
        yield tuple()


def _get_nested_key(nested_dict, nested_key):
    if not isinstance(nested_dict, dict):
        if nested_key is tuple():
            return nested_dict
        else:
            raise KeyError("Key is too deeply nested.")
    else:
        return _get_nested_key(nested_dict[nested_key[0]], nested_key[1:])


def _iterate_nested_items(nested_dict):
    if isinstance(nested_dict, dict):
        for k, inner_d in nested_dict.iteritems():
            for inner_k, v in _iterate_nested_items(inner_d):
                yield ((k,) + inner_k), v
    else:
        yield tuple(), nested_dict


def flatten_nested_list_dict(nested_dict):
    flat = {}
    for nk, nv in _iterate_nested_items(nested_dict):
        for i in xrange(1,len(nk)+1):
            a = flat.setdefault(nk[:i],[])
            a += nv
    for k in nested_dict.iterkeys():
        if (k,) in flat:
            flat[k] = flat[k,]
    return flat


def _deinterlace(interlaced):
    d = {}
    for frame in interlaced:
        for k, v in frame.iteritems():
            d.setdefault(k,[]).append(v)
    return d


def _update_deinterlaced_dict(d, e):
    for k, v in e.iteritems():
        a = d.setdefault(k,[])
        a += v


class FlatDict(object):
    
    def __init__(self, nested_dict):
        self._d = {}
        self._keys_tree = {}
        for nk, nv in _iterate_nested_items(nested_dict):
            denv = _deinterlace(nv)
            for i in xrange(1,len(nk)+1):
                a = self._d.setdefault(nk[:i],{})
                _update_deinterlaced_dict(a, denv)
            for i in xrange(0,len(nk)):
                a = self._keys_tree.setdefault(nk[:i],set())
                a.add(nk[i])
        for k in nested_dict.iterkeys():
            if (k,) in self._d:
                self._d[k] = self._d[k,]
            if (k,) in self._keys_tree:
                self._keys_tree[k] = self._keys_tree[k,]

    def keys(self, nk=()):
        return iter(self._keys_tree[nk])

    def items(self, nk=()):
        for k in self.keys(nk):
            yield k, self[nk + (k,)]

    def __getitem__(self, nk):
        return self._d[nk]
