
def _iterate_nested_items(nested_dict):
    if isinstance(nested_dict, dict):
        for k, inner_d in nested_dict.iteritems():
            for inner_k, v in _iterate_nested_items(inner_d):
                yield ((k,) + inner_k), v
    else:
        yield tuple(), nested_dict


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

    def __init__(self, nested_dict, deinterlaced_keys=None):
        self._d = {}
        self._keys_tree = {}
        if deinterlaced_keys is not None:
            self._deinterlaced_keys = deinterlaced_keys
            self._has_default = True
        else:
            self._deinterlaced_keys = []
            self._has_default = False
        for nk, nv in _iterate_nested_items(nested_dict):
            denv = _deinterlace(nv)
            for i in xrange(1,len(nk)+1):
                a = self._d.setdefault(nk[:i],
                    {k: [] for k in self._deinterlaced_keys})
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
        if self._has_default:
            return self._d.get(nk, {k: [] for k in self._deinterlaced_keys})
        else:
            return self._d[nk]