from sequencing.calling.models import CalledAlleles


def get_mono_mutations_dict(srs, calling_scheme, confidence_threshold=0.01):
    d = dict()
    for sr in srs:
        d[sr] = dict()
        for ca in CalledAlleles.objects.filter(
                calling_scheme=calling_scheme,
                histogram__sample_reads=sr).select_subclasses():
            if ca.confidence > confidence_threshold:
                continue
            d[sr][ca.microsatellite] = ca.genotypes.allele1
    return d

