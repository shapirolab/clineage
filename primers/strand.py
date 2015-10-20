__author__ = 'ofirr'

class StrandBaseMixin(object):
    @property
    def ref_sequence(self):
        raise NotImplementedError()

    @property
    def sequence(self):
        raise NotImplementedError()

class StrandPlusMixin(StrandBaseMixin):
    @property
    def ref_sequence(self):
        return self.sequence

    @property
    def sequence(self):
        return self.ref_sequence

class StrandMinusMixin(StrandBaseMixin):
    @property
    def ref_sequence(self):
        return self.sequence.rev_comp()

    @property
    def sequence(self):
        return self.ref_sequence.rev_comp()
