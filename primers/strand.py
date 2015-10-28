__author__ = 'ofirr'

class BaseStrandMixin(object):
    @property
    def ref_sequence(self):
        raise NotImplementedError()

    @property
    def sequence(self):
        raise NotImplementedError()

class PlusStrandMixin(BaseStrandMixin):
    @property
    def ref_sequence(self):
        return self.sequence

    @property
    def sequence(self):
        return self.ref_sequence

class MinusStrandMixin(BaseStrandMixin):
    @property
    def ref_sequence(self):
        return self.sequence.rev_comp()

    @property
    def sequence(self):
        return self.ref_sequence.rev_comp()
