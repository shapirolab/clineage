from sequencing.calling.hist import Histogram


class MonoSimulatedHistogram(Histogram):
    def __init__(self, ms_len, simulation_cycle, simulated_hist):
        assert isinstance(ms_len, int)
        assert isinstance(simulation_cycle, int)
        self._ms_len = ms_len
        self._sim_cyc = simulation_cycle
        super(MonoSimulatedHistogram, self).__init__(simulated_hist)

    def __repr__(self):
        return '{}\n{}'.format(self.allele_frozenset, super().__repr__())

    @property
    def simulation_cycle(self):
        return self._sim_cyc

    @property
    def ms_len(self):
        return self._ms_len

    @property
    def allele_frozenset(self):
        return frozenset([self.ms_len])


class MultiSimulatedHistogram(Histogram):
    def __init__(self, ms_lens, simulation_cycle, simulated_hist):
        assert isinstance(ms_lens, frozenset)
        assert isinstance(simulation_cycle, int)
        self._ms_lens = ms_lens
        self._sim_cyc = simulation_cycle
        super(MultiSimulatedHistogram, self).__init__(simulated_hist)

    def __repr__(self):
        return '{}\n{}'.format(self.allele_frozenset, super().__repr__())

    @property
    def simulation_cycle(self):
        return self._sim_cyc

    @property
    def ms_lens(self):
        return self._ms_lens

    @property
    def allele_frozenset(self):
        return self.ms_lens


class ProportionalMultiSimulatedHistogram(Histogram):
    def __init__(self, ms_lens_and_proportions, simulation_cycle, simulated_hist):
        assert isinstance(ms_lens_and_proportions, frozenset)
        assert isinstance(simulation_cycle, int)
        self._ms_lens_and_proportions = ms_lens_and_proportions
        self._sim_cyc = simulation_cycle
        self._alleles_to_proportions = {a: p for a, p in ms_lens_and_proportions if p > 0}
        assert sum(self._alleles_to_proportions.values()) == 1
        super(ProportionalMultiSimulatedHistogram, self).__init__(simulated_hist)

    def __repr__(self):
        return '{}\n{}'.format(self.alleles_to_proportions, super().__repr__())

    @property
    def simulation_cycle(self):
        return self._sim_cyc

    @property
    def alleles_to_proportions(self):
        return self._alleles_to_proportions

    @property
    def ms_lens_and_proportions(self):
        return self._alleles_to_proportions.items()

    @property
    def allele_frozenset(self):
        return frozenset(self._alleles_to_proportions.keys())
