from sequencing.calling.hist import Histogram


class MonoSimulatedHistogram(Histogram):
    def __init__(self, ms_len, simulation_cycle, simulated_hist):
        assert isinstance(ms_len, int)
        assert isinstance(simulation_cycle, int)
        self._ms_len = ms_len
        self._sim_cyc = simulation_cycle
        super(MonoSimulatedHistogram, self).__init__(simulated_hist)

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
        super(ProportionalMultiSimulatedHistogram, self).__init__(simulated_hist)

    @property
    def simulation_cycle(self):
        return self._sim_cyc

    @property
    def ms_lens_and_proportions(self):
        return self._ms_lens_and_proportions

    @property
    def allele_frozenset(self):
        return frozenset([a for a, p in self.ms_lens_and_proportions])
