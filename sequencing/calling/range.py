import itertools
from sequencing.calling.proportion_generators import proportions_generator, filtered_proportions_generator


class BaseAlleleNumber(object):
    @property
    def allele_number(self):
        raise NotImplemented


class MSLengthBoundsMixin(object):

    @property
    def ms_len_bounds(self):
        raise NotImplemented


class AllelesRangeMixin(MSLengthBoundsMixin):

    @property
    def alleles(self):
        """
        Alleles space generator
        """
        yield from range(*self.ms_len_bounds)


class MultiAlleleMixin(AllelesRangeMixin, BaseAlleleNumber):

    @property
    def alleles(self):
        for comb in itertools.chain(
            *[itertools.combinations(
                range(*self.ms_len_bounds), allele_number
            ) for allele_number in range(1, self.allele_number+1)]
        ):
            yield frozenset(comb)


class MultiAlleleStrictNMixin(AllelesRangeMixin, BaseAlleleNumber):

    @property
    def alleles(self):
        for comb in itertools.combinations(
                range(*self.ms_len_bounds),
                self.allele_number
        ):
            yield frozenset(comb)


# class ProportionsStepMixin(object):
#
#     @property
#     def proportion_step(self):
#         raise NotImplemented


# class ProportionsRangeMixin(BaseAlleleNumber, ProportionsStepMixin):
class ProportionsRangeMixin(BaseAlleleNumber):

    @property
    def proportions(self):
        yield from proportions_generator(number_of_items=self.allele_number, step=self.proportion_step)


class ProportionsBoundsMixin(object):

    @property
    def proportion_bounds(self):
        raise NotImplemented


class BoundProportionsRangeMixin(ProportionsRangeMixin, ProportionsBoundsMixin):

    @property
    def proportions(self):
        yield from filtered_proportions_generator(
            number_of_items=self.allele_number,
            step=self.proportion_step,
            proportion_bounds=self.proportion_bounds
        )


class ProportionalAllelesMixin(MultiAlleleStrictNMixin):

    @property
    def alleles(self):
        """
        Alleles and proportions space generator
        """

        for alleles, proportions in itertools.product(
            MultiAlleleStrictNMixin.alleles.fget(self),
            self.proportions,
        ):
            yield frozenset(zip(alleles, proportions))
