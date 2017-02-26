import itertools
import decimal


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


def proportions_generator(number_of_items, step):
    """
    Generate the space of proportions that sum to 1 for k items.
    >>> list(proportions_generator(items=2, step=decimal.Decimal('0.2')))
    [(Decimal('0.0'), Decimal('1.0')),
     (Decimal('0.2'), Decimal('0.8')),
     (Decimal('0.4'), Decimal('0.6')),
     (Decimal('0.6'), Decimal('0.4')),
     (Decimal('0.8'), Decimal('0.2')),
     (Decimal('1.0'), Decimal('0.0'))]
    >>> list(proportions_generator(items=3, step=decimal.Decimal('0.5')))
    [(Decimal('0.0'), Decimal('0.0'), Decimal('1.0')),
     (Decimal('0.0'), Decimal('0.5'), Decimal('0.5')),
     (Decimal('0.0'), Decimal('1.0'), Decimal('0.0')),
     (Decimal('0.5'), Decimal('0.0'), Decimal('0.5')),
     (Decimal('0.5'), Decimal('0.5'), Decimal('0.0')),
     (Decimal('1.0'), Decimal('0.0'), Decimal('0.0'))]
    """
    assert isinstance(step, decimal.Decimal)
    number_of_steps = int(1 / step)

    perms = filter(
        lambda x: sum(x) <= number_of_steps,
        itertools.product(
            range(0, number_of_steps+1),
            repeat=number_of_items - 1  # We use items-1 so that we can add the last item as 1-sum(rest of the items)
        )
    )
    for i in perms:
        full_tuple = i + (number_of_steps - sum(i),)  # add the 1-sum item to the tuple
        yield tuple([decimal.Decimal(i) * step for i in full_tuple])  # yield the proportions tuple


def filtered_proportions_generator(number_of_items, step, proportion_bounds):
    """
    Generate the space of proportions that sum to 1 for k items.
    >>> list(proportions_generator(items=2, step=decimal.Decimal('0.2'), proportion_bounds=(decimal.Decimal('0.4'),decimal.Decimal('0.6'))))
    [(Decimal('0.0'), Decimal('1.0')),
     (Decimal('0.4'), Decimal('0.6')),
     (Decimal('0.6'), Decimal('0.4')),
     (Decimal('1.0'), Decimal('0.0'))]
    """
    lower_bound, upper_bound = proportion_bounds
    assert isinstance(lower_bound, decimal.Decimal)
    assert isinstance(upper_bound, decimal.Decimal)
    yield from filter(
        lambda pt: not sum(1 for p in pt if p != decimal.Decimal(0) and p!=decimal.Decimal(1) and (p<lower_bound or p>upper_bound)),
        proportions_generator(
            number_of_items=number_of_items,
            step=step
        ))


class ProportionsRangeMixin(BaseAlleleNumber):

    @property
    def proportion_step(self):
        raise NotImplemented

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


class ProportionalAllelesMixin(MultiAlleleMixin, ProportionsRangeMixin):
    @property
    def alleles(self):
        """
        Alleles and proportions space generator
        """

        yield from itertools.product(
            self.proportions,
            MultiAlleleMixin.alleles.fget(self)
        )


class BoundProportionalAllelesMixin(ProportionalAllelesMixin):
    pass

