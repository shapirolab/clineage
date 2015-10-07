from django.test import TestCase
from misc.models import Taxa
from genomes.models import Assembly, Chromosome
import os


class LinappModelTestCase(TestCase):
    seq = "TTAAGTAACATCAGCCAAGCA"
    start = 1700000
    stop = start + len(seq) - 1


    @classmethod
    def setUpClass(cls):
        t = Taxa.objects.create(
            name="Homo sapiens",
            taxonomy_id=7,
            rank="Specie",
            friendly_name="Human",
        )

        a = Assembly.objects.create(
            taxa=t,
            name="Human #19",
            friendly_name="hg19",
        )

        cls.cx = Chromosome.objects.create(
            name="X",
            assembly=a,
            sequence_length=1
        )

        cls.cx.sequence_length = os.path.getsize(cls.cx.get_abs_path())
        cls.cx.save()

    def test_getdna(self):
        self.assertEquals(self.seq, self.cx.getdna(self.start, self.stop))

    def test_locate(self):
        self.assertEquals((self.start, self.stop), self.cx.locate(self.start, self.stop, self.seq))
        self.assertEquals((self.start, self.stop), self.cx.locate(self.start-5, self.stop-5, self.seq))
        self.assertEquals((self.start, self.stop), self.cx.locate(self.start+5, self.stop+5, self.seq))
        self.assertRaises(ValueError, self.cx.locate, self.start+50, self.stop+50, self.seq)

