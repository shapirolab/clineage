import os
from Bio import SeqIO

os.chdir('/net/mraid11/export/data/dcsoft/home/LINEAGE/GenomesData/Mouse/mm10/')
mm10 = Assembly.objects.all()[2]
for chr in chr_list:
    f = SeqIO.read(open('Mus_musculus.GRCm38.69.dna.chromosome.{}.fa'.format(chr)), "fasta")
    with open('chr{}.txt'.format(chr),'wb') as t_file:
        t_file.write(f.seq.upper().tostring())
        seq_len = len(f.seq.upper().tostring())
    chr_object, created = Chromosome.objects.get_or_create(assembly=mm10, name=chr, sequence_length=seq_len, cyclic=False)


chr1 = Chromosome.objects.filter(assembly=hg19).get(name='1')
chr1.get_abs_path()


# mm9
chr_list = ['10',
'18',
'6',
'11',
'19',
'7',
'12',
'1',
'13',
'8',
'Y',
'2',
'14',
'3',
'9',
'15',
'16',
'4',
'17',
'5',]

# mm10
chr_list = ['1',
'2',
'3',
'4',
'5',
'6',
'7',
'8',
'9',
'10',
'11',
'12',
'13',
'14',
'15',
'16',
'17',
'18',
'19',
'Y',
'X',]
