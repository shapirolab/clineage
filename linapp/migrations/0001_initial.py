# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import mptt.fields
import linapp.models
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('contenttypes', '0002_remove_content_type_name'),
    ]

    operations = [
        migrations.CreateModel(
            name='Algorithm',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('version', models.CharField(max_length=50)),
                ('developers', models.ManyToManyField(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='AlgorithmParameter',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('type', models.IntegerField()),
                ('algorithm', models.ManyToManyField(related_name=b'parameters', to='linapp.Algorithm')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='AlgorithmRun',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('runname', models.CharField(max_length=50, null=True)),
                ('timestamp', models.DateTimeField()),
                ('status', models.IntegerField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='AlgorithmRunParameters',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.CharField(max_length=50)),
                ('parameter', models.ForeignKey(to='linapp.AlgorithmParameter')),
                ('run', models.ForeignKey(to='linapp.AlgorithmRun')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='AlgorithmType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('input', models.CharField(max_length=50)),
                ('output', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Assembly',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('friendly_name', models.CharField(max_length=50)),
            ],
            options={
                'verbose_name_plural': 'Assemblies',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Cell',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('comment', models.TextField(null=True, blank=True)),
                ('classification', models.CharField(max_length=50, null=True, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CellContent',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50, null=True, blank=True)),
                ('seq_ready', models.BooleanField(default=False)),
                ('comment', models.TextField()),
                ('cell', models.ForeignKey(to='linapp.Cell')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CellContentType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CellTreeNode',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('distance', models.IntegerField()),
                ('lft', models.PositiveIntegerField(editable=False, db_index=True)),
                ('rght', models.PositiveIntegerField(editable=False, db_index=True)),
                ('tree_id', models.PositiveIntegerField(editable=False, db_index=True)),
                ('level', models.PositiveIntegerField(editable=False, db_index=True)),
                ('cell', models.ForeignKey(to='linapp.Cell')),
                ('parent', mptt.fields.TreeForeignKey(related_name=b'children', blank=True, to='linapp.CellTreeNode', null=True)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Chromosome',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('sequence_length', models.IntegerField(null=True)),
                ('cyclic', models.BooleanField()),
                ('assembly', models.ForeignKey(to='linapp.Assembly')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Coordinates',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('x', models.DecimalField(max_digits=10, decimal_places=4)),
                ('y', models.DecimalField(max_digits=10, decimal_places=4)),
                ('z', models.DecimalField(max_digits=10, decimal_places=4)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='CorrectedRawData',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('file', models.FilePathField()),
                ('creation', models.ForeignKey(related_name=b'crd', to='linapp.AlgorithmRun')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='DM',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('distance', models.IntegerField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Experiment',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('created_date', models.DateField(auto_now_add=True)),
                ('description', models.TextField()),
                ('is_public', models.BooleanField(default=False)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ExperimentFile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('title', models.CharField(max_length=50)),
                ('file_name', models.CharField(max_length=50)),
                ('file', models.FileField(upload_to=linapp.models.expriment_file_path)),
                ('upload_date', models.DateField()),
                ('description', models.TextField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ExperimentLog',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('date', models.DateTimeField(auto_now=True)),
                ('comment', models.TextField()),
                ('experiment', models.ForeignKey(related_name=b'comments', to='linapp.Experiment')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ExperimentUser',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('experiment', models.ForeignKey(to='linapp.Experiment')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Extraction',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('date', models.DateTimeField(null=True, blank=True)),
                ('comment', models.TextField(null=True, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ExtractionEvent',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=100)),
                ('comment', models.TextField(null=True, blank=True)),
                ('date', models.DateTimeField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='FACSMarker',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='FailedTargetValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('comment', models.TextField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='FileContext',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='GeneticBackground',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='GenSig',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.TextField()),
                ('creation', models.ForeignKey(related_name=b'gensigs', to='linapp.AlgorithmRun')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Individual',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('sex', models.CharField(max_length=1, choices=[(b'M', b'Male'), (b'F', b'Female')])),
                ('born', models.DateTimeField(null=True, blank=True)),
                ('comment', models.TextField(null=True, blank=True)),
                ('sacrificed', models.DateTimeField(null=True, blank=True)),
                ('background', models.ForeignKey(blank=True, to='linapp.GeneticBackground', null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='LineageRole',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('read', models.BooleanField(default=True)),
                ('write', models.BooleanField(default=False)),
                ('delete', models.BooleanField(default=False)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Location',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Machine',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('machineid', models.CharField(max_length=50)),
                ('name', models.CharField(max_length=100, null=True, blank=True)),
                ('ip', models.IPAddressField(null=True, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='MachineType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('company', models.CharField(max_length=50)),
                ('model', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Organ',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Panel',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Plate',
            fields=[
                ('code', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=200, blank=True)),
                ('barcode', models.CharField(max_length=20, blank=True)),
                ('timestamp', models.DateField(null=True, blank=True)),
                ('state', models.CharField(max_length=20, blank=True)),
                ('lastusedwell', models.CharField(default=b'A1', max_length=4)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PlateContext',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('description', models.CharField(max_length=30, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PlatePlastica',
            fields=[
                ('code', models.AutoField(serialize=False, primary_key=True)),
                ('description', models.CharField(max_length=30, blank=True)),
                ('rows', models.IntegerField(default=8)),
                ('columns', models.IntegerField(default=12)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PlateStorage',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('inner_location', models.CharField(max_length=100, blank=True)),
                ('notes', models.CharField(max_length=250, blank=True)),
                ('plate', models.ForeignKey(to='linapp.Plate')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PlateType',
            fields=[
                ('code', models.AutoField(serialize=False, primary_key=True)),
                ('friendly', models.CharField(max_length=100)),
                ('context', models.ForeignKey(to='linapp.PlateContext', null=True)),
                ('plastic', models.ForeignKey(to='linapp.PlatePlastica', null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PrimersMultiplex',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=20)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='PrimerTail',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('tail', models.CharField(max_length=50, null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Protocol',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('initials', models.CharField(max_length=10)),
                ('name', models.CharField(max_length=100)),
                ('abstract', models.TextField()),
                ('fulldescription', models.TextField()),
                ('kit', models.CharField(max_length=100, null=True, blank=True)),
                ('file', models.FilePathField(null=True, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ProtocolType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='RawData',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('file', models.FilePathField(null=True, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='RestrictionSiteType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('sequence', models.CharField(max_length=50)),
                ('cut_delta', models.IntegerField()),
                ('sticky_bases', models.IntegerField()),
                ('sequence_len', models.PositiveIntegerField()),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='SampleComposition',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='SampleLocation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('well', models.CharField(db_index=True, max_length=3, blank=True)),
                ('object_id', models.PositiveIntegerField(db_index=True)),
                ('volume', models.DecimalField(null=True, max_digits=10, decimal_places=3, blank=True)),
                ('concentration', models.DecimalField(null=True, max_digits=10, decimal_places=5, blank=True)),
                ('timestamp', models.DateTimeField(auto_now=True)),
                ('content_type', models.ForeignKey(to='contenttypes.ContentType')),
                ('plate', models.ForeignKey(to='linapp.Plate')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='SampleStatus',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
                'verbose_name': 'Sample status',
                'verbose_name_plural': 'Samples status',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='SamplingEvent',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('date', models.DateField(null=True, blank=True)),
                ('comment', models.TextField(null=True, blank=True)),
                ('attachment', models.FileField(null=True, upload_to=linapp.models.sampling_event_path, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='LaserCapture',
            fields=[
                ('samplingevent_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='linapp.SamplingEvent')),
                ('coordinates', models.ForeignKey(blank=True, to='linapp.Coordinates', null=True)),
            ],
            options={
            },
            bases=('linapp.samplingevent',),
        ),
        migrations.CreateModel(
            name='FACS',
            fields=[
                ('samplingevent_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='linapp.SamplingEvent')),
                ('marker', models.ForeignKey(to='linapp.FACSMarker')),
            ],
            options={
            },
            bases=('linapp.samplingevent',),
        ),
        migrations.CreateModel(
            name='CellSelector',
            fields=[
                ('samplingevent_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='linapp.SamplingEvent')),
                ('coordinates', models.ForeignKey(blank=True, to='linapp.Coordinates', null=True)),
            ],
            options={
            },
            bases=('linapp.samplingevent',),
        ),
        migrations.CreateModel(
            name='Sequence',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('length', models.IntegerField()),
                ('sequence', models.TextField()),
                ('hash', models.CharField(unique=True, max_length=32)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='SequenceDistribution',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('count', models.IntegerField()),
                ('failed', models.ForeignKey(to='linapp.FailedTargetValue', null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Sequencing',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=100)),
                ('date', models.DateField()),
                ('data', models.ForeignKey(related_name=b'sequencing_event', blank=True, to='linapp.RawData', null=True)),
                ('machine', models.ForeignKey(to='linapp.Machine')),
                ('protocol', models.ForeignKey(to='linapp.Protocol')),
                ('samples', models.ManyToManyField(to='linapp.CellContent')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='StorageBox',
            fields=[
                ('code', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=100, blank=True)),
                ('barcode', models.CharField(max_length=20, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='StorageType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=100, blank=True)),
                ('temperature', models.DecimalField(null=True, max_digits=5, decimal_places=1, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Target',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('start_pos', models.IntegerField(db_index=True)),
                ('end_pos', models.IntegerField(db_index=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='SNP',
            fields=[
                ('target_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='linapp.Target')),
                ('mutation', models.CharField(max_length=10)),
                ('modified', models.CharField(max_length=10)),
            ],
            options={
            },
            bases=('linapp.target',),
        ),
        migrations.CreateModel(
            name='RestrictionSite',
            fields=[
                ('target_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='linapp.Target')),
                ('restriction_type', models.ForeignKey(to='linapp.RestrictionSiteType')),
            ],
            options={
            },
            bases=('linapp.target',),
        ),
        migrations.CreateModel(
            name='Primer',
            fields=[
                ('target_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='linapp.Target')),
                ('strand', models.CharField(max_length=1, null=True, choices=[(b'+', b'Plus'), (b'-', b'Minus')])),
                ('sequence', models.ForeignKey(to='linapp.Sequence')),
                ('tail', models.ForeignKey(to='linapp.PrimerTail', null=True)),
            ],
            options={
            },
            bases=('linapp.target',),
        ),
        migrations.CreateModel(
            name='Microsatellite',
            fields=[
                ('target_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='linapp.Target')),
                ('repeat_unit_len', models.PositiveIntegerField()),
                ('repeat_unit_type', models.CharField(max_length=50)),
                ('repeat_number', models.DecimalField(null=True, max_digits=5, decimal_places=1)),
            ],
            options={
            },
            bases=('linapp.target',),
        ),
        migrations.CreateModel(
            name='TargetAnalysis',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('creation', models.ForeignKey(related_name=b'targetsdistributions', to='linapp.AlgorithmRun')),
                ('distribution', models.ForeignKey(to='linapp.SequenceDistribution')),
                ('sequencingdata', models.ForeignKey(to='linapp.CorrectedRawData')),
                ('target', models.ForeignKey(to='linapp.Target')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TargetEnrichment',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('amplicon', models.CharField(max_length=500)),
                ('passed_validation', models.NullBooleanField()),
                ('validation_date', models.DateField(null=True, blank=True)),
                ('comment', models.CharField(max_length=50, null=True, blank=True)),
                ('chromosome', models.ForeignKey(to='linapp.Chromosome')),
                ('left', models.ForeignKey(related_name=b'left_primer', to='linapp.Primer')),
                ('partner', models.ManyToManyField(to=settings.AUTH_USER_MODEL)),
                ('right', models.ForeignKey(related_name=b'right_primer', to='linapp.Primer')),
                ('targets', models.ManyToManyField(related_name=b'primer_pair', to='linapp.Target')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TargetEnrichmentFailureType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('description', models.TextField(null=True, blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TargetEnrichmentType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('protocol', models.ForeignKey(to='linapp.Protocol', null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TargetType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TargetVariant',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.PositiveIntegerField()),
                ('creation', models.ForeignKey(related_name=b'targetsvariants', to='linapp.AlgorithmRun')),
                ('distribution', models.ManyToManyField(to='linapp.SequenceDistribution')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='TargetVariantType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Taxa',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('taxonomy_id', models.IntegerField()),
                ('rank', models.CharField(max_length=50)),
                ('parent', models.IntegerField(null=True, blank=True)),
                ('friendly_name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Tissue',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='UserProfile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('institute', models.CharField(max_length=50)),
                ('comment', models.TextField()),
                ('user', models.OneToOneField(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='UserReport',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('creation_date', models.DateField(auto_now_add=True)),
                ('cells', models.ManyToManyField(to='linapp.Cell')),
                ('individual', models.ManyToManyField(to='linapp.Individual')),
                ('partner', models.ForeignKey(to=settings.AUTH_USER_MODEL, null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='targetvariant',
            name='type',
            field=models.ForeignKey(to='linapp.TargetVariantType'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='type',
            field=models.ForeignKey(to='linapp.TargetEnrichmentType'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='validation_failure',
            field=models.ForeignKey(to='linapp.TargetEnrichmentFailureType', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='target',
            name='chromosome',
            field=models.ForeignKey(to='linapp.Chromosome'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='target',
            name='partner',
            field=models.ManyToManyField(to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='target',
            name='referencevalue',
            field=models.ForeignKey(to='linapp.Sequence'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='target',
            name='type',
            field=models.ForeignKey(to='linapp.TargetType'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='storagebox',
            name='storage_type',
            field=models.ForeignKey(to='linapp.StorageType'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='sequencedistribution',
            name='target',
            field=models.ManyToManyField(to='linapp.Target', through='linapp.TargetAnalysis'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='sequencedistribution',
            name='value',
            field=models.ForeignKey(to='linapp.Sequence'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='samplingevent',
            name='extraction',
            field=models.ForeignKey(to='linapp.Extraction'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='samplingevent',
            name='user',
            field=models.ForeignKey(blank=True, to=settings.AUTH_USER_MODEL, null=True),
            preserve_default=True,
        ),
        migrations.AlterIndexTogether(
            name='samplelocation',
            index_together=set([('content_type', 'object_id')]),
        ),
        migrations.AddField(
            model_name='rawdata',
            name='sequencing',
            field=models.ForeignKey(to='linapp.Sequencing'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='rawdata',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='protocol',
            name='type',
            field=models.ForeignKey(to='linapp.ProtocolType'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='primersmultiplex',
            name='primers',
            field=models.ManyToManyField(to='linapp.TargetEnrichment'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='platestorage',
            name='storageBox',
            field=models.ForeignKey(to='linapp.StorageBox'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='plate',
            name='type',
            field=models.ForeignKey(to='linapp.PlateType'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='panel',
            name='targets',
            field=models.ManyToManyField(related_name=b'panels', to='linapp.TargetEnrichment'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='machine',
            name='type',
            field=models.ForeignKey(to='linapp.MachineType'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='individual',
            name='location',
            field=models.ForeignKey(blank=True, to='linapp.Location', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='individual',
            name='partner',
            field=models.ForeignKey(blank=True, to=settings.AUTH_USER_MODEL, null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='individual',
            name='taxa',
            field=models.ForeignKey(to='linapp.Taxa'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='gensig',
            name='variants',
            field=models.ManyToManyField(to='linapp.TargetVariant'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='extractionevent',
            name='individual',
            field=models.ForeignKey(to='linapp.Individual'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='extractionevent',
            name='location',
            field=models.ForeignKey(blank=True, to='linapp.Location', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='extractionevent',
            name='user_documented',
            field=models.ForeignKey(related_name=b'+', to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='extractionevent',
            name='user_performed',
            field=models.ForeignKey(related_name=b'+', to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='extraction',
            name='extraction_event',
            field=models.ForeignKey(to='linapp.ExtractionEvent'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='extraction',
            name='organ',
            field=models.ForeignKey(to='linapp.Organ'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='extraction',
            name='tissue',
            field=models.ForeignKey(to='linapp.Tissue'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='experimentuser',
            name='role',
            field=models.ForeignKey(to='linapp.LineageRole'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='experimentuser',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='experimentfile',
            name='context',
            field=models.ForeignKey(to='linapp.FileContext'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='experimentfile',
            name='experiment',
            field=models.ForeignKey(to='linapp.Experiment'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='experimentfile',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='experiment',
            name='users',
            field=models.ManyToManyField(related_name=b'experiments', through='linapp.ExperimentUser', to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='dm',
            name='cell1',
            field=models.ForeignKey(related_name=b'+', to='linapp.GenSig'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='dm',
            name='cell2',
            field=models.ForeignKey(related_name=b'+', to='linapp.GenSig'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='dm',
            name='creation',
            field=models.ForeignKey(related_name=b'dms', to='linapp.AlgorithmRun'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='correctedrawdata',
            name='rawdata',
            field=models.ForeignKey(to='linapp.RawData'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cellcontent',
            name='panel',
            field=models.ForeignKey(blank=True, to='linapp.Panel', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cellcontent',
            name='parent',
            field=models.ForeignKey(blank=True, to='linapp.CellContent', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cellcontent',
            name='protocol',
            field=models.ForeignKey(blank=True, to='linapp.Protocol', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cellcontent',
            name='type',
            field=models.ForeignKey(to='linapp.CellContentType'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cellcontent',
            name='user',
            field=models.ForeignKey(blank=True, to=settings.AUTH_USER_MODEL, null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cell',
            name='composition',
            field=models.ForeignKey(to='linapp.SampleComposition'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cell',
            name='experiment',
            field=models.ManyToManyField(related_name=b'cells', to='linapp.Experiment'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cell',
            name='individual',
            field=models.ForeignKey(to='linapp.Individual'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cell',
            name='sampling',
            field=models.ForeignKey(to='linapp.SamplingEvent', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='cell',
            name='status',
            field=models.ForeignKey(blank=True, to='linapp.SampleStatus', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='assembly',
            name='taxa',
            field=models.ForeignKey(to='linapp.Taxa'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='algorithmrun',
            name='ExtraFiles',
            field=models.ManyToManyField(to='linapp.ExperimentFile'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='algorithmrun',
            name='algorithm',
            field=models.ForeignKey(related_name=b'runs', to='linapp.Algorithm'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='algorithmrun',
            name='parameters',
            field=models.ManyToManyField(to='linapp.AlgorithmParameter', through='linapp.AlgorithmRunParameters'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='algorithmrun',
            name='user',
            field=models.ForeignKey(to=settings.AUTH_USER_MODEL),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='algorithm',
            name='type',
            field=models.ForeignKey(to='linapp.AlgorithmType'),
            preserve_default=True,
        ),
    ]
