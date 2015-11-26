# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0002_split_to_apps'),
        ('misc', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    state_ops = [
        migrations.CreateModel(
            name='Assembly',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('friendly_name', models.CharField(max_length=50)),
                ('taxa', models.ForeignKey(to='misc.Taxa')),
            ],
            options={
                'verbose_name_plural': 'Assemblies',
            },
        ),
        migrations.CreateModel(
            name='Chromosome',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('sequence_length', models.IntegerField(null=True)),
                ('cyclic', models.BooleanField()),
                ('assembly', models.ForeignKey(to='genomes.Assembly')),
            ],
        ),
        migrations.CreateModel(
            name='PrimerTail',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('tail', models.CharField(max_length=50, null=True)),
            ],
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
        ),
        migrations.CreateModel(
            name='Sequence',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('length', models.IntegerField()),
                ('sequence', models.TextField()),
                ('hash', models.CharField(unique=True, max_length=32)),
            ],
        ),
        migrations.CreateModel(
            name='Target',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('start_pos', models.IntegerField(db_index=True)),
                ('end_pos', models.IntegerField(db_index=True)),
            ],
        ),
        migrations.CreateModel(
            name='TargetEnrichment',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('amplicon', models.CharField(max_length=500)),
                ('passed_validation', models.NullBooleanField()),
                ('validation_date', models.DateField(null=True, blank=True)),
                ('comment', models.CharField(max_length=50, null=True, blank=True)),
                ('chromosome', models.ForeignKey(to='genomes.Chromosome')),
                ('partner', models.ManyToManyField(to=settings.AUTH_USER_MODEL, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='TargetEnrichmentFailureType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('description', models.TextField(null=True, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='TargetEnrichmentType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('protocol', models.ForeignKey(to='linapp.Protocol', null=True)),
            ],
        ),
        migrations.CreateModel(
            name='TargetType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='Microsatellite',
            fields=[
                ('target_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='genomes.Target')),
                ('repeat_unit_len', models.PositiveIntegerField()),
                ('repeat_unit_type', models.CharField(max_length=50)),
                ('repeat_number', models.DecimalField(null=True, max_digits=5, decimal_places=1)),
            ],
            bases=('genomes.target',),
        ),
        migrations.CreateModel(
            name='Primer',
            fields=[
                ('target_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='genomes.Target')),
                ('strand', models.CharField(max_length=1, null=True, choices=[(b'+', b'Plus'), (b'-', b'Minus')])),
                ('sequence', models.ForeignKey(to='genomes.Sequence')),
                ('tail', models.ForeignKey(to='genomes.PrimerTail', null=True)),
            ],
            bases=('genomes.target',),
        ),
        migrations.CreateModel(
            name='RestrictionSite',
            fields=[
                ('target_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='genomes.Target')),
                ('restriction_type', models.ForeignKey(to='genomes.RestrictionSiteType')),
            ],
            bases=('genomes.target',),
        ),
        migrations.CreateModel(
            name='SNP',
            fields=[
                ('target_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='genomes.Target')),
                ('mutation', models.CharField(max_length=10)),
                ('modified', models.CharField(max_length=10)),
            ],
            bases=('genomes.target',),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='targets',
            field=models.ManyToManyField(related_name='primer_pair', null=True, to='genomes.Target', blank=True),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='type',
            field=models.ForeignKey(to='genomes.TargetEnrichmentType'),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='validation_failure',
            field=models.ForeignKey(to='genomes.TargetEnrichmentFailureType', null=True),
        ),
        migrations.AddField(
            model_name='target',
            name='chromosome',
            field=models.ForeignKey(to='genomes.Chromosome'),
        ),
        migrations.AddField(
            model_name='target',
            name='partner',
            field=models.ManyToManyField(to=settings.AUTH_USER_MODEL, null=True),
        ),
        migrations.AddField(
            model_name='target',
            name='referencevalue',
            field=models.ForeignKey(to='genomes.Sequence'),
        ),
        migrations.AddField(
            model_name='target',
            name='type',
            field=models.ForeignKey(to='genomes.TargetType'),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='left',
            field=models.ForeignKey(related_name='left_primer', to='genomes.Primer'),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='right',
            field=models.ForeignKey(related_name='right_primer', to='genomes.Primer'),
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
        ),
    ]
