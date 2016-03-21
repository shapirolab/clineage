# -*- coding: utf-8 -*-


from django.db import migrations, models
from django.conf import settings
import primers.strand


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('genomes', '0002_create_dnaslice'),
        ('linapp', '0003_add_fks_after_split'),
    ]

    operations = [
        migrations.CreateModel(
            name='RestrictionEnzyme',
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
            name='RestrictionSite',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('enzyme', models.ForeignKey(related_name='sites', to='planning.RestrictionEnzyme')),
                ('slice', models.ForeignKey(to='genomes.DNASlice')),
            ],
        ),
        migrations.CreateModel(
            name='Target',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='TargetEnrichment',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('chromosome', models.ForeignKey(to='genomes.Chromosome')),
                ('planning_version', models.IntegerField()),
            ],
        ),
        migrations.CreateModel(
            name='UGSMinus',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('slice', models.ForeignKey(to='genomes.DNASlice')),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.MinusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='UGSPlus',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('slice', models.ForeignKey(to='genomes.DNASlice')),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.PlusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='Microsatellite',
            fields=[
                ('target_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='planning.Target')),
                ('repeat_unit_len', models.PositiveIntegerField()),
                ('repeat_unit_type', models.CharField(max_length=50)),
                ('repeat_number', models.DecimalField(null=True, max_digits=5, decimal_places=1)),
            ],
            bases=('planning.target',),
        ),
        migrations.CreateModel(
            name='SNP',
            fields=[
                ('target_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='planning.Target')),
                ('mutation', models.CharField(max_length=10, null=True)),
                ('modified', models.CharField(max_length=10, null=True)),
            ],
            bases=('planning.target',),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='left',
            field=models.ForeignKey(to='planning.UGSPlus'),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='partner',
            field=models.ManyToManyField(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='right',
            field=models.ForeignKey(to='planning.UGSMinus'),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='targets',
            field=models.ManyToManyField(to='planning.Target'),
        ),
        migrations.AddField(
            model_name='target',
            name='partner',
            field=models.ManyToManyField(to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='target',
            name='slice',
            field=models.ForeignKey(to='genomes.DNASlice'),
        ),
        # Temporary fields to assist with the data migration.
        migrations.AddField(
            model_name='ugsplus',
            name='old_primer',
            field=models.ForeignKey(to='linapp.Primer', null=True),
        ),
        migrations.AddField(
            model_name='ugsminus',
            name='old_primer',
            field=models.ForeignKey(to='linapp.Primer', null=True),
        ),
        migrations.AddField(
            model_name='target',
            name='old_target',
            field=models.ForeignKey(to='linapp.Target', null=True),
        ),
        migrations.AddField(
            model_name='targetenrichment',
            name='old_te',
            field=models.ForeignKey(to='linapp.TargetEnrichment', null=True),
        ),
    ]
