# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('workflows', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('lib_prep', '0003_kill'),
    ]

    state_ops = [
        migrations.CreateModel(
            name='MachineType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('company', models.CharField(max_length=50)),
                ('model', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='Machine',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('machineid', models.CharField(max_length=50)),
                ('name', models.CharField(max_length=100, null=True, blank=True)),
                ('ip', models.IPAddressField(null=True, blank=True)),
                ('type', models.ForeignKey(to='runs.MachineType')),
            ],
        ),
        migrations.CreateModel(
            name='NGSRun',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(unique=True, max_length=100)),
                ('date', models.DateField()),
                ('machine', models.ForeignKey(to='runs.Machine')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
        ),
        migrations.AddField(
            model_name='ngsrun',
            name='directory',
            field=models.FilePathField(null=True),
        ),
        migrations.CreateModel(
            name='DemultiplexingScheme',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('description', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='DemultiplexedReads',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('directory', models.FilePathField(null=True)),
                ('ngs_run', models.ForeignKey(to='runs.NGSRun')),
                ('demux_scheme', models.ForeignKey(to='runs.DemultiplexingScheme')),
            ],
        ),
        migrations.CreateModel(
            name='MergingScheme',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('description', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='MergedReads',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('directory', models.FilePathField(null=True)),
                ('demux_reads', models.ForeignKey(to='runs.DemultiplexedReads')),
                ('merge_scheme', models.ForeignKey(to='runs.MergingScheme')),
            ],
        ),
    ]
