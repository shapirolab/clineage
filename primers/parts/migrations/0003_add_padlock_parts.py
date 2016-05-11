# -*- coding: utf-8 -*-


from django.db import migrations, models
import primers.strand


class Migration(migrations.Migration):

    dependencies = [
        ('parts', '0002_remove_temp_fks'),
    ]

    operations = [
        migrations.CreateModel(
            name='Backbone',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='PadlockAmplificationMinusPrimerPart1',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.MinusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='PadlockAmplificationMinusPrimerPart2',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.MinusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='PadlockAmplificationPlusPrimerPart1',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.PlusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='PadlockAmplificationPlusPrimerPart2',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.PlusStrandMixin, primers.strand.BaseStrandMixin),
        ),
    ]
