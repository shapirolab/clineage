# -*- coding: utf-8 -*-


from django.db import migrations, models
import primers.strand


class Migration(migrations.Migration):

    dependencies = [
        ('parts', '0002_remove_temp_fks'),
    ]

    operations = [
        migrations.CreateModel(
            name='PadlockAmplificationMinusPrimer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.MinusStrandMixin, primers.strand.BaseStrandMixin),
        ),
        migrations.CreateModel(
            name='PadlockAmplificationPlusPrimer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('_sequence', models.CharField(max_length=250)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model, primers.strand.PlusStrandMixin, primers.strand.BaseStrandMixin),
        ),
    ]
