# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('parts', '0001_initial'),
        ('workflows', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='BarcodePair',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('left', models.ForeignKey(to='parts.DNABarcode1')),
                ('right', models.ForeignKey(to='parts.DNABarcode2')),
            ],
        ),
        migrations.CreateModel(
            name='WorkFlowCell',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('barcodes', models.ForeignKey(to='workflows.BarcodePair')),
                ('content', models.ForeignKey(to='workflows.CellContent')),
            ],
        ),
    ]
