# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('runs', '0002_auto_20160215_1652'),
        ('parts', '0003_padlockamplificationminusprimer_padlockamplificationplusprimer'),
        ('workflows', '0003_auto_20160215_1652'),
    ]

    operations = [
        migrations.AddField(
            model_name='ngsrun',
            name='libraries',
            field=models.ManyToManyField(to='workflows.Library'),
        ),
        migrations.AddField(
            model_name='ngskit',
            name='reading_adaptor1',
            field=models.ForeignKey(to='parts.IlluminaReadingAdaptor1'),
        ),
        migrations.AddField(
            model_name='ngskit',
            name='reading_adaptor2',
            field=models.ForeignKey(to='parts.IlluminaReadingAdaptor2'),
        ),
        migrations.AddField(
            model_name='mergedreads',
            name='demux_read',
            field=models.ForeignKey(to='runs.DemultiplexedReads'),
        ),
        migrations.AddField(
            model_name='mergedreads',
            name='merge_scheme',
            field=models.ForeignKey(to='runs.MergingScheme'),
        ),
        migrations.AddField(
            model_name='demultiplexing',
            name='demux_scheme',
            field=models.ForeignKey(to='runs.DemultiplexingScheme'),
        ),
        migrations.AddField(
            model_name='demultiplexing',
            name='ngs_run',
            field=models.ForeignKey(to='runs.NGSRun'),
        ),
        migrations.AddField(
            model_name='demultiplexedreads',
            name='barcoded_content',
            field=models.ForeignKey(to='workflows.BarcodedContent'),
        ),
        migrations.AddField(
            model_name='demultiplexedreads',
            name='demux',
            field=models.ForeignKey(to='runs.Demultiplexing'),
        ),
        migrations.AddField(
            model_name='demultiplexedreads',
            name='library',
            field=models.ForeignKey(to='workflows.Library'),
        ),
        migrations.AddField(
            model_name='ngsrun',
            name='kit',
            field=models.ForeignKey(to='runs.NGSKit', null=True),
        ),
    ]
