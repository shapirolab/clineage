# -*- coding: utf-8 -*-


import django
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('parts', '0003_add_padlock_parts'),
        ('planning', '0003_auto_20160215_1652'),
        ('synthesis', '0002_remove_temp_fks'),
    ]

    operations = [
        migrations.CreateModel(
            name='OM6OligomixStock',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('left_amp_primer_part1', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='parts.PadlockAmplificationPlusPrimerPart1')),
                ('left_amp_primer_part2', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='parts.PadlockAmplificationPlusPrimerPart2')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='OM6Padlock',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('umi_length', models.PositiveSmallIntegerField()),
                ('backbone', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='parts.Backbone')),
                ('ira1', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='parts.IlluminaReadingAdaptor1')),
                ('ira2', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='parts.IlluminaReadingAdaptor2')),
                ('left_ugs', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='planning.UGSPlus')),
                ('right_ugs', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='planning.UGSMinus')),
            ],
        ),
        migrations.AddField(
            model_name='om6oligomixstock',
            name='padlock',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='synthesis.OM6Padlock'),
        ),
        migrations.AddField(
            model_name='om6oligomixstock',
            name='restriction_enzyme',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='planning.RestrictionEnzyme'),
        ),
        migrations.AddField(
            model_name='om6oligomixstock',
            name='right_amp_primer_part1',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='parts.PadlockAmplificationMinusPrimerPart1'),
        ),
        migrations.AddField(
            model_name='om6oligomixstock',
            name='right_amp_primer_part2',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='parts.PadlockAmplificationMinusPrimerPart2'),
        ),
    ]
