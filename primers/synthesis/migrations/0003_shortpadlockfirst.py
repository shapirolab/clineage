# -*- coding: utf-8 -*-


from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('parts', '0003_padlockamplificationminusprimer_padlockamplificationplusprimer'),
        ('planning', '0003_auto_20160215_1652'),
        ('synthesis', '0002_remove_temp_fks'),
    ]

    operations = [
        migrations.CreateModel(
            name='ShortPadlockFirst',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('irac1', models.ForeignKey(to='parts.IlluminaReadingAdaptor1Cuts')),
                ('irac2', models.ForeignKey(to='parts.IlluminaReadingAdaptor2Cuts')),
                ('left_amp_primer', models.ForeignKey(to='parts.PadlockAmplificationPlusPrimer')),
                ('left_ugs', models.ForeignKey(to='planning.UGSPlus')),
                ('restriction_enzyme', models.ForeignKey(to='planning.RestrictionEnzyme')),
                ('right_amp_primer', models.ForeignKey(to='parts.PadlockAmplificationMinusPrimer')),
                ('right_ugs', models.ForeignKey(to='planning.UGSMinus')),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
