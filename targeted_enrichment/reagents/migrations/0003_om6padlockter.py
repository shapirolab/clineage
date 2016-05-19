# -*- coding: utf-8 -*-


from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('planning', '0003_auto_20160215_1652'),
        ('synthesis', '0003_om6padlock'),
        ('reagents', '0002_remove_temp_fks'),
    ]

    operations = [
        migrations.CreateModel(
            name='OM6PadlockTERBase',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('passed_validation', models.NullBooleanField()),
                ('validation_date', models.DateField(blank=True, null=True)),
                ('comment', models.CharField(blank=True, max_length=50, null=True)),
                ('te', models.ForeignKey(to='planning.TargetEnrichment')),
                ('validation_failure', models.ForeignKey(to='reagents.TargetEnrichmentFailureType', null=True)),
                ('old_adam_te_pk', models.PositiveIntegerField(null=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='OM6PadlockTER',
            fields=[
                ('om6padlockterbase_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='reagents.OM6PadlockTERBase')),
                ('padlock', models.ForeignKey(to='synthesis.OM6Padlock')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='OM6PadlockTERDeprecated',
            fields=[
                ('om6padlockterbase_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='reagents.OM6PadlockTERBase')),
                ('padlock', models.ForeignKey(to='synthesis.OM6PadlockDeprecated')),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
