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
            name='OM6PadlockTER',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('passed_validation', models.NullBooleanField()),
                ('validation_date', models.DateField(blank=True, null=True)),
                ('comment', models.CharField(blank=True, max_length=50, null=True)),
                ('padlock', models.ForeignKey(to='synthesis.OM6Padlock')),
                ('te', models.ForeignKey(to='planning.TargetEnrichment')),
                ('validation_failure', models.ForeignKey(to='reagents.TargetEnrichmentFailureType', null=True)),
                ('old_adam_te_pk', models.PositiveIntegerField(null=True)),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
