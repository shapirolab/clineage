# -*- coding: utf-8 -*-


from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('genomes', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='DNASlice',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('start_pos', models.IntegerField(db_index=True)),
                ('end_pos', models.IntegerField(db_index=True)),
                ('_sequence', models.CharField(max_length=300, default=None, null=True)),
                ('chromosome', models.ForeignKey(to='genomes.Chromosome')),
            ],
        ),
    ]
