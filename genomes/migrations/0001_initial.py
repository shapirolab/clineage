# -*- coding: utf-8 -*-


from django.db import migrations, models
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0002_split_to_apps'),
        ('misc', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    state_ops = [
        migrations.CreateModel(
            name='Assembly',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('friendly_name', models.CharField(max_length=50)),
                ('taxa', models.ForeignKey(to='misc.Taxa')),
            ],
            options={
                'verbose_name_plural': 'Assemblies',
            },
        ),
        migrations.CreateModel(
            name='Chromosome',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('sequence_length', models.IntegerField(null=True)),
                ('cyclic', models.BooleanField()),
                ('assembly', models.ForeignKey(to='genomes.Assembly')),
            ],
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
        ),
    ]
