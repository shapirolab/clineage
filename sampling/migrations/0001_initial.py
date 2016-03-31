# -*- coding: utf-8 -*-


from django.db import migrations, models
import sampling.models
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0002_split_to_apps'),
        ('misc', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    state_ops = [
        migrations.CreateModel(
            name='Cell',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('comment', models.TextField(null=True, blank=True)),
                ('classification', models.CharField(max_length=50, null=True, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='Coordinates',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('x', models.DecimalField(max_digits=10, decimal_places=4)),
                ('y', models.DecimalField(max_digits=10, decimal_places=4)),
                ('z', models.DecimalField(max_digits=10, decimal_places=4)),
            ],
        ),
        migrations.CreateModel(
            name='Extraction',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('date', models.DateTimeField(null=True, blank=True)),
                ('comment', models.TextField(null=True, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='ExtractionEvent',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=100)),
                ('comment', models.TextField(null=True, blank=True)),
                ('date', models.DateTimeField()),
            ],
        ),
        migrations.CreateModel(
            name='FACSMarker',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='GeneticBackground',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='Individual',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('sex', models.CharField(max_length=1, choices=[(b'M', b'Male'), (b'F', b'Female')])),
                ('born', models.DateTimeField(null=True, blank=True)),
                ('comment', models.TextField(null=True, blank=True)),
                ('sacrificed', models.DateTimeField(null=True, blank=True)),
                ('background', models.ForeignKey(blank=True, to='sampling.GeneticBackground', null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Location',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='Organ',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='SampleComposition',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='SampleStatus',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
            options={
                'verbose_name': 'Sample status',
                'verbose_name_plural': 'Samples status',
            },
        ),
        migrations.CreateModel(
            name='SamplingEvent',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
                ('date', models.DateField(null=True, blank=True)),
                ('comment', models.TextField(null=True, blank=True)),
                ('attachment', models.FileField(null=True, upload_to=sampling.models.sampling_event_path, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='Tissue',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='CellSelector',
            fields=[
                ('samplingevent_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='sampling.SamplingEvent')),
                ('coordinates', models.ForeignKey(blank=True, to='sampling.Coordinates', null=True)),
            ],
            bases=('sampling.samplingevent',),
        ),
        migrations.CreateModel(
            name='FACS',
            fields=[
                ('samplingevent_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='sampling.SamplingEvent')),
                ('marker', models.ForeignKey(to='sampling.FACSMarker')),
            ],
            bases=('sampling.samplingevent',),
        ),
        migrations.CreateModel(
            name='LaserCapture',
            fields=[
                ('samplingevent_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='sampling.SamplingEvent')),
                ('coordinates', models.ForeignKey(blank=True, to='sampling.Coordinates', null=True)),
            ],
            bases=('sampling.samplingevent',),
        ),
        migrations.AddField(
            model_name='samplingevent',
            name='extraction',
            field=models.ForeignKey(to='sampling.Extraction'),
        ),
        migrations.AddField(
            model_name='samplingevent',
            name='user',
            field=models.ForeignKey(blank=True, to=settings.AUTH_USER_MODEL, null=True),
        ),
        migrations.AddField(
            model_name='individual',
            name='location',
            field=models.ForeignKey(blank=True, to='sampling.Location', null=True),
        ),
        migrations.AddField(
            model_name='individual',
            name='partner',
            field=models.ForeignKey(blank=True, to=settings.AUTH_USER_MODEL, null=True),
        ),
        migrations.AddField(
            model_name='individual',
            name='taxa',
            field=models.ForeignKey(to='misc.Taxa'),
        ),
        migrations.AddField(
            model_name='extractionevent',
            name='individual',
            field=models.ForeignKey(to='sampling.Individual'),
        ),
        migrations.AddField(
            model_name='extractionevent',
            name='location',
            field=models.ForeignKey(blank=True, to='sampling.Location', null=True),
        ),
        migrations.AddField(
            model_name='extractionevent',
            name='user_documented',
            field=models.ForeignKey(related_name='+', to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='extractionevent',
            name='user_performed',
            field=models.ForeignKey(related_name='+', to=settings.AUTH_USER_MODEL),
        ),
        migrations.AddField(
            model_name='extraction',
            name='extraction_event',
            field=models.ForeignKey(to='sampling.ExtractionEvent'),
        ),
        migrations.AddField(
            model_name='extraction',
            name='organ',
            field=models.ForeignKey(to='sampling.Organ'),
        ),
        migrations.AddField(
            model_name='extraction',
            name='tissue',
            field=models.ForeignKey(to='sampling.Tissue'),
        ),
        migrations.AddField(
            model_name='cell',
            name='composition',
            field=models.ForeignKey(to='sampling.SampleComposition'),
        ),
        migrations.AddField(
            model_name='cell',
            name='individual',
            field=models.ForeignKey(to='sampling.Individual'),
        ),
        migrations.AddField(
            model_name='cell',
            name='sampling',
            field=models.ForeignKey(to='sampling.SamplingEvent', null=True),
        ),
        migrations.AddField(
            model_name='cell',
            name='status',
            field=models.ForeignKey(blank=True, to='sampling.SampleStatus', null=True),
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
        ),
    ]
