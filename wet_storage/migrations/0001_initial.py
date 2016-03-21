# -*- coding: utf-8 -*-


from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0002_split_to_apps'),
    ]

    state_ops = [
        migrations.CreateModel(
            name='Plate',
            fields=[
                ('code', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=200, blank=True)),
                ('barcode', models.CharField(max_length=20, blank=True)),
                ('timestamp', models.DateField(null=True, blank=True)),
                ('state', models.CharField(max_length=20, blank=True)),
                ('lastusedwell', models.CharField(default=b'A1', max_length=4)),
            ],
        ),
        migrations.CreateModel(
            name='PlateContext',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('description', models.CharField(max_length=30, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='PlatePlastica',
            fields=[
                ('code', models.AutoField(serialize=False, primary_key=True)),
                ('description', models.CharField(max_length=30, blank=True)),
                ('rows', models.IntegerField(default=8)),
                ('columns', models.IntegerField(default=12)),
            ],
        ),
        migrations.CreateModel(
            name='PlateStorage',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('inner_location', models.CharField(max_length=100, blank=True)),
                ('notes', models.CharField(max_length=250, blank=True)),
                ('plate', models.ForeignKey(to='wet_storage.Plate')),
            ],
        ),
        migrations.CreateModel(
            name='PlateType',
            fields=[
                ('code', models.AutoField(serialize=False, primary_key=True)),
                ('friendly', models.CharField(max_length=100)),
                ('context', models.ForeignKey(to='wet_storage.PlateContext', null=True)),
                ('plastic', models.ForeignKey(to='wet_storage.PlatePlastica', null=True)),
            ],
        ),
        migrations.CreateModel(
            name='SampleLocation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('well', models.CharField(db_index=True, max_length=3, blank=True)),
                ('object_id', models.PositiveIntegerField(db_index=True)),
                ('volume', models.DecimalField(null=True, max_digits=10, decimal_places=3, blank=True)),
                ('concentration', models.DecimalField(null=True, max_digits=10, decimal_places=5, blank=True)),
                ('timestamp', models.DateTimeField(auto_now=True)),
                ('content_type', models.ForeignKey(to='contenttypes.ContentType')),
                ('plate', models.ForeignKey(to='wet_storage.Plate')),
            ],
        ),
        migrations.CreateModel(
            name='StorageBox',
            fields=[
                ('code', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=100, blank=True)),
                ('barcode', models.CharField(max_length=20, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='StorageType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=100, blank=True)),
                ('temperature', models.DecimalField(null=True, max_digits=5, decimal_places=1, blank=True)),
            ],
        ),
        migrations.AddField(
            model_name='storagebox',
            name='storage_type',
            field=models.ForeignKey(to='wet_storage.StorageType'),
        ),
        migrations.AddField(
            model_name='platestorage',
            name='storageBox',
            field=models.ForeignKey(to='wet_storage.StorageBox'),
        ),
        migrations.AddField(
            model_name='plate',
            name='type',
            field=models.ForeignKey(to='wet_storage.PlateType'),
        ),
        migrations.AlterIndexTogether(
            name='samplelocation',
            index_together=set([('content_type', 'object_id')]),
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
        ),
    ]
