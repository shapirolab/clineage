# -*- coding: utf-8 -*-


from django.db import migrations, models

from utils.content_type_operation import AlterContentType


class Migration(migrations.Migration):

    dependencies = [
        ('sampling', '0001_initial'),
        ('parts', '0003_padlockamplificationminusprimer_padlockamplificationplusprimer'),
        ('multiplexes', '0002_pcr1multiplexcollection'),
        ('workflows', '0002_subclass_cell_content_protocol'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='CellContent',
            new_name='AmplifiedContent',
        ),
        AlterContentType(
            from_model='cellcontent',
            to_app='workflows',
            to_model='amplifiedcontent',
        ),
        migrations.RemoveField(
            model_name='AmplifiedContent',
            name='type',
        ),
        migrations.RemoveField(
            model_name='AmplifiedContent',
            name='user',
        ),
        migrations.CreateModel(
            name='BarcodedContent',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
        ),
        migrations.CreateModel(
            name='BarcodePair',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('left', models.ForeignKey(to='parts.DNABarcode1')),
                ('right', models.ForeignKey(to='parts.DNABarcode2')),
            ],
        ),
        migrations.CreateModel(
            name='Library',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
        ),
        migrations.CreateModel(
            name='MagicalPCR1BarcodedContent',
            fields=[
                ('barcodedcontent_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='workflows.BarcodedContent')),
                ('content', models.ForeignKey(to='workflows.AmplifiedContent')),
            ],
            bases=('workflows.barcodedcontent',),
        ),
        migrations.CreateModel(
            name='MagicalPCR1Library',
            fields=[
                ('library_ptr', models.OneToOneField(parent_link=True, auto_created=True, primary_key=True, serialize=False, to='workflows.Library')),
                ('mpx_collection', models.ForeignKey(to='multiplexes.PCR1MultiplexCollection')),
            ],
            bases=('workflows.library',),
        ),
        migrations.DeleteModel(
            name='CellContentType',
        ),
        migrations.AddField(
            model_name='barcodedcontent',
            name='barcodes',
            field=models.ForeignKey(to='workflows.BarcodePair'),
        ),
        migrations.AddField(
            model_name='magicalpcr1barcodedcontent',
            name='library',
            field=models.ForeignKey(to='workflows.MagicalPCR1Library'),
        ),
    ]
