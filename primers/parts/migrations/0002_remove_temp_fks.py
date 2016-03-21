# -*- coding: utf-8 -*-


from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0004_split_targets_etc'),
        ('parts', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='illuminareadingadaptor1cuts',
            name='old_tail',
        ),
        migrations.RemoveField(
            model_name='illuminareadingadaptor2cuts',
            name='old_tail',
        ),
    ]
