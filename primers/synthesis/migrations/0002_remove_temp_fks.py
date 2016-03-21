# -*- coding: utf-8 -*-


from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0004_split_targets_etc'),
        ('synthesis', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='pcr1minusprimer',
            name='old_primer',
        ),
        migrations.RemoveField(
            model_name='pcr1plusprimer',
            name='old_primer',
        ),
        migrations.RemoveField(
            model_name='pcr1withcompanytagminusprimer',
            name='old_primer',
        ),
        migrations.RemoveField(
            model_name='pcr1withcompanytagplusprimer',
            name='old_primer',
        ),
        migrations.RemoveField(
            model_name='targetednotailminusprimer',
            name='old_primer',
        ),
        migrations.RemoveField(
            model_name='targetednotailplusprimer',
            name='old_primer',
        ),
    ]
