# -*- coding: utf-8 -*-


from django.db import migrations, models
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('linapp', '0002_split_to_apps'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    state_ops = [
        migrations.CreateModel(
            name='CellContentType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='CellContent',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=50, null=True, blank=True)),
                ('comment', models.TextField()),
                ('cell', models.ForeignKey(to='sampling.Cell')),
                ('protocol', models.ForeignKey(blank=True, to='linapp.Protocol', null=True)),
                ('type', models.ForeignKey(to='workflows.CellContentType')),
                ('user', models.ForeignKey(blank=True, to=settings.AUTH_USER_MODEL, null=True)),
            ],
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=state_ops,
        ),
    ]
