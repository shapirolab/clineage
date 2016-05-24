# -*- coding: utf-8 -*-


from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('multiplexes', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='PCR1Panel',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=20)),
                ('mpxs', models.ManyToManyField(to='multiplexes.PCR1Multiplex')),
            ],
        ),
    ]
