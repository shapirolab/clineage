from django.db import migrations, models
import django.db.models.deletion

rename_view_contains = "rename table genomes_dnaslice_contains to genomes_dnaslice_contains_view;"
rename_view_overlaps = "rename table genomes_dnaslice_overlaps to genomes_dnaslice_overlaps_view;"
undo_rename_view_contains = "rename table genomes_dnaslice_contains_view to genomes_dnaslice_contains;"
undo_rename_view_overlaps = "rename table genomes_dnaslice_overlaps_view to genomes_dnaslice_overlaps;"

create_overlaps_table = """create table genomes_dnaslice_overlaps(
                            slice1_id int(11) not null,
                            slice2_id int(11) not null
                            );"""

create_contains_table = """ create table genomes_dnaslice_contains(
                            inner_id int(11) not null,
                            outer_id int(11) not null
                            ) ; """

populate_overlaps = """ insert into genomes_dnaslice_overlaps
                        SELECT
                          s1.id AS slice1_id, s2.id AS slice2_id
                        FROM
                          genomes_dnaslice s1
                        JOIN
                          genomes_dnaslice s2 ON (s1.chromosome_id = s2.chromosome_id)
                        WHERE ((s1.start_pos <= s2.end_pos)
                        AND (s1.end_pos >= s2.start_pos)
                        AND (s1.id <> s2.id));"""

populate_contains = """ insert into genomes_dnaslice_contains
                        SELECT
                          s1.id AS inner_id, s2.id AS outer_id
                        FROM
                          genomes_dnaslice s1
                        JOIN
                          genomes_dnaslice s2 ON (s1.chromosome_id = s2.chromosome_id)
                        WHERE ((s1.start_pos >= s2.start_pos)
                        AND (s1.end_pos <= s2.end_pos)
                        AND (s1.id <> s2.id));"""
class Migration(migrations.Migration):
    dependencies = [
        ('genomes', '0008_delete_restrictionsites_dnaslices')
    ]

    operations = [
        migrations.RunSQL(
            sql=rename_view_contains,
            reverse_sql=undo_rename_view_contains,
        ),
        migrations.RunSQL(
            sql=rename_view_overlaps,
            reverse_sql=undo_rename_view_overlaps,
        ),

        migrations.RunSQL(
            sql = create_contains_table,
            reverse_sql = "drop table genomes_dnaslice_contains;",
        ),
        migrations.RunSQL(
            sql=create_overlaps_table,
            reverse_sql = "drop table genomes_dnaslice_overlaps;"
        ),


        migrations.RunSQL(
            sql = populate_overlaps,
            reverse_sql = "truncate table genomes_dnaslice_overlaps;",
        ),

        migrations.RunSQL(
            sql = populate_contains,
            reverse_sql = "truncate table genomes_dnaslice_contains;",
        ),

        migrations.RunSQL(
            sql = "alter table genomes_dnaslice_overlaps add foreign key g_dnaslice_overlaps_slice1_id_fk (slice1_id) references genomes_dnaslice(id) on delete cascade on update cascade;",
            reverse_sql = "alter table genomes_dnaslice_overlaps drop foreign key g_dnaslice_overlaps_slice1_id_fk;",
        ),

        migrations.RunSQL(
            sql = "alter table genomes_dnaslice_overlaps add foreign key g_dnaslice_overlaps_slice2_id_fk (slice2_id) references genomes_dnaslice(id) on delete cascade on update cascade;",
            reverse_sql = "alter table genomes_dnaslice_overlaps drop foreign key g_dnaslice_overlaps_slice2_id_fk;",
        ),

        migrations.RunSQL(
            sql = "alter table genomes_dnaslice_contains add foreign key g_dnaslice_contains_in_id_fk (inner_id) references genomes_dnaslice(id) on delete cascade on update cascade;",
            reverse_sql = "alter table genomes_dnaslice_contains drop foreign key g_dnaslice_contains_in_id_fk;",
        ),

        migrations.RunSQL(
            sql = "alter table genomes_dnaslice_contains add foreign key g_dnaslice_contains_out_id_fk (outer_id) references genomes_dnaslice(id) on delete cascade on update cascade;",
            reverse_sql = "alter table genomes_dnaslice_contains drop foreign key g_dnaslice_contains_out_id_fk;"
        )
    ]
