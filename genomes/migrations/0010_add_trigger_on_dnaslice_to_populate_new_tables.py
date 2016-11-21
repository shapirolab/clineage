from django.db import migrations, models
import django.db.models.deletion


trigger_create = """create trigger genome_dnasalice_insert_trg after insert on genomes_dnaslice
                    for each row
                    begin
                        declare new_id int(11);
                        declare new_chromosome_id int(11);
                        declare new_start_pos int(11);
                        declare new_end_pos int(11);
                        declare slice1 int(11);
                        declare slice2 int(11);
                        set new_id = NEW.id;
                         set new_chromosome_id = new.chromosome_id;
                         set new_start_pos = new.start_pos;
                         set new_end_pos = new.end_pos;

                        insert into genomes_dnaslice_contains (
                            select s1.id AS inner_id,new_id AS outer_id from
                                genomes_dnaslice s1 where
                                s1.chromosome_id = new_chromosome_id and
                                s1.start_pos >= new_start_pos and
                                s1.end_pos <= new_end_pos and
                                s1.id <> new_id
                            );

                        /* symetrical check*/
                        insert into genomes_dnaslice_contains (
                            select new_id AS inner_id,s1.id AS outer_id from
                                genomes_dnaslice s1 where
                                s1.chromosome_id = new_chromosome_id and
                                new_start_pos >= s1.start_pos and
                                new_end_pos <= s1.end_pos and
                                s1.id <> new_id
                            );


                        insert into genomes_dnaslice_overlaps(
                            select s1.id , new_id from genomes_dnaslice s1 where
                                    s1.chromosome_id = new_chromosome_id and
                                    s1.start_pos <= new_end_pos and
                                    s1.end_pos >= new_start_pos and
                                    s1.id <> new_id);

                        /*symetry*/
                        insert into genomes_dnaslice_overlaps(
                            select new_id , s1.id from genomes_dnaslice s1 where
                                    s1.chromosome_id = new_chromosome_id and
                                    s1.start_pos <= new_end_pos and
                                    s1.end_pos >= new_start_pos and
                                    s1.id <> new_id);

                     end;"""

class Migration(migrations.Migration):
    dependencies = [
        ('genomes', '0009_rename_views_and_create_tables')
    ]

    operations = [
        migrations.RunSQL(
            sql = trigger_create,
            reverse_sql = "drop trigger genome_dnasalice_insert_trg;",
        )
    ]