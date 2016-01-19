# -*- coding: utf-8 -*-

def get_or_create_slice(old_target, apps, schema_editor):
    """
    Get or create a slice for a given old target.
    """
    db_alias = schema_editor.connection.alias
    DNASlice = apps.get_model("genomes","DNASlice")
    slice, created = DNASlice.objects.using(db_alias).get_or_create(
        **{k: old_target.__dict__[k] for k in [
            "chromosome_id",  # Now we don't need a select_related.
            "start_pos",
            "end_pos",
    ]})
    return slice
