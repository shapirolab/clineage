# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'TargetEnrichment.chromosome'
        db.add_column(u'linapp_targetenrichment', 'chromosome',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Chromosome'], null=True),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting field 'TargetEnrichment.chromosome'
        db.delete_column(u'linapp_targetenrichment', 'chromosome_id')


    models = {
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Group']", 'symmetrical': 'False', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'linapp.algorithm': {
            'Meta': {'object_name': 'Algorithm'},
            'developers': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.User']", 'symmetrical': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.AlgorithmType']"}),
            'version': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.algorithmparameter': {
            'Meta': {'object_name': 'AlgorithmParameter'},
            'algorithm': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'parameters'", 'symmetrical': 'False', 'to': u"orm['linapp.Algorithm']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'type': ('django.db.models.fields.IntegerField', [], {})
        },
        u'linapp.algorithmrun': {
            'ExtraFiles': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['linapp.ExperimentFile']", 'symmetrical': 'False'}),
            'Meta': {'object_name': 'AlgorithmRun'},
            'algorithm': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'runs'", 'to': u"orm['linapp.Algorithm']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'parameters': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['linapp.AlgorithmParameter']", 'through': u"orm['linapp.AlgorithmRunParameters']", 'symmetrical': 'False'}),
            'runname': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True'}),
            'status': ('django.db.models.fields.IntegerField', [], {}),
            'timestamp': ('django.db.models.fields.DateTimeField', [], {}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"})
        },
        u'linapp.algorithmrunparameters': {
            'Meta': {'object_name': 'AlgorithmRunParameters'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'parameter': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.AlgorithmParameter']"}),
            'run': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.AlgorithmRun']"}),
            'value': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.algorithmtype': {
            'Meta': {'object_name': 'AlgorithmType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'input': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'output': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.assembly': {
            'Meta': {'object_name': 'Assembly'},
            'friendly_name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'taxa': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Taxa']"})
        },
        u'linapp.cell': {
            'Meta': {'object_name': 'Cell'},
            'comment': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'composition': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.SampleComposition']"}),
            'experiment': ('django.db.models.fields.related.ManyToManyField', [], {'blank': 'True', 'related_name': "'cells'", 'null': 'True', 'symmetrical': 'False', 'to': u"orm['linapp.Experiment']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'sampling': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.SamplingEvent']"}),
            'status': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.SampleStatus']", 'null': 'True', 'blank': 'True'})
        },
        u'linapp.cellcontent': {
            'Meta': {'object_name': 'CellContent'},
            'cell': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Cell']"}),
            'comment': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            'panel': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Panel']", 'null': 'True', 'blank': 'True'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.CellContent']", 'null': 'True', 'blank': 'True'}),
            'protocol': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Protocol']", 'null': 'True', 'blank': 'True'}),
            'seq_ready': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.CellContentType']"}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'linapp.cellcontenttype': {
            'Meta': {'object_name': 'CellContentType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.cellselector': {
            'Meta': {'object_name': 'CellSelector', '_ormbases': [u'linapp.SamplingEvent']},
            'coordinates': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Coordinates']", 'null': 'True', 'blank': 'True'}),
            u'samplingevent_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['linapp.SamplingEvent']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'linapp.celltreenode': {
            'Meta': {'object_name': 'CellTreeNode'},
            'cell': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Cell']"}),
            'distance': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            u'level': ('django.db.models.fields.PositiveIntegerField', [], {'db_index': 'True'}),
            u'lft': ('django.db.models.fields.PositiveIntegerField', [], {'db_index': 'True'}),
            'parent': ('mptt.fields.TreeForeignKey', [], {'blank': 'True', 'related_name': "'children'", 'null': 'True', 'to': u"orm['linapp.CellTreeNode']"}),
            u'rght': ('django.db.models.fields.PositiveIntegerField', [], {'db_index': 'True'}),
            u'tree_id': ('django.db.models.fields.PositiveIntegerField', [], {'db_index': 'True'})
        },
        u'linapp.chromosome': {
            'Meta': {'object_name': 'Chromosome'},
            'assembly': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Assembly']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'sequence_length': ('django.db.models.fields.IntegerField', [], {'null': 'True'})
        },
        u'linapp.coordinates': {
            'Meta': {'object_name': 'Coordinates'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'x': ('django.db.models.fields.DecimalField', [], {'max_digits': '10', 'decimal_places': '4'}),
            'y': ('django.db.models.fields.DecimalField', [], {'max_digits': '10', 'decimal_places': '4'}),
            'z': ('django.db.models.fields.DecimalField', [], {'max_digits': '10', 'decimal_places': '4'})
        },
        u'linapp.correctedrawdata': {
            'Meta': {'object_name': 'CorrectedRawData'},
            'creation': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'crd'", 'to': u"orm['linapp.AlgorithmRun']"}),
            'file': ('django.db.models.fields.FilePathField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'rawdata': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.RawData']"})
        },
        u'linapp.dm': {
            'Meta': {'object_name': 'DM'},
            'cell1': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'+'", 'to': u"orm['linapp.GenSig']"}),
            'cell2': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'+'", 'to': u"orm['linapp.GenSig']"}),
            'creation': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'dms'", 'to': u"orm['linapp.AlgorithmRun']"}),
            'distance': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'linapp.experiment': {
            'Meta': {'object_name': 'Experiment'},
            'created_date': ('django.db.models.fields.DateField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_public': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'experiments'", 'symmetrical': 'False', 'through': u"orm['linapp.ExperimentUser']", 'to': u"orm['auth.User']"})
        },
        u'linapp.experimentfile': {
            'Meta': {'object_name': 'ExperimentFile'},
            'context': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.FileContext']"}),
            'description': ('django.db.models.fields.TextField', [], {}),
            'experiment': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Experiment']"}),
            'file': ('django.db.models.fields.files.FileField', [], {'max_length': '100'}),
            'file_name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'upload_date': ('django.db.models.fields.DateField', [], {}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"})
        },
        u'linapp.experimentlog': {
            'Meta': {'object_name': 'ExperimentLog'},
            'comment': ('django.db.models.fields.TextField', [], {}),
            'date': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'experiment': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'comments'", 'to': u"orm['linapp.Experiment']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"})
        },
        u'linapp.experimentuser': {
            'Meta': {'object_name': 'ExperimentUser'},
            'experiment': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Experiment']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'role': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.LineageRole']"}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"})
        },
        u'linapp.extraction': {
            'Meta': {'object_name': 'Extraction'},
            'comment': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'date': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'extraction_event': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.ExtractionEvent']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'organ': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Organ']"}),
            'tissue': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Tissue']"})
        },
        u'linapp.extractionevent': {
            'Meta': {'object_name': 'ExtractionEvent'},
            'comment': ('django.db.models.fields.TextField', [], {}),
            'date': ('django.db.models.fields.DateTimeField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'individual': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Individual']"}),
            'location': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Location']", 'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'user_documented': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'+'", 'to': u"orm['auth.User']"}),
            'user_performed': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'+'", 'to': u"orm['auth.User']"})
        },
        u'linapp.facs': {
            'Meta': {'object_name': 'FACS', '_ormbases': [u'linapp.SamplingEvent']},
            'marker': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.FACSMarker']"}),
            u'samplingevent_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['linapp.SamplingEvent']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'linapp.facsmarker': {
            'Meta': {'object_name': 'FACSMarker'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.failedtargetvalue': {
            'Meta': {'object_name': 'FailedTargetValue'},
            'comment': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'linapp.filecontext': {
            'Meta': {'object_name': 'FileContext'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.geneticbackground': {
            'Meta': {'object_name': 'GeneticBackground'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.gensig': {
            'Meta': {'object_name': 'GenSig'},
            'creation': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'gensigs'", 'to': u"orm['linapp.AlgorithmRun']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {}),
            'variants': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['linapp.TargetVariant']", 'symmetrical': 'False'})
        },
        u'linapp.individual': {
            'Meta': {'object_name': 'Individual'},
            'background': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.GeneticBackground']", 'null': 'True', 'blank': 'True'}),
            'born': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'comment': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'location': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Location']", 'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'sacrificed': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'sex': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'taxa': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Taxa']"})
        },
        u'linapp.lasercapture': {
            'Meta': {'object_name': 'LaserCapture', '_ormbases': [u'linapp.SamplingEvent']},
            'coordinates': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Coordinates']", 'null': 'True', 'blank': 'True'}),
            u'samplingevent_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['linapp.SamplingEvent']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'linapp.lineagerole': {
            'Meta': {'object_name': 'LineageRole'},
            'delete': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'read': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'write': ('django.db.models.fields.BooleanField', [], {'default': 'False'})
        },
        u'linapp.location': {
            'Meta': {'object_name': 'Location'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.machine': {
            'Meta': {'object_name': 'Machine'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'ip': ('django.db.models.fields.IPAddressField', [], {'max_length': '15', 'null': 'True', 'blank': 'True'}),
            'machineid': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.MachineType']"})
        },
        u'linapp.machinetype': {
            'Meta': {'object_name': 'MachineType'},
            'company': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.microsatellite': {
            'Meta': {'object_name': 'Microsatellite', '_ormbases': [u'linapp.Target']},
            'repeat_number': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'repeat_type': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'repeat_unit': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            u'target_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['linapp.Target']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'linapp.organ': {
            'Meta': {'object_name': 'Organ'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.panel': {
            'Meta': {'object_name': 'Panel'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'targets': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'panels'", 'symmetrical': 'False', 'to': u"orm['linapp.TargetEnrichment']"})
        },
        u'linapp.plate': {
            'Meta': {'object_name': 'Plate'},
            'barcode': ('django.db.models.fields.CharField', [], {'max_length': '20', 'blank': 'True'}),
            'code': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'lastusedwell': ('django.db.models.fields.CharField', [], {'default': "'A1'", 'max_length': '4'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'state': ('django.db.models.fields.CharField', [], {'max_length': '20', 'blank': 'True'}),
            'timestamp': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.PlateType']"})
        },
        u'linapp.platecontext': {
            'Meta': {'object_name': 'PlateContext'},
            'description': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'linapp.plateplastica': {
            'Meta': {'object_name': 'PlatePlastica'},
            'code': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'columns': ('django.db.models.fields.IntegerField', [], {'default': '12'}),
            'description': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'rows': ('django.db.models.fields.IntegerField', [], {'default': '8'})
        },
        u'linapp.platestorage': {
            'Meta': {'object_name': 'PlateStorage'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'inner_location': ('django.db.models.fields.CharField', [], {'max_length': '100', 'blank': 'True'}),
            'notes': ('django.db.models.fields.CharField', [], {'max_length': '250', 'blank': 'True'}),
            'plate': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Plate']"}),
            'storageBox': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.StorageBox']"})
        },
        u'linapp.platetype': {
            'Meta': {'object_name': 'PlateType'},
            'code': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'context': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.PlateContext']", 'null': 'True'}),
            'friendly': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'plastic': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.PlatePlastica']", 'null': 'True'})
        },
        u'linapp.primer': {
            'Meta': {'object_name': 'Primer', '_ormbases': [u'linapp.Target']},
            'sequence': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Sequence']"}),
            u'target_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['linapp.Target']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'linapp.primersmultiplex': {
            'Meta': {'object_name': 'PrimersMultiplex'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'primers': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['linapp.TargetEnrichment']", 'symmetrical': 'False'})
        },
        u'linapp.protocol': {
            'Meta': {'object_name': 'Protocol'},
            'abstract': ('django.db.models.fields.TextField', [], {}),
            'file': ('django.db.models.fields.FilePathField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'fulldescription': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'initials': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'kit': ('django.db.models.fields.CharField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.ProtocolType']"})
        },
        u'linapp.protocoltype': {
            'Meta': {'object_name': 'ProtocolType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'linapp.rawdata': {
            'Meta': {'object_name': 'RawData'},
            'file': ('django.db.models.fields.FilePathField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sequencing': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Sequencing']"}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"})
        },
        u'linapp.samplecomposition': {
            'Meta': {'object_name': 'SampleComposition'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.samplelocation': {
            'Meta': {'object_name': 'SampleLocation'},
            'concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '5', 'blank': 'True'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'object_id': ('django.db.models.fields.PositiveIntegerField', [], {}),
            'plate': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Plate']"}),
            'timestamp': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'volume': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '3', 'blank': 'True'}),
            'well': ('django.db.models.fields.CharField', [], {'max_length': '3', 'blank': 'True'})
        },
        u'linapp.samplestatus': {
            'Meta': {'object_name': 'SampleStatus'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.samplingevent': {
            'Meta': {'object_name': 'SamplingEvent'},
            'attachment': ('django.db.models.fields.files.FileField', [], {'max_length': '100', 'null': 'True', 'blank': 'True'}),
            'comment': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'extraction': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Extraction']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True', 'blank': 'True'})
        },
        u'linapp.sequence': {
            'Meta': {'object_name': 'Sequence'},
            'hash': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '32'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.IntegerField', [], {}),
            'sequence': ('django.db.models.fields.TextField', [], {})
        },
        u'linapp.sequencedistribution': {
            'Meta': {'object_name': 'SequenceDistribution'},
            'count': ('django.db.models.fields.IntegerField', [], {}),
            'failed': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.FailedTargetValue']", 'null': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'target': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['linapp.Target']", 'through': u"orm['linapp.TargetAnalysis']", 'symmetrical': 'False'}),
            'value': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Sequence']"})
        },
        u'linapp.sequencing': {
            'Meta': {'object_name': 'Sequencing'},
            'data': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'sequencing_event'", 'null': 'True', 'to': u"orm['linapp.RawData']"}),
            'date': ('django.db.models.fields.DateField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'machine': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Machine']"}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '100'}),
            'protocol': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Protocol']"}),
            'samples': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['linapp.CellContent']", 'symmetrical': 'False'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"})
        },
        u'linapp.storagebox': {
            'Meta': {'object_name': 'StorageBox'},
            'barcode': ('django.db.models.fields.CharField', [], {'max_length': '20', 'blank': 'True'}),
            'code': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100', 'blank': 'True'}),
            'storage_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.StorageType']"})
        },
        u'linapp.storagetype': {
            'Meta': {'object_name': 'StorageType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100', 'blank': 'True'}),
            'temperature': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '1', 'blank': 'True'})
        },
        u'linapp.target': {
            'Meta': {'object_name': 'Target'},
            'assembly': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Assembly']"}),
            'chromosome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Chromosome']"}),
            'end_pos': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'referencevalue': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Sequence']"}),
            'start_pos': ('django.db.models.fields.IntegerField', [], {}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.TargetType']"})
        },
        u'linapp.targetanalysis': {
            'Meta': {'object_name': 'TargetAnalysis'},
            'creation': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'targetsdistributions'", 'to': u"orm['linapp.AlgorithmRun']"}),
            'distribution': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.SequenceDistribution']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sequencingdata': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.CorrectedRawData']"}),
            'target': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Target']"})
        },
        u'linapp.targetenrichment': {
            'Meta': {'object_name': 'TargetEnrichment'},
            'chromosome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Chromosome']", 'null': 'True'}),
            'comment': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'left': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'left_primer'", 'to': u"orm['linapp.Primer']"}),
            'passed_validation': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'}),
            'right': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'right_primer'", 'to': u"orm['linapp.Primer']"}),
            'targets': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['linapp.Target']", 'null': 'True', 'blank': 'True'}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.TargetEnrichmentType']"}),
            'validation_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'validation_failure': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.TargetEnrichmentFailureType']", 'null': 'True'})
        },
        u'linapp.targetenrichmentfailuretype': {
            'Meta': {'object_name': 'TargetEnrichmentFailureType'},
            'description': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.targetenrichmenttype': {
            'Meta': {'object_name': 'TargetEnrichmentType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'protocol': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Protocol']", 'null': 'True'})
        },
        u'linapp.targettype': {
            'Meta': {'object_name': 'TargetType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.targetvariant': {
            'Meta': {'object_name': 'TargetVariant'},
            'creation': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'targetsvariants'", 'to': u"orm['linapp.AlgorithmRun']"}),
            'distribution': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['linapp.SequenceDistribution']", 'symmetrical': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.TargetVariantType']"}),
            'value': ('django.db.models.fields.PositiveIntegerField', [], {})
        },
        u'linapp.targetvarianttype': {
            'Meta': {'object_name': 'TargetVariantType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.taxa': {
            'Meta': {'object_name': 'Taxa'},
            'friendly_name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'parent': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'rank': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'taxonomy_id': ('django.db.models.fields.IntegerField', [], {})
        },
        u'linapp.tissue': {
            'Meta': {'object_name': 'Tissue'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'linapp.userprofile': {
            'Meta': {'object_name': 'UserProfile'},
            'comment': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'institute': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['auth.User']", 'unique': 'True'})
        }
    }

    complete_apps = ['linapp']