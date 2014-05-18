# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'LineageRole'
        db.create_table(u'linapp_lineagerole', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('read', self.gf('django.db.models.fields.BooleanField')(default=True)),
            ('write', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('delete', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'linapp', ['LineageRole'])

        # Adding model 'UserProfile'
        db.create_table(u'linapp_userprofile', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['auth.User'], unique=True)),
            ('institute', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('comment', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'linapp', ['UserProfile'])

        # Adding model 'Experiment'
        db.create_table(u'linapp_experiment', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('created_date', self.gf('django.db.models.fields.DateField')(auto_now_add=True, blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')()),
            ('is_public', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'linapp', ['Experiment'])

        # Adding model 'ExperimentLog'
        db.create_table(u'linapp_experimentlog', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('experiment', self.gf('django.db.models.fields.related.ForeignKey')(related_name='comments', to=orm['linapp.Experiment'])),
            ('date', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
            ('comment', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'linapp', ['ExperimentLog'])

        # Adding model 'FileContext'
        db.create_table(u'linapp_filecontext', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['FileContext'])

        # Adding model 'ExperimentFile'
        db.create_table(u'linapp_experimentfile', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('title', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('experiment', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Experiment'])),
            ('file_name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('file', self.gf('django.db.models.fields.files.FileField')(max_length=100)),
            ('context', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.FileContext'])),
            ('upload_date', self.gf('django.db.models.fields.DateField')()),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('description', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'linapp', ['ExperimentFile'])

        # Adding model 'ExperimentUser'
        db.create_table(u'linapp_experimentuser', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('experiment', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Experiment'])),
            ('role', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.LineageRole'])),
        ))
        db.send_create_signal(u'linapp', ['ExperimentUser'])

        # Adding model 'ProtocolType'
        db.create_table(u'linapp_protocoltype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100)),
        ))
        db.send_create_signal(u'linapp', ['ProtocolType'])

        # Adding model 'SampleComposition'
        db.create_table(u'linapp_samplecomposition', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['SampleComposition'])

        # Adding model 'CellContentType'
        db.create_table(u'linapp_cellcontenttype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['CellContentType'])

        # Adding model 'TargetType'
        db.create_table(u'linapp_targettype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['TargetType'])

        # Adding model 'TargetVariantType'
        db.create_table(u'linapp_targetvarianttype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['TargetVariantType'])

        # Adding model 'SampleStatus'
        db.create_table(u'linapp_samplestatus', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['SampleStatus'])

        # Adding model 'TargetEnrichmentFailureType'
        db.create_table(u'linapp_targetenrichmentfailuretype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('description', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['TargetEnrichmentFailureType'])

        # Adding model 'Sequence'
        db.create_table(u'linapp_sequence', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('length', self.gf('django.db.models.fields.IntegerField')()),
            ('sequence', self.gf('django.db.models.fields.TextField')()),
            ('hash', self.gf('django.db.models.fields.CharField')(unique=True, max_length=32)),
        ))
        db.send_create_signal(u'linapp', ['Sequence'])

        # Adding model 'Taxa'
        db.create_table(u'linapp_taxa', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('taxonomy_id', self.gf('django.db.models.fields.IntegerField')()),
            ('rank', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('parent', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('friendly_name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['Taxa'])

        # Adding model 'GeneticBackground'
        db.create_table(u'linapp_geneticbackground', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['GeneticBackground'])

        # Adding model 'Organ'
        db.create_table(u'linapp_organ', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['Organ'])

        # Adding model 'Tissue'
        db.create_table(u'linapp_tissue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['Tissue'])

        # Adding model 'Assembly'
        db.create_table(u'linapp_assembly', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('taxa', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Taxa'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('friendly_name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['Assembly'])

        # Adding model 'Chromosome'
        db.create_table(u'linapp_chromosome', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('assembly', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Assembly'])),
            ('sequence_length', self.gf('django.db.models.fields.IntegerField')(null=True)),
        ))
        db.send_create_signal(u'linapp', ['Chromosome'])

        # Adding model 'Protocol'
        db.create_table(u'linapp_protocol', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('initials', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('abstract', self.gf('django.db.models.fields.TextField')()),
            ('fulldescription', self.gf('django.db.models.fields.TextField')()),
            ('kit', self.gf('django.db.models.fields.CharField')(max_length=100, null=True, blank=True)),
            ('type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.ProtocolType'])),
            ('file', self.gf('django.db.models.fields.FilePathField')(max_length=100, null=True, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['Protocol'])

        # Adding model 'Target'
        db.create_table(u'linapp_target', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.TargetType'])),
            ('chromosome', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Chromosome'])),
            ('start_pos', self.gf('django.db.models.fields.IntegerField')()),
            ('end_pos', self.gf('django.db.models.fields.IntegerField')()),
            ('referencevalue', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Sequence'])),
        ))
        db.send_create_signal(u'linapp', ['Target'])

        # Adding model 'PrimerTail'
        db.create_table(u'linapp_primertail', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('tail', self.gf('django.db.models.fields.CharField')(max_length=50, null=True)),
        ))
        db.send_create_signal(u'linapp', ['PrimerTail'])

        # Adding model 'Primer'
        db.create_table(u'linapp_primer', (
            (u'target_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['linapp.Target'], unique=True, primary_key=True)),
            ('strand', self.gf('django.db.models.fields.CharField')(max_length=1, null=True)),
            ('sequence', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Sequence'])),
            ('tail', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.PrimerTail'], null=True)),
        ))
        db.send_create_signal(u'linapp', ['Primer'])

        # Adding model 'Microsatellite'
        db.create_table(u'linapp_microsatellite', (
            (u'target_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['linapp.Target'], unique=True, primary_key=True)),
            ('repeat_type', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('repeat_unit', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('repeat_number', self.gf('django.db.models.fields.PositiveIntegerField')()),
        ))
        db.send_create_signal(u'linapp', ['Microsatellite'])

        # Adding model 'TargetEnrichmentType'
        db.create_table(u'linapp_targetenrichmenttype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('protocol', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Protocol'], null=True)),
        ))
        db.send_create_signal(u'linapp', ['TargetEnrichmentType'])

        # Adding model 'TargetEnrichment'
        db.create_table(u'linapp_targetenrichment', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.TargetEnrichmentType'])),
            ('chromosome', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Chromosome'])),
            ('left', self.gf('django.db.models.fields.related.ForeignKey')(related_name='left_primer', to=orm['linapp.Primer'])),
            ('right', self.gf('django.db.models.fields.related.ForeignKey')(related_name='right_primer', to=orm['linapp.Primer'])),
            ('amplicon', self.gf('django.db.models.fields.CharField')(max_length=500)),
            ('passed_validation', self.gf('django.db.models.fields.NullBooleanField')(null=True, blank=True)),
            ('validation_failure', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.TargetEnrichmentFailureType'], null=True)),
            ('validation_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('comment', self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['TargetEnrichment'])

        # Adding M2M table for field targets on 'TargetEnrichment'
        m2m_table_name = db.shorten_name(u'linapp_targetenrichment_targets')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('targetenrichment', models.ForeignKey(orm[u'linapp.targetenrichment'], null=False)),
            ('target', models.ForeignKey(orm[u'linapp.target'], null=False))
        ))
        db.create_unique(m2m_table_name, ['targetenrichment_id', 'target_id'])

        # Adding model 'Coordinates'
        db.create_table(u'linapp_coordinates', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('x', self.gf('django.db.models.fields.DecimalField')(max_digits=10, decimal_places=4)),
            ('y', self.gf('django.db.models.fields.DecimalField')(max_digits=10, decimal_places=4)),
            ('z', self.gf('django.db.models.fields.DecimalField')(max_digits=10, decimal_places=4)),
        ))
        db.send_create_signal(u'linapp', ['Coordinates'])

        # Adding model 'FACSMarker'
        db.create_table(u'linapp_facsmarker', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['FACSMarker'])

        # Adding model 'Panel'
        db.create_table(u'linapp_panel', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['Panel'])

        # Adding M2M table for field targets on 'Panel'
        m2m_table_name = db.shorten_name(u'linapp_panel_targets')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('panel', models.ForeignKey(orm[u'linapp.panel'], null=False)),
            ('targetenrichment', models.ForeignKey(orm[u'linapp.targetenrichment'], null=False))
        ))
        db.create_unique(m2m_table_name, ['panel_id', 'targetenrichment_id'])

        # Adding model 'Location'
        db.create_table(u'linapp_location', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['Location'])

        # Adding model 'AlgorithmType'
        db.create_table(u'linapp_algorithmtype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('input', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('output', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['AlgorithmType'])

        # Adding model 'Algorithm'
        db.create_table(u'linapp_algorithm', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.AlgorithmType'])),
            ('version', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['Algorithm'])

        # Adding M2M table for field developers on 'Algorithm'
        m2m_table_name = db.shorten_name(u'linapp_algorithm_developers')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('algorithm', models.ForeignKey(orm[u'linapp.algorithm'], null=False)),
            ('user', models.ForeignKey(orm[u'auth.user'], null=False))
        ))
        db.create_unique(m2m_table_name, ['algorithm_id', 'user_id'])

        # Adding model 'AlgorithmParameter'
        db.create_table(u'linapp_algorithmparameter', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('type', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'linapp', ['AlgorithmParameter'])

        # Adding M2M table for field algorithm on 'AlgorithmParameter'
        m2m_table_name = db.shorten_name(u'linapp_algorithmparameter_algorithm')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('algorithmparameter', models.ForeignKey(orm[u'linapp.algorithmparameter'], null=False)),
            ('algorithm', models.ForeignKey(orm[u'linapp.algorithm'], null=False))
        ))
        db.create_unique(m2m_table_name, ['algorithmparameter_id', 'algorithm_id'])

        # Adding model 'AlgorithmRun'
        db.create_table(u'linapp_algorithmrun', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('algorithm', self.gf('django.db.models.fields.related.ForeignKey')(related_name='runs', to=orm['linapp.Algorithm'])),
            ('runname', self.gf('django.db.models.fields.CharField')(max_length=50, null=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('timestamp', self.gf('django.db.models.fields.DateTimeField')()),
            ('status', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'linapp', ['AlgorithmRun'])

        # Adding M2M table for field ExtraFiles on 'AlgorithmRun'
        m2m_table_name = db.shorten_name(u'linapp_algorithmrun_ExtraFiles')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('algorithmrun', models.ForeignKey(orm[u'linapp.algorithmrun'], null=False)),
            ('experimentfile', models.ForeignKey(orm[u'linapp.experimentfile'], null=False))
        ))
        db.create_unique(m2m_table_name, ['algorithmrun_id', 'experimentfile_id'])

        # Adding model 'AlgorithmRunParameters'
        db.create_table(u'linapp_algorithmrunparameters', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('run', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.AlgorithmRun'])),
            ('parameter', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.AlgorithmParameter'])),
            ('value', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['AlgorithmRunParameters'])

        # Adding model 'Individual'
        db.create_table(u'linapp_individual', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('taxa', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Taxa'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('sex', self.gf('django.db.models.fields.CharField')(max_length=1)),
            ('born', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('comment', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('background', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.GeneticBackground'], null=True, blank=True)),
            ('location', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Location'], null=True, blank=True)),
            ('sacrificed', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['Individual'])

        # Adding model 'ExtractionEvent'
        db.create_table(u'linapp_extractionevent', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('individual', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Individual'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('comment', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('date', self.gf('django.db.models.fields.DateTimeField')()),
            ('location', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Location'], null=True, blank=True)),
            ('user_performed', self.gf('django.db.models.fields.related.ForeignKey')(related_name='+', to=orm['auth.User'])),
            ('user_documented', self.gf('django.db.models.fields.related.ForeignKey')(related_name='+', to=orm['auth.User'])),
        ))
        db.send_create_signal(u'linapp', ['ExtractionEvent'])

        # Adding model 'Extraction'
        db.create_table(u'linapp_extraction', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('extraction_event', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.ExtractionEvent'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('date', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('organ', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Organ'])),
            ('tissue', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Tissue'])),
            ('comment', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['Extraction'])

        # Adding model 'SamplingEvent'
        db.create_table(u'linapp_samplingevent', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('extraction', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Extraction'])),
            ('date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'], null=True, blank=True)),
            ('comment', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('attachment', self.gf('django.db.models.fields.files.FileField')(max_length=100, null=True, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['SamplingEvent'])

        # Adding model 'FACS'
        db.create_table(u'linapp_facs', (
            (u'samplingevent_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['linapp.SamplingEvent'], unique=True, primary_key=True)),
            ('marker', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.FACSMarker'])),
        ))
        db.send_create_signal(u'linapp', ['FACS'])

        # Adding model 'LaserCapture'
        db.create_table(u'linapp_lasercapture', (
            (u'samplingevent_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['linapp.SamplingEvent'], unique=True, primary_key=True)),
            ('coordinates', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Coordinates'], null=True, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['LaserCapture'])

        # Adding model 'CellSelector'
        db.create_table(u'linapp_cellselector', (
            (u'samplingevent_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['linapp.SamplingEvent'], unique=True, primary_key=True)),
            ('coordinates', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Coordinates'], null=True, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['CellSelector'])

        # Adding model 'Cell'
        db.create_table(u'linapp_cell', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sampling', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.SamplingEvent'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('composition', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.SampleComposition'])),
            ('status', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.SampleStatus'], null=True, blank=True)),
            ('comment', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['Cell'])

        # Adding M2M table for field experiment on 'Cell'
        m2m_table_name = db.shorten_name(u'linapp_cell_experiment')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('cell', models.ForeignKey(orm[u'linapp.cell'], null=False)),
            ('experiment', models.ForeignKey(orm[u'linapp.experiment'], null=False))
        ))
        db.create_unique(m2m_table_name, ['cell_id', 'experiment_id'])

        # Adding model 'CellContent'
        db.create_table(u'linapp_cellcontent', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('parent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.CellContent'], null=True, blank=True)),
            ('cell', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Cell'])),
            ('panel', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Panel'], null=True, blank=True)),
            ('type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.CellContentType'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50, null=True, blank=True)),
            ('protocol', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Protocol'], null=True, blank=True)),
            ('seq_ready', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'], null=True, blank=True)),
            ('comment', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'linapp', ['CellContent'])

        # Adding model 'MachineType'
        db.create_table(u'linapp_machinetype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('company', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('model', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'linapp', ['MachineType'])

        # Adding model 'Machine'
        db.create_table(u'linapp_machine', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('machineid', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100, null=True, blank=True)),
            ('type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.MachineType'])),
            ('ip', self.gf('django.db.models.fields.IPAddressField')(max_length=15, null=True, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['Machine'])

        # Adding model 'Sequencing'
        db.create_table(u'linapp_sequencing', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('data', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='sequencing_event', null=True, to=orm['linapp.RawData'])),
            ('name', self.gf('django.db.models.fields.CharField')(unique=True, max_length=100)),
            ('machine', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Machine'])),
            ('protocol', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Protocol'])),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
            ('date', self.gf('django.db.models.fields.DateField')()),
        ))
        db.send_create_signal(u'linapp', ['Sequencing'])

        # Adding M2M table for field samples on 'Sequencing'
        m2m_table_name = db.shorten_name(u'linapp_sequencing_samples')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('sequencing', models.ForeignKey(orm[u'linapp.sequencing'], null=False)),
            ('cellcontent', models.ForeignKey(orm[u'linapp.cellcontent'], null=False))
        ))
        db.create_unique(m2m_table_name, ['sequencing_id', 'cellcontent_id'])

        # Adding model 'RawData'
        db.create_table(u'linapp_rawdata', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sequencing', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Sequencing'])),
            ('file', self.gf('django.db.models.fields.FilePathField')(max_length=100, null=True, blank=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'])),
        ))
        db.send_create_signal(u'linapp', ['RawData'])

        # Adding model 'CorrectedRawData'
        db.create_table(u'linapp_correctedrawdata', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('rawdata', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.RawData'])),
            ('file', self.gf('django.db.models.fields.FilePathField')(max_length=100)),
            ('creation', self.gf('django.db.models.fields.related.ForeignKey')(related_name='crd', to=orm['linapp.AlgorithmRun'])),
        ))
        db.send_create_signal(u'linapp', ['CorrectedRawData'])

        # Adding model 'FailedTargetValue'
        db.create_table(u'linapp_failedtargetvalue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('comment', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'linapp', ['FailedTargetValue'])

        # Adding model 'SequenceDistribution'
        db.create_table(u'linapp_sequencedistribution', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Sequence'])),
            ('count', self.gf('django.db.models.fields.IntegerField')()),
            ('failed', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.FailedTargetValue'], null=True)),
        ))
        db.send_create_signal(u'linapp', ['SequenceDistribution'])

        # Adding model 'TargetAnalysis'
        db.create_table(u'linapp_targetanalysis', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('target', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Target'])),
            ('distribution', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.SequenceDistribution'])),
            ('creation', self.gf('django.db.models.fields.related.ForeignKey')(related_name='targetsdistributions', to=orm['linapp.AlgorithmRun'])),
            ('sequencingdata', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.CorrectedRawData'])),
        ))
        db.send_create_signal(u'linapp', ['TargetAnalysis'])

        # Adding model 'TargetVariant'
        db.create_table(u'linapp_targetvariant', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('creation', self.gf('django.db.models.fields.related.ForeignKey')(related_name='targetsvariants', to=orm['linapp.AlgorithmRun'])),
            ('type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.TargetVariantType'])),
            ('value', self.gf('django.db.models.fields.PositiveIntegerField')()),
        ))
        db.send_create_signal(u'linapp', ['TargetVariant'])

        # Adding M2M table for field distribution on 'TargetVariant'
        m2m_table_name = db.shorten_name(u'linapp_targetvariant_distribution')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('targetvariant', models.ForeignKey(orm[u'linapp.targetvariant'], null=False)),
            ('sequencedistribution', models.ForeignKey(orm[u'linapp.sequencedistribution'], null=False))
        ))
        db.create_unique(m2m_table_name, ['targetvariant_id', 'sequencedistribution_id'])

        # Adding model 'GenSig'
        db.create_table(u'linapp_gensig', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('creation', self.gf('django.db.models.fields.related.ForeignKey')(related_name='gensigs', to=orm['linapp.AlgorithmRun'])),
            ('value', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'linapp', ['GenSig'])

        # Adding M2M table for field variants on 'GenSig'
        m2m_table_name = db.shorten_name(u'linapp_gensig_variants')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('gensig', models.ForeignKey(orm[u'linapp.gensig'], null=False)),
            ('targetvariant', models.ForeignKey(orm[u'linapp.targetvariant'], null=False))
        ))
        db.create_unique(m2m_table_name, ['gensig_id', 'targetvariant_id'])

        # Adding model 'DM'
        db.create_table(u'linapp_dm', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('cell1', self.gf('django.db.models.fields.related.ForeignKey')(related_name='+', to=orm['linapp.GenSig'])),
            ('cell2', self.gf('django.db.models.fields.related.ForeignKey')(related_name='+', to=orm['linapp.GenSig'])),
            ('distance', self.gf('django.db.models.fields.IntegerField')()),
            ('creation', self.gf('django.db.models.fields.related.ForeignKey')(related_name='dms', to=orm['linapp.AlgorithmRun'])),
        ))
        db.send_create_signal(u'linapp', ['DM'])

        # Adding model 'CellTreeNode'
        db.create_table(u'linapp_celltreenode', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('cell', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Cell'])),
            ('parent', self.gf('mptt.fields.TreeForeignKey')(blank=True, related_name='children', null=True, to=orm['linapp.CellTreeNode'])),
            ('distance', self.gf('django.db.models.fields.IntegerField')()),
            (u'lft', self.gf('django.db.models.fields.PositiveIntegerField')(db_index=True)),
            (u'rght', self.gf('django.db.models.fields.PositiveIntegerField')(db_index=True)),
            (u'tree_id', self.gf('django.db.models.fields.PositiveIntegerField')(db_index=True)),
            (u'level', self.gf('django.db.models.fields.PositiveIntegerField')(db_index=True)),
        ))
        db.send_create_signal(u'linapp', ['CellTreeNode'])

        # Adding model 'StorageType'
        db.create_table(u'linapp_storagetype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100, blank=True)),
            ('temperature', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=5, decimal_places=1, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['StorageType'])

        # Adding model 'StorageBox'
        db.create_table(u'linapp_storagebox', (
            ('code', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('storage_type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.StorageType'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100, blank=True)),
            ('barcode', self.gf('django.db.models.fields.CharField')(max_length=20, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['StorageBox'])

        # Adding model 'PlateContext'
        db.create_table(u'linapp_platecontext', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=30, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['PlateContext'])

        # Adding model 'PlatePlastica'
        db.create_table(u'linapp_plateplastica', (
            ('code', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('description', self.gf('django.db.models.fields.CharField')(max_length=30, blank=True)),
            ('rows', self.gf('django.db.models.fields.IntegerField')(default=8)),
            ('columns', self.gf('django.db.models.fields.IntegerField')(default=12)),
        ))
        db.send_create_signal(u'linapp', ['PlatePlastica'])

        # Adding model 'PlateType'
        db.create_table(u'linapp_platetype', (
            ('code', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('friendly', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('context', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.PlateContext'], null=True)),
            ('plastic', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.PlatePlastica'], null=True)),
        ))
        db.send_create_signal(u'linapp', ['PlateType'])

        # Adding model 'Plate'
        db.create_table(u'linapp_plate', (
            ('code', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.PlateType'])),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('barcode', self.gf('django.db.models.fields.CharField')(max_length=20, blank=True)),
            ('timestamp', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('state', self.gf('django.db.models.fields.CharField')(max_length=20, blank=True)),
            ('lastusedwell', self.gf('django.db.models.fields.CharField')(default='A1', max_length=4)),
        ))
        db.send_create_signal(u'linapp', ['Plate'])

        # Adding model 'PlateStorage'
        db.create_table(u'linapp_platestorage', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('storageBox', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.StorageBox'])),
            ('plate', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Plate'])),
            ('inner_location', self.gf('django.db.models.fields.CharField')(max_length=100, blank=True)),
            ('notes', self.gf('django.db.models.fields.CharField')(max_length=250, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['PlateStorage'])

        # Adding model 'SampleLocation'
        db.create_table(u'linapp_samplelocation', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('plate', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['linapp.Plate'])),
            ('well', self.gf('django.db.models.fields.CharField')(max_length=3, blank=True)),
            ('content_type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['contenttypes.ContentType'])),
            ('object_id', self.gf('django.db.models.fields.PositiveIntegerField')()),
            ('volume', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=10, decimal_places=3, blank=True)),
            ('concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=10, decimal_places=5, blank=True)),
            ('timestamp', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
        ))
        db.send_create_signal(u'linapp', ['SampleLocation'])

        # Adding model 'PrimersMultiplex'
        db.create_table(u'linapp_primersmultiplex', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=20)),
        ))
        db.send_create_signal(u'linapp', ['PrimersMultiplex'])

        # Adding M2M table for field primers on 'PrimersMultiplex'
        m2m_table_name = db.shorten_name(u'linapp_primersmultiplex_primers')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('primersmultiplex', models.ForeignKey(orm[u'linapp.primersmultiplex'], null=False)),
            ('targetenrichment', models.ForeignKey(orm[u'linapp.targetenrichment'], null=False))
        ))
        db.create_unique(m2m_table_name, ['primersmultiplex_id', 'targetenrichment_id'])


    def backwards(self, orm):
        # Deleting model 'LineageRole'
        db.delete_table(u'linapp_lineagerole')

        # Deleting model 'UserProfile'
        db.delete_table(u'linapp_userprofile')

        # Deleting model 'Experiment'
        db.delete_table(u'linapp_experiment')

        # Deleting model 'ExperimentLog'
        db.delete_table(u'linapp_experimentlog')

        # Deleting model 'FileContext'
        db.delete_table(u'linapp_filecontext')

        # Deleting model 'ExperimentFile'
        db.delete_table(u'linapp_experimentfile')

        # Deleting model 'ExperimentUser'
        db.delete_table(u'linapp_experimentuser')

        # Deleting model 'ProtocolType'
        db.delete_table(u'linapp_protocoltype')

        # Deleting model 'SampleComposition'
        db.delete_table(u'linapp_samplecomposition')

        # Deleting model 'CellContentType'
        db.delete_table(u'linapp_cellcontenttype')

        # Deleting model 'TargetType'
        db.delete_table(u'linapp_targettype')

        # Deleting model 'TargetVariantType'
        db.delete_table(u'linapp_targetvarianttype')

        # Deleting model 'SampleStatus'
        db.delete_table(u'linapp_samplestatus')

        # Deleting model 'TargetEnrichmentFailureType'
        db.delete_table(u'linapp_targetenrichmentfailuretype')

        # Deleting model 'Sequence'
        db.delete_table(u'linapp_sequence')

        # Deleting model 'Taxa'
        db.delete_table(u'linapp_taxa')

        # Deleting model 'GeneticBackground'
        db.delete_table(u'linapp_geneticbackground')

        # Deleting model 'Organ'
        db.delete_table(u'linapp_organ')

        # Deleting model 'Tissue'
        db.delete_table(u'linapp_tissue')

        # Deleting model 'Assembly'
        db.delete_table(u'linapp_assembly')

        # Deleting model 'Chromosome'
        db.delete_table(u'linapp_chromosome')

        # Deleting model 'Protocol'
        db.delete_table(u'linapp_protocol')

        # Deleting model 'Target'
        db.delete_table(u'linapp_target')

        # Deleting model 'PrimerTail'
        db.delete_table(u'linapp_primertail')

        # Deleting model 'Primer'
        db.delete_table(u'linapp_primer')

        # Deleting model 'Microsatellite'
        db.delete_table(u'linapp_microsatellite')

        # Deleting model 'TargetEnrichmentType'
        db.delete_table(u'linapp_targetenrichmenttype')

        # Deleting model 'TargetEnrichment'
        db.delete_table(u'linapp_targetenrichment')

        # Removing M2M table for field targets on 'TargetEnrichment'
        db.delete_table(db.shorten_name(u'linapp_targetenrichment_targets'))

        # Deleting model 'Coordinates'
        db.delete_table(u'linapp_coordinates')

        # Deleting model 'FACSMarker'
        db.delete_table(u'linapp_facsmarker')

        # Deleting model 'Panel'
        db.delete_table(u'linapp_panel')

        # Removing M2M table for field targets on 'Panel'
        db.delete_table(db.shorten_name(u'linapp_panel_targets'))

        # Deleting model 'Location'
        db.delete_table(u'linapp_location')

        # Deleting model 'AlgorithmType'
        db.delete_table(u'linapp_algorithmtype')

        # Deleting model 'Algorithm'
        db.delete_table(u'linapp_algorithm')

        # Removing M2M table for field developers on 'Algorithm'
        db.delete_table(db.shorten_name(u'linapp_algorithm_developers'))

        # Deleting model 'AlgorithmParameter'
        db.delete_table(u'linapp_algorithmparameter')

        # Removing M2M table for field algorithm on 'AlgorithmParameter'
        db.delete_table(db.shorten_name(u'linapp_algorithmparameter_algorithm'))

        # Deleting model 'AlgorithmRun'
        db.delete_table(u'linapp_algorithmrun')

        # Removing M2M table for field ExtraFiles on 'AlgorithmRun'
        db.delete_table(db.shorten_name(u'linapp_algorithmrun_ExtraFiles'))

        # Deleting model 'AlgorithmRunParameters'
        db.delete_table(u'linapp_algorithmrunparameters')

        # Deleting model 'Individual'
        db.delete_table(u'linapp_individual')

        # Deleting model 'ExtractionEvent'
        db.delete_table(u'linapp_extractionevent')

        # Deleting model 'Extraction'
        db.delete_table(u'linapp_extraction')

        # Deleting model 'SamplingEvent'
        db.delete_table(u'linapp_samplingevent')

        # Deleting model 'FACS'
        db.delete_table(u'linapp_facs')

        # Deleting model 'LaserCapture'
        db.delete_table(u'linapp_lasercapture')

        # Deleting model 'CellSelector'
        db.delete_table(u'linapp_cellselector')

        # Deleting model 'Cell'
        db.delete_table(u'linapp_cell')

        # Removing M2M table for field experiment on 'Cell'
        db.delete_table(db.shorten_name(u'linapp_cell_experiment'))

        # Deleting model 'CellContent'
        db.delete_table(u'linapp_cellcontent')

        # Deleting model 'MachineType'
        db.delete_table(u'linapp_machinetype')

        # Deleting model 'Machine'
        db.delete_table(u'linapp_machine')

        # Deleting model 'Sequencing'
        db.delete_table(u'linapp_sequencing')

        # Removing M2M table for field samples on 'Sequencing'
        db.delete_table(db.shorten_name(u'linapp_sequencing_samples'))

        # Deleting model 'RawData'
        db.delete_table(u'linapp_rawdata')

        # Deleting model 'CorrectedRawData'
        db.delete_table(u'linapp_correctedrawdata')

        # Deleting model 'FailedTargetValue'
        db.delete_table(u'linapp_failedtargetvalue')

        # Deleting model 'SequenceDistribution'
        db.delete_table(u'linapp_sequencedistribution')

        # Deleting model 'TargetAnalysis'
        db.delete_table(u'linapp_targetanalysis')

        # Deleting model 'TargetVariant'
        db.delete_table(u'linapp_targetvariant')

        # Removing M2M table for field distribution on 'TargetVariant'
        db.delete_table(db.shorten_name(u'linapp_targetvariant_distribution'))

        # Deleting model 'GenSig'
        db.delete_table(u'linapp_gensig')

        # Removing M2M table for field variants on 'GenSig'
        db.delete_table(db.shorten_name(u'linapp_gensig_variants'))

        # Deleting model 'DM'
        db.delete_table(u'linapp_dm')

        # Deleting model 'CellTreeNode'
        db.delete_table(u'linapp_celltreenode')

        # Deleting model 'StorageType'
        db.delete_table(u'linapp_storagetype')

        # Deleting model 'StorageBox'
        db.delete_table(u'linapp_storagebox')

        # Deleting model 'PlateContext'
        db.delete_table(u'linapp_platecontext')

        # Deleting model 'PlatePlastica'
        db.delete_table(u'linapp_plateplastica')

        # Deleting model 'PlateType'
        db.delete_table(u'linapp_platetype')

        # Deleting model 'Plate'
        db.delete_table(u'linapp_plate')

        # Deleting model 'PlateStorage'
        db.delete_table(u'linapp_platestorage')

        # Deleting model 'SampleLocation'
        db.delete_table(u'linapp_samplelocation')

        # Deleting model 'PrimersMultiplex'
        db.delete_table(u'linapp_primersmultiplex')

        # Removing M2M table for field primers on 'PrimersMultiplex'
        db.delete_table(db.shorten_name(u'linapp_primersmultiplex_primers'))


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
            'comment': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
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
            'strand': ('django.db.models.fields.CharField', [], {'max_length': '1', 'null': 'True'}),
            'tail': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.PrimerTail']", 'null': 'True'}),
            u'target_ptr': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['linapp.Target']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'linapp.primersmultiplex': {
            'Meta': {'object_name': 'PrimersMultiplex'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'primers': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['linapp.TargetEnrichment']", 'symmetrical': 'False'})
        },
        u'linapp.primertail': {
            'Meta': {'object_name': 'PrimerTail'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'tail': ('django.db.models.fields.CharField', [], {'max_length': '50', 'null': 'True'})
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
            'amplicon': ('django.db.models.fields.CharField', [], {'max_length': '500'}),
            'chromosome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['linapp.Chromosome']"}),
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