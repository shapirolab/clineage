#!/usr/bin/env python
# encoding: utf8

from django.views.decorators.csrf import csrf_exempt
from linapp.models import *
from spyne.server.django import DjangoApplication
from spyne.model.primitive import String, Integer, AnyDict
from spyne.model.complex import Iterable
from spyne.service import ServiceBase
from spyne.interface.wsdl import Wsdl11
from spyne.protocol.soap import Soap11
from spyne.application import Application
from spyne.decorator import rpc
from spyne.error import *
from Bio.SeqUtils.MeltingTemp import Tm_staluc

class HelloWorldService(ServiceBase):
    @rpc(String, Integer, _returns=Iterable(String))
    def say_hello(ctx, name, times):
        for i in xrange(times):
            yield 'Hello, %s' % name

hello_world_service = csrf_exempt(DjangoApplication(Application([HelloWorldService],
    'spyne.examples.django',
    in_protocol=Soap11(),
    out_protocol=Soap11(),
    interface=Wsdl11(),
)))

class CLineageWebServices(ServiceBase):
    @rpc(Iterable(String), _returns=Iterable(Iterable(String)))
    def get_targets_data(ctx, target_names):
        '''
        Elaborate query for targets.
        MATLAB example: get_targets_data(service_obj, struct('string',{{'Seq05944'}}))

        output columns:
        # Target name
        # Target: MS/Other Mutation
        # Basic Unit size
        # Expected Number of repeats
        # Basic Unit Type
        # Chromosome
        # Length MS
        # Primer sequence -  Left
        # Primer Tm -  Left
        # Primer sequence -  Right
        # Primer Tm -  Right
        # Validation status
        # Target location on Chromosome - start
        # Target location on Chromosome - end
        # Amplicon location on Chromosome - start
        # Amplicon location on Chromosome - end
        # Mpx groups names
        # number of primer-pairs in the multiplex group

        <b>Parameters:</b>
        @param target_names list of target names to query.
        @return full target data in table format (see columns)
        '''
        tails_pcr = TargetEnrichmentType.objects.get(name='PCR_with_tails')
        targets = [target for target in target_names]
        for tgt in Target.objects.filter(name__in=targets):
            for te in tgt.targetenrichment_set.all().filter(type=tails_pcr):
                for mpx in te.primersmultiplex_set.all():
                    try:
                        ms = tgt.microsatellite
                        repeat_type = str(ms.repeat_type)
                        repeat_number = str(ms.repeat_number)
                        repeat_unit = str(ms.repeat_unit)
                    except Microsatellite.DoesNotExist:
                        repeat_type = ''
                        repeat_number = ''
                        repeat_unit = ''
                    yield [tgt.name,
                           tgt.type.name,  # Target: MS/Other Mutation
                           str(repeat_type),  # Basic Unit size
                           str(repeat_number),  # Expected Number of repeats
                           str(repeat_unit),  # Basic Unit Type
                           tgt.chromosome.name,  # Chromosome
                           str(tgt.end_pos-tgt.start_pos),  # Length MS
                           te.left.sequence.sequence,  # Primer sequence -  Left
                           str(Tm_staluc(te.left.referencevalue.sequence)),  # Primer Tm -  Left
                           te.right.sequence.sequence,  # Primer sequence -  Right
                           str(Tm_staluc(te.right.referencevalue.sequence)),  # Primer Tm -  Right
                           str(te.passed_validation),
                           str(tgt.start_pos),  # Target location on Chromosome - start
                           str(tgt.end_pos),  # Target location on Chromosome - end
                           str(te.left.start_pos),  # Amplicon location on Chromosome - start
                           str(te.right.end_pos),  # Amplicon location on Chromosome - end
                           mpx.name,  # Mpx groups names
                           str(len(mpx.primers.all()))]

    @rpc(String, _returns=Iterable(Iterable(String)))
    def get_targets_by_panel(ctx, panel_name):
        '''
        Elaborate query for targets.
        MATLAB example: get_targets_by_panel(service_obj, struct('string','Seq05944'))

        output columns:
        # Target name
        # Target: MS/Other Mutation
        # Basic Unit size
        # Expected Number of repeats
        # Basic Unit Type
        # Chromosome
        # Length MS
        # Primer sequence -  Left
        # Primer Tm -  Left
        # Primer sequence -  Right
        # Primer Tm -  Right
        # Validation status
        # Target location on Chromosome - start
        # Target location on Chromosome - end
        # Amplicon location on Chromosome - start
        # Amplicon location on Chromosome - end
        # Mpx groups names
        # number of primer-pairs in the multiplex group

        <b>Parameters:</b>
        @param target_names list of target names to query.
        @return full target data in table format (see columns)
        '''
        print panel_name
        try:
            panel = Panel.objects.get(name=panel_name)
        except Panel.DoesNotExist:
            yield ['No such panel']
            return

        for te in panel.targets.all():
            for tgt in te.targets.all():
                for mpx in te.primersmultiplex_set.all():
                    for loc in mpx.physical_locations.all():
                        try:
                            ms = tgt.microsatellite
                            repeat_type = str(ms.repeat_type)
                            repeat_number = str(ms.repeat_number)
                            repeat_unit = str(ms.repeat_unit)
                        except Microsatellite.DoesNotExist:
                            repeat_type = ''
                            repeat_number = ''
                            repeat_unit = ''
                        yield [tgt.name,
                               tgt.type.name,  # Target: MS/Other Mutation
                               str(repeat_type),  # Basic Unit size
                               str(repeat_number),  # Expected Number of repeats
                               str(repeat_unit),  # Basic Unit Type
                               tgt.chromosome.name,  # Chromosome
                               str(tgt.end_pos-tgt.start_pos),  # Length MS
                               te.left.sequence.sequence,  # Primer sequence -  Left
                               str(Tm_staluc(te.left.referencevalue.sequence)),  # Primer Tm -  Left
                               te.right.sequence.sequence,  # Primer sequence -  Right
                               str(Tm_staluc(te.right.referencevalue.sequence)),  # Primer Tm -  Right
                               str(te.passed_validation),
                               str(tgt.start_pos),  # Target location on Chromosome - start
                               str(tgt.end_pos),  # Target location on Chromosome - end
                               str(te.left.start_pos),  # Amplicon location on Chromosome - start
                               str(te.right.end_pos),  # Amplicon location on Chromosome - end
                               str(mpx.name),  # Mpx groups names
                               str(len(mpx.primers.all())),
                               str(loc.plate.name),
                               str(loc.well)]

    @rpc(String, String, Integer, Integer, _returns=String)
    def get_genomic_sequence(ctx, assembly_name, chromosome_name, start_index, end_index):
        '''
        Chromosome sequence query. Like genome browser DNA view.
        <b>Parameters</b>
        @param assembly_name Assembly name in short format only e.g. 'mm9', 'hg19'...
        @param chromosome_name Chromosome name in string format only e.g. 'X', '3'...
        @param start_index Start index in integer format
        @param end_index End index in integer format
        @return Genomic sequence between indices
        '''
        print type(assembly_name), type(chromosome_name), type(start_index), type(end_index)
        try:
            assembly = Assembly.objects.get(friendly_name=assembly_name)
        except Assembly.DoesNotExist:
            return ResourceNotFoundError(faultstring='Assembly \'%s\' was not found.'%assembly_name)
        try:
            chromosome = Chromosome.objects.get(assembly=assembly, name=chromosome_name)
        except Chromosome.DoesNotExist:
            return ResourceNotFoundError(faultstring='Chromosome \'%s\' was not found.'%chromosome_name)
        return chromosome.getdna(start_index, end_index)


# @rpc(Iterable(String), _returns=Iterable(Iterable(String)))
#     def save_panel(ctx, target_names):
#         '''
#         Elaborate query for targets.
#         MATLAB example: get_targets_data(service_obj, struct('string',{{'Seq05944'}}))
#
#         output columns:
#         # Target name
#
#
#         <b>Parameters:</b>
#         @param target_names list of target names to query.
#         @return full target data in table format (see columns)
#         '''
#         targets = [target for target in target_names]
#         for tgt in Target.objects.filter(name__in=targets):
#             for te in tgt.targetenrichment_set.all():
#                 for mpx in te.primersmultiplex_set.all():
#                     try:
#                         ms = tgt.microsatellite
#                         repeat_type = str(ms.repeat_type)
#                         repeat_number = str(ms.repeat_number)
#                         repeat_unit = str(ms.repeat_unit)
#                     except Microsatellite.DoesNotExist:
#                         repeat_type = ''
#                         repeat_number = ''
#                         repeat_unit = ''
#                     yield [tgt.name,
#                            tgt.type.name,  # Target: MS/Other Mutation
#                            str(repeat_type),  # Basic Unit size
#                            str(repeat_number),  # Expected Number of repeats
#                            str(repeat_unit),  # Basic Unit Type
#                            tgt.chromosome.name,  # Chromosome
#                            str(tgt.end_pos-tgt.start_pos),  # Length MS
#                            te.left.sequence.sequence,  # Primer sequence -  Left
#                            str(Tm_staluc(te.left.sequence.sequence)),  # Primer Tm -  Left
#                            te.right.sequence.sequence,  # Primer sequence -  Right
#                            str(Tm_staluc(te.right.sequence.sequence)),  # Primer Tm -  Right
#                            str(te.passed_validation),
#                            str(tgt.start_pos),  # Target location on Chromosome - start
#                            str(tgt.end_pos),  # Target location on Chromosome - end
#                            str(te.left.start_pos),  # Amplicon location on Chromosome - start
#                            str(te.right.end_pos),  # Amplicon location on Chromosome - end
#                            mpx.name,  # Mpx groups names
#                            str(len(mpx.primers.all()))]

niki_service = csrf_exempt(DjangoApplication(Application([CLineageWebServices],
    'spyne.examples.django',
    in_protocol=Soap11(),
    out_protocol=Soap11(),
    interface=Wsdl11(),
)))

    #
# 2.	Cell query -
# This query is by cell name (still there is no official template for the names - in the experiments). The fields are:
# 'If it’s a duplicate of not and with whom' – needs some sort of a system of flags maybe, to combine two duplicates.
# 'Method of extraction'
# 'Method of preparation'
# 'Singe Cell / Bulk'
# 'From which organ'
# 'To which organisms it belongs'
#
    # @rpc(Iterable(String), _returns=Iterable(Iterable(String)))
    # def get_cell_contents_data(ctx, cell_names):
    #     cells = [cell for cell in cell_names]
    #     for cell in Cell.objects.filter(name__in=cells):
    #         for cell_content in cell.cellcontent_set.all():
    #             try:
    #                 cell.sampling.
    #             yield [cell.name,
    #                    cell.sampling]


#__________________________________________________________________________
# class Cell(models.Model):
#     sampling = models.ForeignKey(SamplingEvent)
#     name = models.CharField(max_length=50)
#     experiment = models.ManyToManyField(Experiment, related_name='cells', null=True, blank=True)
#     composition = models.ForeignKey(SampleComposition)#single cell or bulk
#     status = models.ForeignKey(SampleStatus, null=True, blank=True)
#     comment = models.TextField(null=True, blank=True)
#
#     def __unicode__(self):
#         return self.name
#
#     def get_absolute_url(self):
#         return reverse('cell_detail', kwargs={'pk': self.pk})
# ## -------------------------------------------------------------------------------------
# class CellContent(models.Model):  # aka DNA
#     parent = models.ForeignKey('CellContent', null=True, blank=True)
#     cell = models.ForeignKey(Cell)
#     panel = models.ForeignKey(Panel, null=True, blank=True)
#     type = models.ForeignKey(CellContentType)
#     name = models.CharField(max_length=50, null=True, blank=True)
#     protocol = models.ForeignKey(Protocol, null=True, blank=True)
#     seq_ready = models.BooleanField(default=False)
#     user = models.ForeignKey(User, null=True, blank=True)
#     comment = models.TextField()
#     physical_locations = generic.GenericRelation('SampleLocation',
#                                content_type_field='content_type',
#                                object_id_field='object_id')

#__________________________________________________________________________


# class NGS_raw(ServiceBase):
#     @rpc(String, String, _returns=String)
#     def ngs_raw(ctx, sequencing_name, username):
#         try:
#             seq_event = Sequencing.objects.get(name = sequencing_name)
#         except ObjectDoesNotExist:
#             yield 'Error: no such sequencing name %s' %sequencing_name
#         try:
#             usr = User.objects.get(name=username)
#         except ObjectDoesNotExist:
#             yield 'Error: no such username %s' %username
#         RawData.objects.create()
#
# hello_world_service = csrf_exempt(DjangoApplication(Application([HelloWorldService],
#     'spyne.examples.django',
#     in_protocol=Soap11(),
#     out_protocol=Soap11(),
#     interface=Wsdl11(),
# )))


#

# 3.	Experiment query -
# This query will refer to the cell files made by Adam. Basically, I will have the experiment name and it should give me all the cells names and the path to the '.hist' and '.summary' files that are made by Adam program.
# I think it will be the last step in the progress though, since it is something that will need Adam and Tamirs' collaboration.
