import json
from functools import wraps
from datetime import datetime
import csv
import time

from dojango.decorators import json_response, expect_post_request
from .datables import *
from .tabselection import *
from .dataminers import *

from django.shortcuts import render, redirect, get_object_or_404
from django.db.models import Count
from django.http import Http404, HttpResponse
from django.core.urlresolvers import reverse, reverse_lazy
from django.views.generic.edit import CreateView, UpdateView, DeleteView
from django.views.generic.detail import DetailView
from django.template import loader

from clineage import settings
from utils.wells import num2abc
from utils.user_cells_report import user_cells_table_values, get_partner_report
from linapp.forms import PlateInputForm, MultipleCellForm
from wet_storage.models import SampleLocation, Plate, PlateStorage, PlatePlastica
from sampling.models import FACS, LaserCapture, SampleComposition
from lib_prep.workflows.models import CellContentType, CellContent
from misc.models import Taxa
from genomes.models import Assembly
from targeted_enrichment.planning.models import Microsatellite, TargetEnrichment
from utils.plate_string_description import *

@wraps(loader.render_to_string)
def render_to_response(*args,**kwargs):
    """
    Override bug in django.shortcuts.
    """
    #FIXME
    content_type = kwargs.pop("content_type",None)
    status = kwargs.pop("status",None)
    return HttpResponse(loader.render_to_string(*args,**kwargs),content_type=content_type,status=status)


def homepage(request, tab='individuals'):
    # from utils.parse_multiplexes import *
    # test_multiplex()
    #try:
        #publicexp = Experiment.objects.filter(is_public=True)
        #privateexp = Experiment.objects.filter(users=request.user)
    #except Experiment.DoesNotExist:
        #raise Http404
    selectedtab = gethomepagetab(tab)
    individualstable = getindividualstable()
    if individualstable.willHandle(request):
        return individualstable.handleRequest(request)
    extractioneventstable = getextractioneventstable()
    if extractioneventstable.willHandle(request):
        return extractioneventstable.handleRequest(request)
    extractionstable = getextractionstable()
    if extractionstable.willHandle(request):
        return extractionstable.handleRequest(request)
    samplingeventstable = getsamplingeventstable()
    if samplingeventstable.willHandle(request):
        return samplingeventstable.handleRequest(request)
    samplestable = getsamplestable()
    if samplestable.willHandle(request):
        return samplestable.handleRequest(request)
    platestable = getplatestable()
    if platestable.willHandle(request):
        return platestable.handleRequest(request)
    return render_to_response(
        'homepage.html', context={
            'selectedtab': selectedtab,
            'individualstable': individualstable,
            'extractionstable': extractionstable,
            'extractioneventstable': extractioneventstable,
            'samplingeventstable': samplingeventstable,
            'samplestable': samplestable,
            'platestable': platestable
        },request=request)



class JsonFormMixin(object):
    def form_valid(self, form):
        html = super(JsonFormMixin, self).form_valid(form)
        return HttpResponse(json.dumps("ok"), content_type="application/json")


    def form_invalid(self, form):
        return HttpResponse(json.dumps(form.errors), status=401,
                            content_type="application/json")


class IndividualCreate(JsonFormMixin, CreateView):
    model = Individual
    fields = "__all__"

class IndividualUpdate(JsonFormMixin, UpdateView):
    model = Individual
    fields = "__all__"

class IndividualDetail(DetailView):
    model = Individual

class IndividualDelete(DeleteView):
    model = Individual

class ExtractionEventCreate(JsonFormMixin, CreateView):
    model = ExtractionEvent
    fields = "__all__"

class ExtractionEventUpdate(JsonFormMixin, UpdateView):
    model = ExtractionEvent
    fields = "__all__"

class ExtractionEventDetail(DetailView):
    model = ExtractionEvent

class ExtractionEventDelete(DeleteView):
    model = ExtractionEvent

class ExtractionCreate(JsonFormMixin, CreateView):
    model = Extraction
    fields = "__all__"

class ExtractionDetail(DetailView):
    model = Extraction

class ExtractionUpdate(JsonFormMixin, UpdateView):
    model = Extraction
    fields = "__all__"

class ExtractionDelete(DeleteView):
    model = Extraction

class SamplingEventCreate(JsonFormMixin, CreateView):
    model = SamplingEvent
    fields = "__all__"

class SamplingEventUpdate(JsonFormMixin, UpdateView):
    model = SamplingEvent
    fields = "__all__"

class SamplingEventDetail(DetailView):
    model = SamplingEvent

class SamplingEventDelete(DeleteView):
    model = SamplingEvent

class FACSCreate(JsonFormMixin, CreateView):
    model = FACS
    fields = "__all__"

class FACSUpdate(JsonFormMixin, UpdateView):
    model = FACS
    fields = "__all__"

class FACSDetail(DetailView):
    model = FACS

class FACSDelete(DeleteView):
    model = FACS

class LaserCaptureCreate(JsonFormMixin, CreateView):
    model = LaserCapture
    fields = "__all__"

class LaserCaptureUpdate(JsonFormMixin, UpdateView):
    model = LaserCapture
    fields = "__all__"

class LaserCaptureDetail(DetailView):
    model = LaserCapture

class LaserCaptureDelete(DeleteView):
    model = LaserCapture

class CellCreate(JsonFormMixin, CreateView):
    model = Cell
    fields = "__all__"

class CellUpdate(JsonFormMixin, UpdateView):
    model = Cell
    fields = "__all__"
    template_name_suffix = '_update_form'

class CellDetail(DetailView):
    model = Cell

class CellDelete(DeleteView):
    model = Cell

class CellContentCreate(JsonFormMixin, CreateView):
    model = CellContent
    fields = "__all__"

class CellContentUpdate(JsonFormMixin, UpdateView):
    model = CellContent
    fields = "__all__"

class CellContentDetail(DetailView):
    model = CellContent

class CellContentDelete(DeleteView):
    model = CellContent


class PlateCreate(JsonFormMixin, CreateView):
    model = Plate
    fields = "__all__"

class PlateUpdate(JsonFormMixin, UpdateView):
    model = Plate
    fields = "__all__"

class PlateDelete(DeleteView):
    model = Plate

def multiple_cells_create(request):
    if request.method == 'POST':
        cells_form = MultipleCellForm(request.POST, prefix='multicell')
        if cells_form.is_valid():
            print(cells_form.cleaned_data)
            individual = cells_form.cleaned_data['individual']
            sampling_event = cells_form.cleaned_data['sampling']
            cell_composition = cells_form.cleaned_data['composition']
            cell_status = cells_form.cleaned_data['status']
            new_cells = []
            for i in range(1, cells_form.cleaned_data['copies']+1):
                new_cells.append(Cell(individual=individual,
                                      sampling=sampling_event,
                                      name=cells_form.cleaned_data['cells_name_prefix'] + str(i),
                                      composition=cell_composition,
                                      status=cell_status,
                                      comment=cells_form.cleaned_data['comment'],))
            Cell.objects.bulk_create(new_cells)
            return render(request, 'multiple_cells_form.html', {'form': cells_form})  # TODO: success page
    else:
        cells_form = MultipleCellForm(prefix='multicell')
    return render(request, 'multiple_cells_form.html', {'form': cells_form})


def plate_input(request):
    if request.method == 'POST':
        plate_form = PlateInputForm(request.POST, prefix='platecells')
        if plate_form.is_valid():
            plate_content_type = CellContentType.objects.get()
            individual = plate_form.cleaned_data['individual']
            sampling_event = plate_form.cleaned_data['sampling']
            inserting_user = plate_form.cleaned_data['user']
            cell_content_protocol = plate_form.cleaned_data['protocol']

            existing_plate = plate_form.cleaned_data['existing_plate']
            if existing_plate:
                # todo:improve DRY
                try:
                    existing_plate_location = PlateStorage.objects.get(plate=existing_plate)
                    existing_plate_location.storageBox = plate_form.cleaned_data['storage_box']
                    existing_plate_location.inner_location = plate_form.cleaned_data['inner_location']
                    existing_plate_location.notes = plate_form.cleaned_data['notes']
                    existing_plate_location.save()
                except PlateStorage.DoesNotExists:
                    new_plate_location = PlateStorage.objects.create(
                        storageBox=plate_form.cleaned_data['storage_box'],
                        plate=existing_plate,
                        inner_location=plate_form.cleaned_data['inner_location'],
                        notes=plate_form.cleaned_data['notes']
                    )
                plate = existing_plate
            elif plate_form.cleaned_data['plate_type'] and plate_form.cleaned_data['plate_name']:
                plate_type = plate_form.cleaned_data['plate_type']
                new_plate = Plate.objects.create(
                    name=plate_form.cleaned_data['plate_name'],
                    type=plate_type,
                    timestamp=plate_form.cleaned_data['timestamp']
                )
                new_plate_location = PlateStorage.objects.create(
                    storageBox=plate_form.cleaned_data['storage_box'],
                    plate=new_plate,
                    inner_location=plate_form.cleaned_data['inner_location'],
                    notes=plate_form.cleaned_data['notes']
                )
                plate = new_plate

            wells = plate_parser(plate_form.cleaned_data['samples_in_wells'])
            for index, cell_value in wells:
                well_composition = well_value_to_composition(cell_value)
                if well_composition:
                    new_cell = Cell.objects.create(
                        individual=individual,
                        sampling=sampling_event,
                        name=plate_form.cleaned_data['cells_name_prefix'] + index2str(index),
                        composition=well_composition,
                        comment=plate_form.cleaned_data['comment'] + '\r\ncomposition comment:%s' % cell_value,
                    )
                    new_cell_content = CellContent.objects.create(
                        cell=new_cell,
                        type=plate_content_type,
                        name=plate_form.cleaned_data['cells_name_prefix'] + index2str(index),
                        protocol=cell_content_protocol,
                        #seq_ready=False,
                        user=inserting_user,
                        comment=plate_form.cleaned_data['comment']
                    )
                    try:
                        # TODO: check if updating the pre-existing well's content is the desired behaviour
                        existing_location = SampleLocation.objects.get(plate=plate, well=index2str(index))
                        existing_location.sample = new_cell_content
                        existing_location.timestamp = plate_form.cleaned_data['timestamp']
                        existing_location.save()
                    except SampleLocation.DoesNotExist:
                        #plate.place_reagent(new_cell_content)
                        SampleLocation.objects.create(
                            plate=plate,
                            well=index2str(index),
                            reagent=new_cell_content,
                            timestamp=plate_form.cleaned_data['timestamp']
                        )
            return render(request, 'plate_cells_form.html', {'form': plate_form})
    else:
        plate_form = PlateInputForm(prefix='platecells')
    return render(request, 'plate_cells_form.html', {'form': plate_form})


def plate_input_with_names(request):
    if request.method == 'POST':
        plate_form = PlateInputForm(request.POST, prefix='platecells')
        if plate_form.is_valid():
            plate_content_type = CellContentType.objects.get()
            individual = plate_form.cleaned_data['individual']
            sampling_event = plate_form.cleaned_data['sampling']
            inserting_user = plate_form.cleaned_data['user']
            cell_content_protocol = plate_form.cleaned_data['protocol']

            existing_plate = plate_form.cleaned_data['existing_plate']
            if existing_plate:
                # todo:improve DRY
                try:
                    existing_plate_location = PlateStorage.objects.get(plate=existing_plate)
                    existing_plate_location.storageBox = plate_form.cleaned_data['storage_box']
                    existing_plate_location.inner_location = plate_form.cleaned_data['inner_location']
                    existing_plate_location.notes = plate_form.cleaned_data['notes']
                    existing_plate_location.save()
                except PlateStorage.DoesNotExists:
                    new_plate_location = PlateStorage.objects.create(
                        storageBox=plate_form.cleaned_data['storage_box'],
                        plate=existing_plate,
                        inner_location=plate_form.cleaned_data['inner_location'],
                        notes=plate_form.cleaned_data['notes']
                    )
                plate = existing_plate
            elif plate_form.cleaned_data['plate_type'] and plate_form.cleaned_data['plate_name']:
                plate_type = plate_form.cleaned_data['plate_type']
                new_plate = Plate.objects.create(
                    name=plate_form.cleaned_data['plate_name'],
                    type=plate_type,
                    timestamp=plate_form.cleaned_data['timestamp']
                )
                new_plate_location = PlateStorage.objects.create(
                    storageBox=plate_form.cleaned_data['storage_box'],
                    plate=new_plate,
                    inner_location=plate_form.cleaned_data['inner_location'],
                    notes=plate_form.cleaned_data['notes']
                )
                plate = new_plate
            else:
                print('invalid choice in plate form')
                raise
            wells = plate_parser(plate_form.cleaned_data['samples_in_wells'])
            for index, cell_value in wells:
                if cell_value:
                    new_cell = Cell.objects.create(
                        individual=individual,
                        sampling=sampling_event,
                        name=cell_value,
                        composition=SampleComposition.objects.get(name='Single Cell'),
                        comment=plate_form.cleaned_data['comment'] + '\r\ncomposition comment:%s' % cell_value,
                    )
                    new_cell_content = CellContent.objects.create(
                        cell=new_cell,
                        type=plate_content_type,
                        name=cell_value,
                        protocol=cell_content_protocol,
                        #seq_ready=False,
                        user=inserting_user,
                        comment=plate_form.cleaned_data['comment']
                    )
                    try:
                        # TODO: check if updating the pre-existing well's content is the desired behaviour
                        existing_location = SampleLocation.objects.get(plate=plate, well=index2str(index))
                        existing_location.sample = new_cell_content
                        existing_location.timestamp = plate_form.cleaned_data['timestamp']
                        existing_location.save()
                    except SampleLocation.DoesNotExist:
                        #plate.place_reagent(new_cell_content)
                        SampleLocation.objects.create(
                            plate=plate,
                            well=index2str(index),
                            reagent=new_cell_content,
                            timestamp=plate_form.cleaned_data['timestamp']
                        )
            return render(request, 'plate_cells_form.html', {'form': plate_form})
    else:
        plate_form = PlateInputForm(prefix='platecells')
    return render(request, 'plate_cells_form.html', {'form': plate_form})




def targets_tdv(request, taxa, assem):
    taxa_object = get_object_or_404(Taxa, taxonomy_id=taxa)
    assembly = get_object_or_404(Assembly, taxa=taxa_object, friendly_name=assem)
    s = 'target_name\trepeat_unit_len\trepeat_number\trepeat_unit\ttarget_start_pos\ttarget_end_pos\ttarget_length\t'
    s += 'fwd_primer\tfwd_start_pos\tfwd_end_pos\trev_primer\trev_start_pos\trev_end_pos\tplate_name\twell\t'
    s += 'passed_validation\tvalidation_failure\tvalidation_date\r\n'
    for enrichment in TargetEnrichment.objects.filter(physical_locationsleft__assembly=assembly)\
            .prefetch_related('targets__microsatellite', 'targets__referencevalue', 'physical_locations',
                              'left__sequence', 'right__sequence'):
        for target in enrichment.targets.all():
            s += target.name + '\t'
            try:
                mstarget = target.microsatellite
                s += str(mstarget.repeat_unit_len) + '\t'
                s += str(mstarget.repeat_number) + '\t'
                s += mstarget.repeat_unit_type + '\t'
            except Microsatellite.DoesNotExist:
                s += '\t\t\t'
            s += str(target.start_pos) + '\t'
            s += str(target.end_pos) + '\t'
            s += str(target.referencevalue.length) + '\t'
            s += enrichment.left.sequence.sequence + '\t'
            s += str(enrichment.left.start_pos) + '\t'
            s += str(enrichment.left.end_pos) + '\t'
            s += enrichment.right.sequence.sequence + '\t'
            s += str(enrichment.right.start_pos) + '\t'
            s += str(enrichment.right.end_pos) + '\t'
            locations = enrichment.physical_locations.all()
            if len(locations) == 1:
                s += locations[0].plate.name + '\t'
                s += locations[0].well + '\t'
            else:
                s += '\t\t'
            s += str(enrichment.passed_validation) + '\t'
            s += str(enrichment.validation_failure_id) + '\t'
            s += str(enrichment.validation_date) + '\t'
            s += '\r\n'
    return HttpResponse(s, content_type="text/plain")


def existing_primer_pairs_tdv(request, taxa, assem):
    taxa_object = get_object_or_404(Taxa, taxonomy_id=taxa)
    assembly = get_object_or_404(Assembly, taxa=taxa_object, friendly_name=assem)
    pcr_with_tails_type = TargetEnrichmentType.objects.get(pk=2)#TEMP
    s = 'target_name\trepeat_unit_len\trepeat_number\trepeat_unit\ttarget_start_pos\ttarget_end_pos\ttarget_length\t'
    s += 'fwd_primer\trev_primer\tplate_name\twell\tpassed_validation\tvalidation_failure\t'
    s += 'validation_date\tmultiplex_group\r\n'
    for enrichment in TargetEnrichment.objects.filter(left__assembly=assembly)\
            .annotate(locations_count=Count('physical_locations'))\
            .filter(locations_count__gt=0)\
            .prefetch_related('targets__microsatellite', 'targets__referencevalue', 'physical_locations',
                              'left__sequence', 'right__sequence')\
            .filter(type=pcr_with_tails_type):
            #.filter(type = pcr_with_tails_type) is temporary
        for target in enrichment.targets.all():
            s += target.name + '\t'
            try:
                mstarget = target.microsatellite
                s += str(mstarget.repeat_unit_len) + '\t'
                s += str(mstarget.repeat_number) + '\t'
                s += mstarget.repeat_unit_type + '\t'
            except Microsatellite.DoesNotExist:
                s += '\t\t\t'
            s += str(target.start_pos) + '\t'
            s += str(target.end_pos) + '\t'
            s += str(target.referencevalue.length) + '\t'
            s += enrichment.left.sequence.sequence + '\t'
            s += enrichment.right.sequence.sequence + '\t'
            locations = enrichment.physical_locations.all()
            if len(locations) == 1:
                s += locations[0].plate.name + '\t'
                s += locations[0].well + '\t'
            else:
                s += '\t\t'
            s += str(enrichment.passed_validation) + '\t'
            s += str(enrichment.validation_failure_id) + '\t'
            s += str(enrichment.validation_date) + '\t'
            s += '\r\n'
    return HttpResponse(s, content_type="text/plain")


def plate_well_selection(request, plate_plastica_id):
    plate_plastica = get_object_or_404(PlatePlastica, pk=plate_plastica_id)
    rows = [num2abc(row) for row in range(1, plate_plastica.rows+1)]
    columns = list(range(1, plate_plastica.columns+1))
    return render_to_response('forms/plate_boolean_form.html', context={'rows': rows, 'columns': columns},request=request)


def partner_cells_table_view(request, partner_name,
                             cell_folder=settings.S_MAIN,
                             individual_name=None,
                             palette_name='hls'):
    response = HttpResponse(content_type='text/csv')
    try:
        p = User.objects.get(username__contains=partner_name)
        response['Content-Disposition'] = 'attachment; filename="{}_cell_data.csv"'.format(partner_name)
        fieldnames = ['CellContent ID',
                      'Cell ID',
                      'Cell Name',
                      'Cell Type',
                      'Cell Group',
                      'Plate',
                      'Well',
                      'Plate Location',
                      'Group Color',
                      'Sampling Event',
                      'Sampling Comment',
                      'FACS Marker',
                      'Organ',
                      'Tissue',
                      'Sample Name',
                      'Sample Comment',
                      'Extraction Event',
                      'Extraction Event Comment',
                      'Individual Name',
                      'Individual Comment',
                      'Gender',
                      'Sequencing File Name',
                      ]
        writer = csv.DictWriter(response, fieldnames=fieldnames)
        writer.writeheader()
        for cell_values in user_cells_table_values(partner_name, individual_name, cell_folder):
                writer.writerow(cell_values)
    except User.DoesNotExist:
        raise Http404("No Partner names matches the given query.")
    return response


def partner_cells_html_view(request, partner_name, individual_name=None, palette_name='hls'):
    # response = HttpResponse(content_type='text/html')
    # response['content_type'] = 'application/liquid; charset=utf-8'
    partner_dict = get_partner_report(partner_name, individual_name, palette_name=palette_name)
    return render_to_response('user_report.html', {'partner_dict': partner_dict})
    # return HttpResponse(response)

