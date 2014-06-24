
from django.conf.urls import patterns, include, url
from linapp.views import *
from django.views.generic.base import RedirectView
from django.core.urlresolvers import reverse

urlpatterns = patterns('',
    url(r'^experiment/(?P<exp_id>\d+)/$', 'linapp.views.experiment'),
    url(r'^experiment/(?P<exp_id>\d+)/tab:(?P<tab>\w+)$', 'linapp.views.experiment'),
    url(r'^experiment/workflow/(?P<exp_id>\d+)$', 'linapp.views.experimentworkflow'),
    url(r'^comment/(?P<exp_id>\d+)/$', 'linapp.views.commentPost'),
    url(r'^postexperiment/(?P<exp_id>\d+)/$', 'linapp.views.experimentPost'),
    url(r'^postfile/(?P<exp_id>\d+)/$', 'linapp.views.filesPost'),
    url(r'^algorithms/(?P<alg_id>\d+)/$', algorithm, name='algorithm-details'),
    url(r'^algorithms/(?P<alg_id>\d+)/tab:(?P<tab>\w+)$', 'linapp.views.algorithm'),
    url(r'^treedata/$', 'linapp.views.testtree'),
    url(r'^algform/$', 'linapp.views.algrunform'),
    url(r'^algform/(?P<alg_id>\d+)$', 'linapp.views.algrunform'),
    url(r'^taxa:(?P<taxa>\d+)/assembly:(?P<assem>\w+)$', 'linapp.views.targets_tdv'),
    url(r'^existing/taxa:(?P<taxa>\d+)/assembly:(?P<assem>\w+)$', 'linapp.views.existing_primer_pairs_tdv'),

    ###############################################################################################
    ##################################Forms urls for dialogs#######################################
    ###############################################################################################

    ## *the relation between the urls of "path_to_object/add" and "path_to_object/<pk>" is crucial
    ##  for the dialog-on-row-click template-implemented client-side functionality since the
    ##  templates use "path_to_object/add"[:-3]+pk as a client side hack.  TODO: reconsider that

    #individual
    url(r'^forms/individual/add$',IndividualCreate.as_view(), name='individual_add'),
    url(r'^details/individual/(?P<pk>\d+)$',IndividualUpdate.as_view(), name='individual_detail'),
    url(r'^forms/individual/(?P<pk>\d+)$',IndividualUpdate.as_view(), name='individual_update'),
    url(r'^forms/individual/delete/(?P<pk>\d+)$',IndividualDelete.as_view(), name='individual_delete'),

    #extraction event
    url(r'^forms/extractionevent/add$',ExtractionEventCreate.as_view(), name='extraction_event_add'),
    url(r'^details/extractionevent/(?P<pk>\d+)$',ExtractionEventUpdate.as_view(), name='extraction_event_detail'),
    url(r'^forms/extractionevent/(?P<pk>\d+)$',ExtractionEventUpdate.as_view(), name='extraction_event_update'),
    url(r'^forms/extractionevent/delete/(?P<pk>\d+)$',ExtractionEventDelete.as_view(), name='extraction_event_delete'),

    #extraction
    url(r'^forms/extraction/add$',ExtractionCreate.as_view(), name='extraction_add'),
    url(r'^details/extraction/(?P<pk>\d+)$',ExtractionUpdate.as_view(), name='extraction_detail'),
    url(r'^forms/extraction/(?P<pk>\d+)$',ExtractionUpdate.as_view(), name='extraction_update'),
    url(r'^forms/extraction/delete/(?P<pk>\d+)$',ExtractionDelete.as_view(), name='extraction_delete'),

    #sampling event
    url(r'^forms/samplingevent/add$',SamplingEventCreate.as_view(), name='samplingevent_add'),
    url(r'^details/samplingevent/(?P<pk>\d+)$',SamplingEventUpdate.as_view(), name='samplingevent_detail'),
    url(r'^forms/samplingevent/(?P<pk>\d+)$',SamplingEventUpdate.as_view(), name='samplingevent_update'),
    url(r'^forms/samplingevent/delete/(?P<pk>\d+)$',SamplingEventDelete.as_view(), name='samplingevent_delete'),
        #FACS
    url(r'^forms/facs/add$',FACSCreate.as_view(), name='facs_add'),
    url(r'^details/facs/(?P<pk>\d+)$',FACSUpdate.as_view(), name='facs_detail'),
    url(r'^forms/facs/(?P<pk>\d+)$',FACSUpdate.as_view(), name='facs_update'),
    url(r'^forms/facs/delete/(?P<pk>\d+)$',FACSDelete.as_view(), name='facs_delete'),
        #laser capture
    url(r'^forms/lasercapture/add$',LaserCaptureCreate.as_view(), name='lasercapture_add'),
    url(r'^details/lasercapture/(?P<pk>\d+)$',LaserCaptureUpdate.as_view(), name='lasercapture_detail'),
    url(r'^forms/lasercapture/(?P<pk>\d+)$',LaserCaptureUpdate.as_view(), name='lasercapture_update'),
    url(r'^forms/lasercapture/delete/(?P<pk>\d+)$',LaserCaptureDelete.as_view(), name='lasercapture_delete'),

    #sample
    url(r'^forms/cell/add$',CellCreate.as_view(), name='cell_add'),
    url(r'^details/cell/(?P<pk>\d+)$',CellUpdate.as_view(), name='cell_detail'),
    url(r'^forms/cell/(?P<pk>\d+)$',CellUpdate.as_view(), name='cell_update'),
    url(r'^forms/cell/delete/(?P<pk>\d+)$',CellDelete.as_view(), name='cell_delete'),

    #multiple cells
    url(r'^forms/cell/add_multiple_cells$', 'linapp.views.multiple_cells_create', name='cell_multiple_add'),
    url(r'^forms/cell/add_cells_plate$', 'linapp.views.plate_input', name='cells_plate_add'),
    url(r'^forms/cell/add_cells_plate_with_names$', 'linapp.views.plate_input_with_names', name='cells_plate_add_with_names'),

    #algorithm
    url(r'^forms/algorithm/add$',AlgorithmCreate.as_view(), name='algorithm_add'),
    url(r'^details/algorithm/(?P<pk>\d+)$',AlgorithmUpdate.as_view(), name='algorithm_detail'),
    url(r'^forms/algorithm/(?P<pk>\d+)$',AlgorithmUpdate.as_view(), name='algorithm_update'),
    url(r'^forms/algorithm/delete/(?P<pk>\d+)$',AlgorithmDelete.as_view(), name='algorithm_delete'),

    #plate
    url(r'^forms/plate/add$',PlateCreate.as_view(), name='plate_add'),
    url(r'^details/plate/(?P<pk>\d+)$',PlateUpdate.as_view(), name='plate_detail'),
    url(r'^forms/plate/(?P<pk>\d+)$',PlateUpdate.as_view(), name='plate_update'),
    url(r'^forms/plate/delete/(?P<pk>\d+)$',PlateDelete.as_view(), name='plate_delete'),

    url(r'^plate_wells_selection/(?P<plate_plastica_id>\d+)$', 'linapp.views.plate_well_selection'),

    # url(r'^cell/$','linapp.views.cellform'),
    # url(r'^cell/(?P<cell_id>\d+)$','linapp.views.cellform'),
    url(r'^experiment/(?P<exp_id>\d+)/member/$', 'linapp.views.memberform'),
    url(r'^experiment/(?P<exp_id>\d+)/member/(?P<mem_id>\d+)$', 'linapp.views.memberform'),
)