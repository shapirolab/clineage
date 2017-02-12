
from django.conf.urls import url

from linapp import views

urlpatterns = [
    url(r'^taxa:(?P<taxa>\d+)/assembly:(?P<assem>\w+)$', views.targets_tdv),
    url(r'^existing/taxa:(?P<taxa>\d+)/assembly:(?P<assem>\w+)$', views.existing_primer_pairs_tdv),

    ###############################################################################################
    ##################################Forms urls for dialogs#######################################
    ###############################################################################################

    ## *the relation between the urls of "path_to_object/add" and "path_to_object/<pk>" is crucial
    ##  for the dialog-on-row-click template-implemented client-side functionality since the
    ##  templates use "path_to_object/add"[:-3]+pk as a client side hack.  TODO: reconsider that

    #individual
    url(r'^forms/individual/add$',views.IndividualCreate.as_view(), name='individual_add'),
    url(r'^details/individual/(?P<pk>\d+)$',views.IndividualUpdate.as_view(), name='individual_detail'),
    url(r'^forms/individual/(?P<pk>\d+)$',views.IndividualUpdate.as_view(), name='individual_update'),
    url(r'^forms/individual/delete/(?P<pk>\d+)$',views.IndividualDelete.as_view(), name='individual_delete'),

    #extraction event
    url(r'^forms/extractionevent/add$',views.ExtractionEventCreate.as_view(), name='extraction_event_add'),
    url(r'^details/extractionevent/(?P<pk>\d+)$',views.ExtractionEventUpdate.as_view(), name='extraction_event_detail'),
    url(r'^forms/extractionevent/(?P<pk>\d+)$',views.ExtractionEventUpdate.as_view(), name='extraction_event_update'),
    url(r'^forms/extractionevent/delete/(?P<pk>\d+)$',views.ExtractionEventDelete.as_view(), name='extraction_event_delete'),

    #extraction
    url(r'^forms/extraction/add$',views.ExtractionCreate.as_view(), name='extraction_add'),
    url(r'^details/extraction/(?P<pk>\d+)$',views.ExtractionUpdate.as_view(), name='extraction_detail'),
    url(r'^forms/extraction/(?P<pk>\d+)$',views.ExtractionUpdate.as_view(), name='extraction_update'),
    url(r'^forms/extraction/delete/(?P<pk>\d+)$',views.ExtractionDelete.as_view(), name='extraction_delete'),

    #sampling event
    url(r'^forms/samplingevent/add$',views.SamplingEventCreate.as_view(), name='samplingevent_add'),
    url(r'^details/samplingevent/(?P<pk>\d+)$',views.SamplingEventUpdate.as_view(), name='samplingevent_detail'),
    url(r'^forms/samplingevent/(?P<pk>\d+)$',views.SamplingEventUpdate.as_view(), name='samplingevent_update'),
    url(r'^forms/samplingevent/delete/(?P<pk>\d+)$',views.SamplingEventDelete.as_view(), name='samplingevent_delete'),
        #FACS
    url(r'^forms/facs/add$',views.FACSCreate.as_view(), name='facs_add'),
    url(r'^details/facs/(?P<pk>\d+)$',views.FACSUpdate.as_view(), name='facs_detail'),
    url(r'^forms/facs/(?P<pk>\d+)$',views.FACSUpdate.as_view(), name='facs_update'),
    url(r'^forms/facs/delete/(?P<pk>\d+)$',views.FACSDelete.as_view(), name='facs_delete'),
        #laser capture
    url(r'^forms/lasercapture/add$',views.LaserCaptureCreate.as_view(), name='lasercapture_add'),
    url(r'^details/lasercapture/(?P<pk>\d+)$',views.LaserCaptureUpdate.as_view(), name='lasercapture_detail'),
    url(r'^forms/lasercapture/(?P<pk>\d+)$',views.LaserCaptureUpdate.as_view(), name='lasercapture_update'),
    url(r'^forms/lasercapture/delete/(?P<pk>\d+)$',views.LaserCaptureDelete.as_view(), name='lasercapture_delete'),

    #sample
    url(r'^forms/cell/add$',views.CellCreate.as_view(), name='cell_add'),
    url(r'^details/cell/(?P<pk>\d+)$',views.CellUpdate.as_view(), name='cell_detail'),
    url(r'^forms/cell/(?P<pk>\d+)$',views.CellUpdate.as_view(), name='cell_update'),
    url(r'^forms/cell/delete/(?P<pk>\d+)$',views.CellDelete.as_view(), name='cell_delete'),

    #multiple cells
    url(r'^forms/cell/add_multiple_cells$', views.multiple_cells_create, name='cell_multiple_add'),
    url(r'^forms/cell/add_cells_plate$', views.plate_input, name='cells_plate_add'),
    url(r'^forms/cell/add_cells_plate_with_names$', views.plate_input_with_names, name='cells_plate_add_with_names'),

    #plate
    url(r'^forms/plate/add$',views.PlateCreate.as_view(), name='plate_add'),
    url(r'^details/plate/(?P<pk>\d+)$',views.PlateUpdate.as_view(), name='plate_detail'),
    url(r'^forms/plate/(?P<pk>\d+)$',views.PlateUpdate.as_view(), name='plate_update'),
    url(r'^forms/plate/delete/(?P<pk>\d+)$',views.PlateDelete.as_view(), name='plate_delete'),

    # collaborators_reports
    # url(r'^csv_view/partner_name:(?P<partner_name>\w+)$', views.partner_cells_table_view),
    url(r'^csv_view/cell_data/partner_name:(?P<partner_name>\w+)/individual_name:(?P<individual_name>\w+)/ngsrun:(?P<ngsrun>\w+)$', views.partner_cells_table_view_db),
    url(r'^csv_view/cell_data/partner_name:(?P<partner_name>\w+)/ngsrun:(?P<ngsrun>\w+)$', views.partner_cells_table_view_db),
    url(r'^csv_view/cell_data/partner_name:(?P<partner_name>\w+)/individual_name:(?P<individual_name>\w+)$', views.partner_cells_table_view_db),
    url(r'^csv_view/cell_data/individual_name:(?P<individual_name>\w+)/ngsrun:(?P<ngsrun>\w+)$', views.partner_cells_table_view_db),

    url(r'^csv_view/cell_data/ngsrun:(?P<ngsrun>\w+)$', views.partner_cells_table_view_db),
    url(r'^csv_view/cell_data/partner_name:(?P<partner_name>\w+)$', views.partner_cells_table_view_db),
    url(r'^csv_view/cell_data/individual_name:(?P<individual_name>\w+)$', views.partner_cells_table_view_db),

    url(r'^csv_view/cell_data/$', views.partner_cells_table_view_db),
    url(r'^csv_view/cell_data$', views.partner_cells_table_view_db),
    # url(r'^csv_view/partner_name:(?P<partner_name>\w+)/individual_name:(?P<individual_name>\w+)$',
    #     views.partner_cells_table_view),
    # url(r'^csv_view/partner_name:(?P<partner_name>\w+)/palette_name:(?P<palette_name>\w+)$', views.partner_cells_table_view),
    # url(r'^csv_view/partner_name:(?P<partner_name>\w+)/individual_name:(?P<individual_name>\w+)/palette_name:(?P<palette_name>\w+)$',
    #     views.partner_cells_table_view),
    #url(r'^csv_view/partner_name:(?P<partner_name>\w+))

    url(r'^html_view/partner_name:(?P<partner_name>\w+)$', views.partner_cells_html_view),
    url(r'^html_view/partner_name:(?P<partner_name>\w+)/individual_name:(?P<individual_name>\w+)$',
        views.partner_cells_html_view),
    url(r'^html_view/partner_name:(?P<partner_name>\w+)/palette_name:(?P<palette_name>\w+)$', views.partner_cells_html_view),
    url(r'^html_view/partner_name:(?P<partner_name>\w+)/individual_name:(?P<individual_name>\w+)/palette_name:(?P<palette_name>\w+)$',
        views.partner_cells_html_view),

]
