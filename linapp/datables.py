#from linapp.models import *
from datable.web.table import Table
from datable.web.storage import Storage
from datable.web.columns import StringColumn, DateTimeColumn
from datable.web.widgets import StringWidget
from django.utils.translation import ugettext as _
from datable.web.extra.widgets import DateTimeLessOrEqual
from datable.web.extra.widgets import ForeignKeyComboBox
from datable.web.extra.widgets import DateTimeGreaterOrEqual
from django.core.urlresolvers import reverse
from django.contrib.auth.models import User

from sampling.models import Individual, Extraction, ExtractionEvent, \
    SamplingEvent, Cell, SampleComposition, SampleStatus, Location, \
    Organ, Tissue
from wet_storage.models import PlateType, Plate

#def getindividualstable(exp=None):
    #if exp:
        #query = exp.individuals()
    #else:
        #query = Individual.objects.order_by('-pk')
def getindividualstable():
    query = Individual.objects.order_by('-pk')
    return Table(
        name='individualstable',
        storage=Storage(
            querySet=query,
            columns=[
                StringColumn('name', width=100),
                StringColumn('taxa', width=100),
                StringColumn('sex', width=50),
                DateTimeColumn('born', width=100),
                StringColumn('comment', width=150),
                StringColumn('background', width=100),
                StringColumn('location', width=100),
                DateTimeColumn('sacrificed', width=100),
                ],
            widgets=[
                StringWidget('name', placeholder=_("Name")),
                DateTimeLessOrEqual('born', paired=True),
                DateTimeGreaterOrEqual('born', paired=True),
            ]
        ),
    )


def getextractioneventstable():
    return Table(
        name='extractioneventstable',
        storage=Storage(
            querySet=ExtractionEvent.objects.order_by('-pk'),
            columns=[
                StringColumn('individual', width=100),
                StringColumn('name', width=150),
                StringColumn('user_performed', width=150),
                StringColumn('user_documented', width=150),
                DateTimeColumn('date', width=150),
                StringColumn('location', width=150),
                StringColumn('comment', width=150),
                ],
            widgets=[
                StringWidget('name', placeholder=_("Name")),
                DateTimeLessOrEqual('date', paired=True, placeholder=_("Extracted after date")),
                DateTimeGreaterOrEqual('date', paired=True, placeholder=_("Extracted before date")),
                StringWidget('comment', placeholder=_("Comment")),
                ForeignKeyComboBox(
                    'location',
                    otherSet=Location.objects.all(),
                    otherField='name',
                    otherFormat='Location: %(name)s',
                    placeholder=_("Extraction event location")),
                ForeignKeyComboBox(
                    'individual',
                    otherSet=Individual.objects.all(),
                    otherField='name',
                    otherFormat='Individual: %(name)s',
                    placeholder=_("individual's name")),
                ForeignKeyComboBox(
                    'user_performed',
                    otherSet=User.objects.all(),
                    otherField='user_performed',
                    otherFormat='User: %(username)s',
                    placeholder=_("user who extracted")),
                ]
        ),
    )


def getextractionstable():
    return Table(
        name='extractionstable',
        storage=Storage(
            querySet=Extraction.objects.order_by('-pk'),
            columns=[
                StringColumn('extraction_event', width=100),
                StringColumn('name', width=150),
                StringColumn('organ', width=150),
                StringColumn('tissue', width=150),
                StringColumn('comment', width=150),
                ],
            widgets=[
                StringWidget('name', placeholder=_("Name")),
                StringWidget('comment', placeholder=_("Comment")),
                ForeignKeyComboBox(
                    'organ',
                    otherSet=Organ.objects.all(),
                    otherField='name',
                    otherFormat='Organ: %(name)s',
                    placeholder=_("Organ")),
                ForeignKeyComboBox(
                    'tissue',
                    otherSet=Tissue.objects.all(),
                    otherField='name',
                    otherFormat='Tissue: %(name)s',
                    placeholder=_("Tissue")),
                ForeignKeyComboBox(
                    'extraction_event',
                    otherSet=Individual.objects.all(),
                    otherField='name',
                    otherFormat='Extraction event: %(name)s',
                    placeholder=_("extraction event")),
                ]
        ),
    )

def getsamplingeventstable():
    return Table(
        name='samplingeventstable',
        storage=Storage(
            querySet=SamplingEvent.objects.order_by('-pk'),
            columns=[
                StringColumn('name', width=100),
                StringColumn('extraction', width=150),
                DateTimeColumn('date', width=150),
                StringColumn('user', width=150),
                StringColumn('comment', width=150),
                ],
            widgets=[
                StringWidget('name', placeholder=_("Name")),
                StringWidget('comment', placeholder=_("Comment")),
                ForeignKeyComboBox(
                    'extraction',
                    otherSet=Extraction.objects.all(),
                    otherField='name',
                    otherFormat='Extraction: %(name)s',
                    placeholder=_("source extraction")),
                ForeignKeyComboBox(
                    'user',
                    otherSet=User.objects.all(),
                    otherField='username',
                    otherFormat='User: %(username)s',
                    placeholder=_("user who sampled")),
                ]
        ),
    )
def getsamplestable():
    return Table(
        name='samplestable',
        storage=Storage(
            querySet=Cell.objects.order_by('-pk').select_related(),
            columns=[
                StringColumn('name', width=100),
                StringColumn('sampling', width=150),
                StringColumn('composition', width=150),
                StringColumn('status', width=150),
                StringColumn('comment', width=150),
                ],
            widgets=[
                StringWidget('name', placeholder=_("Name")),
                StringWidget('comment', placeholder=_("Comment")),
                ForeignKeyComboBox(
                    'sampling',
                    otherSet=SamplingEvent.objects.all(),
                    otherField='name',
                    otherFormat='Sampling event: %(name)s',
                    placeholder=_("sampling event")),
                ForeignKeyComboBox(
                    'composition',
                    otherSet=SampleComposition.objects.all(),
                    otherField='name',
                    otherFormat='Composition: %(name)s',
                    placeholder=_("sample composition")),
                ForeignKeyComboBox(
                    'status',
                    otherSet=SampleStatus.objects.all(),
                    otherField='name',
                    otherFormat='Status: %(name)s',
                    placeholder=_("sample status")),
                ]
        ),
    )


#def getpDNAtable(exp):
    #return Table(
        #name='pDNAtable',
        #storage=Storage(
            #querySet=exp.cellscontents().select_related(),
            #columns=[
                #StringColumn('name', width=100),
                #StringColumn('parent', width=100),
                #StringColumn('cell', width=150),
                #StringColumn('panel', width=150),
                #StringColumn('type', width=150),
                #StringColumn('protocol', width=150),
                #StringColumn('user', width=150),
                #StringColumn('comment', width=150),
                #],
        #),
    #)


def getplatestable():
    return Table(
        name='platestable',
        storage=Storage(
            querySet=Plate.objects.all(),
            columns=[
                StringColumn('name', width=100),
                StringColumn('type', width=150),
                StringColumn('barcode', width=150),
                DateTimeColumn('timestamp', width=100),
                StringColumn('state', width=100),
                StringColumn('lastusedwell', width=100),
                ],
            widgets=[
                StringWidget('name', placeholder=_("Name")),
                StringWidget('state', placeholder=_("State")),
                ForeignKeyComboBox(
                    'type',
                    otherSet=PlateType.objects.all(),
                    otherField='friendly',
                    otherFormat='Type: %(friendly)s',
                    placeholder=_("Type")),
                ]
        ),
        objectpath='/CLineage/plate/',
    )


#def getsequencingtable(exp):
    #return Table(
        #name='sequencingtable',
        #storage=Storage(
            #querySet=exp.sequencings(),
            #columns=[
                #StringColumn('sample', width=100),
                #StringColumn('data', width=150),
                #StringColumn('machine', width=150),
                #StringColumn('protocol', width=150),
                #StringColumn('user', width=150),
                #DateTimeColumn('date', width=150),
                #],
        #),
    #)
