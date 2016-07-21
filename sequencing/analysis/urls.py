
from django.conf.urls import url

from sequencing.analysis import views

urlpatterns = [
    url(r'^(?P<ngs_run>[a-zA-Z0-9-_]+)/summary.xlsx$', views.summary_table),
]
