
from django.conf.urls import include, url
from django.views.generic import TemplateView
from django.contrib import admin

import clineage.settings as settings

# Uncomment the next two lines to enable the admin:
#admin.autodiscover()

urlpatterns = [
    url(r'^CLineage/', include('linapp.urls')),

    #Homepage
    url(r'^$', 'linapp.views.homepage'),
    url(r'^tab:(?P<tab>\w+)$', 'linapp.views.homepage'),

    # Login/Logout/Register Account
    url(r'^accounts/', include('accounts.urls')),

    # Dojango path
    url(r'^dojango/', include('dojango.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),

    #TriStateCheckbox tests
    url('^success_page_background/$', TemplateView.as_view(template_name='success_page_background.html')),

    #soap test
    url(r'^soap/', include('soap.urls')),
]

