from django.conf.urls import *

urlpatterns = patterns('',
    (r'^login/$',
     'django.contrib.auth.views.login',
         {'template_name': 'registration/login.html'}),

    (r'^logout/$',
     'django.contrib.auth.views.logout',
         {'template_name': 'registration/logout.html'}),

    (r'^password_change/$',
     'django.contrib.auth.views.password_change',
#         {'template_name': 'registration/password_change_form.html'}
        ),

    (r'^password_change/done/$',
     'django.contrib.auth.views.password_change_done',
#         {'template_name': 'registration/password_change_done.html'}
        ),

    (r'^password_reset/$',
     'django.contrib.auth.views.password_reset',
#         {'template_name': 'registration/password_reset_form.html',
#          'email_template_name': 'accounts/password_reset_email.html'}
        ),

    (r'^password_reset/done/$',
     'django.contrib.auth.views.password_reset_done',
#         {'template_name': 'registration/password_reset_done.html'}
        ),

    (r'^reset/(?P<uidb36>[0-9A-Za-z]+)-(?P<token>.+)/$',
     'django.contrib.auth.views.password_reset_confirm',
#         {'template_name': 'registration/password_reset_confirm.html'}
        ),

    (r'^reset/done/$',
     'django.contrib.auth.views.password_reset_complete',
#         {'template_name': 'registration/password_reset_complete.html'}
        ),

    (r'^register/$',
     'accounts.views.register'
        ),


#
#    (r'^register/done/$',
#     'linapp.views.signup_done',
#         {'template_name': 'accounts/signup_done.html'}),
#
#    (r'^signup/(?P<uidb36>[0-9A-Za-z]+)-(?P<token>.+)/$',
#     'mysite.accounts.views.signup_confirm'),
#
#    (r'^signup/complete/$',
#     'mysite.accounts.views.signup_complete',
#         {'template_name': 'accounts/signup_complete.html'}),
)