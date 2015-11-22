from django.conf.urls import url
from django.contrib.auth import views as auth_views

from accounts.views import register

urlpatterns = [
    url(r'^login/$',
     auth_views.login,
         {'template_name': 'registration/login.html'}, name='login'),

    url(r'^logout/$',
     auth_views.logout,
         {'template_name': 'registration/logout.html'}, name='logout'),

    url(r'^password_change/$',
     auth_views.password_change, name='password_change'
#         {'template_name': 'registration/password_change_form.html'}
        ),

    url(r'^password_change/done/$',
     auth_views.password_change_done, name='password_change_done'
#         {'template_name': 'registration/password_change_done.html'}
        ),

    url(r'^password_reset/$',
     auth_views.password_reset, name='password_reset'
#         {'template_name': 'registration/password_reset_form.html',
#          'email_template_name': 'accounts/password_reset_email.html'}
        ),

    url(r'^password_reset/done/$',
     auth_views.password_reset_done, name='password_reset_done'
#         {'template_name': 'registration/password_reset_done.html'}
        ),

    url(r'^reset/(?P<uidb36>[0-9A-Za-z]+)-(?P<token>.+)/$',
     auth_views.password_reset_confirm, name='password_reset_confirm'
#         {'template_name': 'registration/password_reset_confirm.html'}
        ),

    url(r'^reset/done/$',
     auth_views.password_reset_complete, name='password_reset_complete'
#         {'template_name': 'registration/password_reset_complete.html'}
        ),

    url(r'^register/$',
     register, name='register'
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
]
