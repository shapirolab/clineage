from django.contrib.auth.decorators import login_required
from django.shortcuts import render_to_response, redirect
from django.contrib.auth.forms import UserCreationForm
from django.http import HttpResponseRedirect, Http404

def register(request):
    if request.method == 'POST':
        form = UserCreationForm(request.POST)
        if form.is_valid():
            new_user = form.save()
            return HttpResponseRedirect("/")
    else:
        form = UserCreationForm()
    return render_to_response("registration/register.html", {
        'form': form,
        })

@login_required
def profile(request):
    user_profile = request.user.userprofile
    url = user_profile.url
