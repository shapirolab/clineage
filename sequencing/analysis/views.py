
from django.http import FileResponse, HttpResponseNotFound

def summary_table(request, ngs_run):
    return HttpResponseNotFound(b"Sorry, run does not exist.")
