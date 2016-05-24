import os
import uuid
from django.conf import settings


def get_unique_path(ext=None):
    if ext is None:
        fmt = "{}"
    else:
        fmt = "{{}}.{}".format(ext)
    return os.path.join(settings.DATA_STORE, fmt.format(uuid.uuid4()))
