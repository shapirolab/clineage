import os
import uuid
import contextlib
import shutil

from django.conf import settings


def get_unique_path(ext=None):
    if ext is None:
        fmt = "{}"
    else:
        fmt = "{{}}.{}".format(ext)
    return os.path.join(settings.DATA_STORE, fmt.format(uuid.uuid4()))


@contextlib.contextmanager
def unique_file_cm(ext):
    """
    Context manager for a unique path to a file, which will be deleted if an
    error occurs.
    """
    fname = get_unique_path(ext)
    try:
        yield fname
    except:
        try:
            os.unlink(fname)
        except:
            pass
        raise


@contextlib.contextmanager
def unique_dir_cm():
    """
    Context manager for a unique path to a dir, which will be created, and
    will be recursively deleted if an error occurs.
    """
    dirname = get_unique_path()
    os.mkdir(dirname)
    try:
        yield dirname
    except:
        try:
            shutil.rmtree(dirname)
        except:
            pass
        raise


@contextlib.contextmanager
def unlink(fname):
    try:
        yield fname
    finally:
        os.unlink(fname)
