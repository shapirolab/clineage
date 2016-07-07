import os
import uuid
import contextlib
import shutil

from django.conf import settings


def get_unique_path(ext=None):
    unique = str(uuid.uuid4())
    parts = (unique[:2], unique[2:4], unique[4:6])
    fname = unique[6:]
    if ext is not None:
        fname += ".{}".format(ext)
    subdir = os.path.join(settings.DATA_STORE, *parts)
    os.makedirs(subdir, exist_ok=True)
    return os.path.join(subdir, fname)


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


class NotCreated(Exception):
    pass


def raise_or_create(model, **kwargs):
    """
    get_or_create, returning the created object if it is created, or raising
    a NotCreated with the gotten object otherwise. Useful for the unique_*_cms
    above.
    """
    obj, c = model.objects.get_or_create(**kwargs)
    if c:
        return obj
    else:
        raise NotCreated(obj)


def get_get_or_create(do_raise_or_create, model, **kwargs):
    """
    """
    try:
        return model.objects.get(**kwargs)
    except model.DoesNotExist:
        try:
            return do_raise_or_create()
        except NotCreated as nc:
            return nc.args[0]
