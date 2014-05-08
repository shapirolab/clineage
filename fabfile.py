from fabric.decorators import task
from contextlib import contextmanager
from fabric.api import *
from fabric.contrib.files import upload_template, append
import os.path

env.hosts = [
    '132.76.81.158'
]

PROJ_DIR = os.path.abspath(os.path.dirname(__file__))
CONF_DIR = os.path.abspath(os.path.join(PROJ_DIR, 'conf'))

# env.user = "ofirr"
env.gunicorn_port = 9000
env.code_dir = '~/CLineage/'
env.venv_dir = '~/.virtualenvs/cl/'
env.venv_command = '. {}bin/activate'.format(env.venv_dir)
env.log_dir = '~/logs/clineage/'
env.clone_url = "git@bitbucket.org:ofirr/clineage.git"
env.pidfile = '~/clineage.pid'
# env.pip_version = "1.5.4"


@contextmanager
def virtualenv():
    with cd(env.code_dir):
        with prefix(env.venv_command):
            yield

@task
def host_type():
    run('uname -s')


@task
def git_log():
    with virtualenv():
        run("git log -n 1")


@task
def freeze():
    with virtualenv():
        run("pip freeze")


@task
def reload_app():
    run("kill -HUP `cat %s`" % env.pidfile)


@task
def upgrade_pip():
    with virtualenv():
        run("pip install pip=={}".format(env.pip_version))


@task
def deploy(restart=True):
    # upgrade_pip()
    with virtualenv():
        run("git pull")
        run("pip install -r requirements.txt")
        run("pip install -r requirements-deploy.txt")
        run("python manage.py syncdb --noinput")
        run("python manage.py migrate --merge --noinput")
        run("python manage.py collectstatic --noinput")
        run("git log -n 1 --format='%ai %h' > static/version.txt")
        run("git log -n 1 > static/version-full.txt")
    if restart:
        reload_app()


# @task
# def hard_reload():
#     run("sudo supervisorctl restart opencommunity")


# @task
# def very_hard_reload():
#     run("sudo service supervisor stop")
#     run("sudo service supervisor start")


@task
def log():
    run("tail %s*" % env.log_dir)


@task
def clone_project():
    run("git clone %s %s" % (env.clone_url, env.code_dir))


@task
def create_venv():
    with cd(env.code_dir):
        # run("mkvirtualenv cl")
        run("mkdir -p ~/.virtualenvs/")
        run("cd ~/.virtualenvs/ && virtualenv cl")


@task
def createsuperuser():
    """ Creates a Django superuser for the project """
    with virtualenv():
        run("python manage.py createsuperuser")


# @task
# def supervisor_status():
#     """ Show server's supoervisord status """
#     run("sudo supervisorctl status")


@task
def switch(branch):
    """ fetches all branchs, and checkouts the specified git branch """
    with cd(env.code_dir):
        run('git fetch origin')
        run('git checkout {}'.format(branch))
        deploy()


@task
def showkeys():
    """ Displays authorized public ssh keys for user """
    with hide('stdout'):
        keys = run('cat .ssh/authorized_keys')
    print keys


@task
def push_key(key_file):
    """ Appends an ssh public key file from the specified file
    """
    with open(key_file) as f:
        key_text = f.read()
    append('~/.ssh/authorized_keys', key_text)