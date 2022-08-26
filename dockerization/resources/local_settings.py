import os.path
import os

#CHROMOSOMES_PATH = '/net/mraid11/export/dcstor/Ofir/GenomesData/'
CHROMOSOMES_PATH = '/clineage/hg19/'
DEBUG=True
# S_MAIN = '/net/mraid11/export/dcstor/LINEAGE/Hiseq/NSR2/fastq_human/Output'
S_MAIN = '/data_store'
# DATABASES = {
#     'default': {
#         'ENGINE': 'django.db.backends.mysql', # Add 'postgresql_psycopg2', 'postgresql', 'mysql', 'sqlite3' or 'oracle'.
#         'NAME': 'clineage_prod',                      # Or path to database file if using sqlite3.
#         'USER': 'dcsoft',                      # Not used with sqlite3.
#         'PASSWORD': '164d8ae81bd7e38a163ea2e144114b25',                  # Not used with sqlite3.
#         'HOST': 'db',                      # Set to empty string for localhost. Not used with sqlite3.
#         'PORT': '3306',                      # Set to empty string for default. Not used with sqlite3.
#     }
# }
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql', # Add 'postgresql_psycopg2', 'postgresql', 'mysql', 'sqlite3' or 'oracle'.
        'HOST': os.environ['MYSQL_HOST'],     # Set to empty string for localhost. Not used with sqlite3.
        'PORT': os.environ['MYSQL_PORT'],     # Set to empty string for default. Not used with sqlite3.
        'NAME': os.environ['MYSQL_DATABASE'], # Or path to database file if using sqlite3.
        'USER': os.environ['MYSQL_USER'],     # Not used with sqlite3.
        'PASSWORD': os.environ['MYSQL_PASSWORD'], # Not used with sqlite3.
    }
}

ALLOWED_HOSTS = [
    '*',
]

SECRET_KEY = 'gcd6^pq#l78%_jg3@2jf)csum=_d%5-pa5!3hma4g&e*lxx+m$'

DATA_STORE = '/data_store'
# NOA_MATLAB = r'/home/dcsoft/s/Ofir/noa_matlab/Code/'
NOA_MATLAB = r''
# IGORS_CODE = r'/home/dcsoft/s/Ofir/igor_tree_reconstruction_20170808/'
