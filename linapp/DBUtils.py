from linapp.models import *
from django.contrib.auth.models import User

class DBUtils:
    @staticmethod
    def getrole(usr,exp):
        if usr.is_superuser:
            return LineageRole.emptySuperUserRole()
        else:
            try:
                return ExperimentUser.objects.get(user=usr,experiment=exp).role
            except ExperimentUser.DoesNotExist:
                return LineageRole.emptyRole()

