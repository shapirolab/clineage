
from django.db import migrations, models
from django.db.migrations.operations.base import Operation

class AlterContentType(Operation):

    reduces_to_sql = True
    reversible = True

    def __init__(self,from_model,to_app,to_model):
        self.from_model = from_model
        self.to_app = to_app
        self.to_model = to_model

    def state_forwards(self,app_label,state):
        pass

    def database_forwards(self,app_label,schema_editor,from_state,to_state):
        self._operate(app_label,self.to_app,self.from_model,self.to_model,schema_editor,from_state)

    def _operate(self,from_app,to_app,from_model,to_model,schema_editor,from_state):
        ContentType = from_state.apps.get_model("contenttypes", "ContentType")
        db_alias = schema_editor.connection.alias
        ct, c = ContentType.objects.using(db_alias).get_or_create(app_label=from_app,model=from_model)
        # This is a hack because django ContentType populating is somewhat broken.
        try:
            to_ct = ContentType.objects.using(db_alias).get(app_label=to_app,model=to_model)
        except ContentType.DoesNotExist:
            pass
        else:
            to_ct.delete()
        ct.app_label = to_app
        ct.model = to_model
        ct.save()

    def database_backwards(self,app_label,schema_editor,from_state,to_state):
        self._operate(self.to_app,app_label,self.to_model,self.from_model,schema_editor,from_state)

    def describe(self):
        return "Change a row in the contenttypes table."

