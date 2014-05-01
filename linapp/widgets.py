
from dojango.forms.widgets import *

class RevcoWidget(MultiWidget):
    def __init__(self, attrs=None):
        shelves = [(shelf, shelf) for shelf in ('A', 'B', 'C', 'D')]
        compartments = [(compartment, compartment) for compartment in ('1', '2', '3', '4', '5')]
        draws = [(draw, draw) for draw in ('1', '2', '3', '4', '5')]
        _widgets = (
            Select(attrs=attrs, choices=shelves),
            Select(attrs=attrs, choices=compartments),
            Select(attrs=attrs, choices=draws),
        )
        super(RevcoWidget, self).__init__(_widgets, attrs)

    def decompress(self, value):
        if value:
            return [value.shelf, value.compartment, value.draw]
        return [None, None, None]

    def format_output(self, rendered_widgets):
        return u''.join(rendered_widgets)

    def value_from_datadict(self, data, files, name):
        locationlist = [
            widget.value_from_datadict(data, files, name + '.%s' % i)
            for i, widget in enumerate(self.widgets)]
        print locationlist
        return 'test'


class FreezerWidget(MultiWidget):
    def __init__(self, attrs=None):
        shelves = [(shelf, shelf) for shelf in ('1', '2', '3', '4', '5')]
        _widgets = (
            Select(attrs=attrs, choices=shelves),
        )
        super(FreezerWidget, self).__init__(_widgets, attrs)

    def decompress(self, value):
        if value:
            return [value.shelf]
        return [None, None, None]

    def format_output(self, rendered_widgets):
        return u''.join(rendered_widgets)