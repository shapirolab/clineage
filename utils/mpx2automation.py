

from linapp.models import PrimersMultiplex

def func(multiplexes, mpx_plate):
    primers_mapping = []
    for multiplex in multiplexes:
        for primers_pair in multiplex.primers.all():
            pp_locations = primers_pair.physical_locations.all()
            if len(pp_locations) == 0:
                print 'reference to non existing primer pair %s' % primers_pair.left+'_'+primers_pair.right
                raise
            if len(pp_locations) == 1:
                pp_location = pp_locations[0]
                primers_mapping.append((pp_location.plate.name, pp_location.well, mpx_plate.plate.name, mpx_plate.well))
            else:
                #handle multiple locations
                raise
    primers_mapping.sort(key=lambda tup: tup[0])
    return primers_mapping