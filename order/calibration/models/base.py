

class ModelParams(object):

    def get_for_x(self, x):
        raise NotImplementedError()

    class Model(object):

        def get_for_cycles(self, cycles):
            raise NotImplementedError()

        class CyclesModel(object):

            def get_hist_for_length(self, length):
                raise NotImplementedError()
