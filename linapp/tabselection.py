


def gethomepagetab(tab): #TODO: change to multiple returns and 404 as endcase
    selectedtab = 'tab1'
    if tab == 'individuals':
        selectedtab = 'tab1'
    if tab == 'extractionevents':
        selectedtab = 'tab2'
    if tab == 'extractions':
        selectedtab = 'tab3'
    if tab == 'samplings':
        selectedtab = 'tab4'
    if tab == 'cells':
        selectedtab = 'tab5'
    if tab == 'algorithms':
        selectedtab = 'tab6'
    if tab == 'plates':
        selectedtab = 'tab7'
    return selectedtab


def getalgorithmtab(tab):
    selectedtab = 'tab1'
    if tab == 'general':
        selectedtab = 'tab1'
    if tab == 'parameters':
        selectedtab = 'tab2'
    if tab == 'runs':
        selectedtab = 'tab3'
    if tab == '?': #'?' is just a stub.
        selectedtab = 'tab4'
    return selectedtab


def getexperimenttab(tab):
    selectedtab = 'tab1'
    if tab == 'general':
        selectedtab = 'tab1'
    if tab == 'samples':
        selectedtab = 'tab2'
    if tab == 'processed':
        selectedtab = 'tab3'
    if tab == 'sequencing':
        selectedtab = 'tab4'
    if tab == 'targets':
        selectedtab = 'tab5'
    if tab == 'ts':
        selectedtab = 'tab6'
    if tab == 'tv':
        selectedtab = 'tab7'
    if tab == 'gs':
        selectedtab = 'tab8'
    if tab == 'dm':
        selectedtab = 'tab9'
    if tab == 'individuals':
        selectedtab = 'tab10'
    if tab == 'messages':
        selectedtab = 'tab11'
    if tab == 'map':
        selectedtab = 'tab12'
    if tab == 'members':
        selectedtab = 'tab13'
    if tab == 'repository':
        selectedtab = 'tab14'
    if tab == 'tree':
        selectedtab = 'tab15'
    return selectedtab