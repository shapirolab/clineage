
def newickify(node): #TODO:This function should be a part of a printing lib
    if node.get_children():
        treestring = '('
        for child in node.get_children():
            treestring += newickify(child) + ','
        treestring = treestring[:-1] + ')'
        if node.is_root_node():
            return treestring + node.cell.name +';'
        return treestring+str(node.cell.name) + ':' + str(node.distance)
    else:
        return str(node.cell.name) + ':' + str(node.distance)

def workflowstore(exp): #TODO:This function should be a part of a printing lib
    individuals = exp.individuals()
    if not individuals:
        return dict()
    indiv = individuals[0] #TODO: support multiple individuals
    extractionrefs = []
    workflow = {'identifier': 'compid',
                'label': 'name',
                'items': [
                        {'compid': 'individual-'+str(indiv.pk), 'name': indiv.name, 'type': 'individual',
                         'children': extractionrefs}
                ]
    }
    for extraction in exp.extractions(): #TODO: filter by experiment ?
        samplerefs = []
        extractionrefs.append({'_reference': 'extraction-'+str(extraction.pk)})
        workflow['items'].append({'compid': 'extraction-'+str(extraction.pk), 'name': extraction.name, 'type': 'extraction',
                                  'children': samplerefs})
        for sample in extraction.samplingevent_set.all():
            cellrefs = []
            samplerefs.append({'_reference': 'sample-'+str(sample.pk)})
            workflow['items'].append({'compid': 'sample-'+str(sample.pk), 'name': sample.name, 'type': 'sample',
                                      'children': cellrefs})
            for cell in sample.cell_set.all():
                cellrefs.append({'_reference': 'cell-'+str(cell.pk)})
                workflow['items'].append({'compid': 'cell-'+str(cell.pk), 'name': cell.name, 'type': 'cell'})
    return workflow