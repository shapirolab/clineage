from Utils.script_header import *
# from Automation.robotic_scripts.getScriptHeader import *


MULTIPLEX_VOLUME = 3


def distribute_multiplex(test_input = [(('sp1','A1'),('dp','B1')),(('sp2','B3'),('dp','B4')),(('sp1','A2'),('dp','C5'))]):
    '''
    gets a list. each item in the list is a tuple.
    each tuple contains 2 tuples:
    first tuple: (src plate, scr well)
    second tupple: (dst plate, dst well).
    prints a roboease script that distributes 2 ul from each first tupple to each second tuple.
    '''

    total_volume = 142.5
    allowedPlates = ['P3', 'P4', 'P5', 'P6', 'P11', 'P13', 'P14', 'P15']
    inputPlates = []
    for tup in test_input:
        inputPlates.append(tup[0][0])
    inputPlates = list(set(inputPlates))
    #first script is for wetting
    script = getScriptHeader()
    script += 'LIST DDW_LIST\n'
    done_wells = []
    for i, tup in enumerate(test_input):#now distributing all multiplexes from src to dst
        volume = MULTIPLEX_VOLUME
        dst_well  = tup[1][1]
        for j, tup2 in enumerate(test_input):
            well = tup2[1][1]
            if j != i and well == dst_well:
                volume+=MULTIPLEX_VOLUME
        if dst_well not in done_wells:
            script+=str(total_volume - volume)+'\n'
            done_wells.append(dst_well)
    script += 'ENDLIST\n'
    script += '# Table Layout (for liquid handling)\n'
    dest_plate_name = test_input[0][1][0]
    script += 'LABWARE ' + dest_plate_name + ' P2  \"96 Well PCR Plate\"\n'  # dest plate is allways in P2
    script += '############  SCRIPT SECTION ###############################\nSCRIPT\n'
    #distributing DDW into detination wells
    script += 'DIST_REAGENT2 DDW '
    done_wells = []
    for i, tup in enumerate(test_input):  # now distributing all multiplexes from src to dst
        platePos = 'P2'
        dstWell = tup[1][1]
        if dstWell not in done_wells:
            script += platePos + ':' + dstWell + '+1'
            done_wells.append(dstWell)
            if i < len(test_input) - 1:
                script += ';'
    script += ' '+'DDW_LIST  AUT_AIR  TIPTYPE:200\n'
    script += '\nENDSCRIPT\n'
    print script
    for i in range(0, len(inputPlates)/len(allowedPlates) + 1):
        if i <= len(inputPlates)/len(allowedPlates):
            last_index = len(allowedPlates)
        else:
            last_index = len(inputPlates) % len(allowedPlates)
        subinputPlates = inputPlates[i*len(allowedPlates):i*len(allowedPlates)+last_index]
        script = getScriptHeader()
        srcPlates_srcPlatesOnTable = []  # creating a list of plates and their positions on table [(plate1,pos1),(plate2,pos2)]
        for i, plate in enumerate(subinputPlates):
            tup = (subinputPlates[i], allowedPlates[i])
            srcPlates_srcPlatesOnTable.append(tup)
        sub_test_input = []
        for i, tup in enumerate(test_input):  # getting the plates that are on the table right now
            for plate in subinputPlates:
                if tup[0][0] == plate:
                    sub_test_input.append(tup)
        script += '# Table Layout (for liquid handling)\n'
        dest_plate_name = test_input[0][1][0]
        script += 'LABWARE ' + dest_plate_name + ' P2  \"96 Well PCR Plate\"\n'  # dest plate is allways in P2
        for i, tup in enumerate(srcPlates_srcPlatesOnTable):
            script += 'LABWARE  ' + tup[0] + ' ' + tup[1] + ' ' + '\"96 Well PCR Plate\"\n'
        script += '############  SCRIPT SECTION ###############################\nSCRIPT\n'
        script += 'TRANSFER_LOCATIONS '
        for i, tup in enumerate(sub_test_input):  # now distributing all multiplexes from src to dst
            srcPlateName = tup[0][0]
            platePos = getPosFromPlate(plateName=srcPlateName,srcPlates_srcPlatesOnTable=srcPlates_srcPlatesOnTable)
            srcWell = tup[0][1]
            script += platePos + ':' + srcWell + '+1'
            if i < len(test_input) - 1:
                script += ';'
        script += '   '
        for i, tup in enumerate(sub_test_input):  # now distributing all multiplexes from src to dst
            platePos = 'P2'
            dstWell = tup[1][1]
            script+=platePos+':'+dstWell+'+1'
            if i < len(test_input) -1:
                script+=';'
        script+=' '+str(MULTIPLEX_VOLUME)+' BOT_BOT_LIQUID  TIPTYPE:20'
        script+='\nENDSCRIPT\n'
        print script


def getPosFromPlate(plateName, srcPlates_srcPlatesOnTable):
    res = None
    for tup in srcPlates_srcPlatesOnTable:
        if tup[0] == plateName:
            res = tup[1]
            break
    return res