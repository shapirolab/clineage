import datetime


def getScriptHeader():
    res ='############  DOCS  ###############################\n\
DOC\n\
DNE Script\n\
Project Name: GFPRun\n\
Created on'+ str(datetime.datetime.now())+''+'\n'\
'ENDDOC\n'\
\
'############  TABLE  ###############################\n\
TABLE TABLE_DNE\n\
\
############  HEADER  ###############################\n\
OPTION UPDATE_VOL_DB\n\
OPTION VERIFY_TABLE\n\
\
#GLOBAL REAGENTS\n\
REAGENT TE T10 1 PIE_AUTBOT 2\n\
REAGENT Lambda_Mix_X5 T10 3 PIE_AUTBOT 2\n\
REAGENT ELN_Mix_X5 T10 5 PIE_AUTBOT 2\n\
REAGENT GIBSON_Mix_X2 T10 15 PIE_AUTBOT 1\n\
\
REAGENT LB_SYBR T10 7 PIE_AUTBOT_SLOW 2\n\
REAGENT LIZ_500 T10 9 PIE_AUTBOT_SLOW 2\n\
\
REAGENT PHOS_MIX T10 11 PIE_AUTBOT 2\n\
\
REAGENT TE_BUF T10 13 PIE_AUTBOT 2\n\
\
REAGENT DDW	BUF12 1 AUT_AIR 8\n\
REAGENT LB BUF12 9 PIE_TROUGH_AUTAIR 8\n\
\
#Wetting Parameters\n\
WET_MIX_VOL = 150\n\
WET_MIX_TIMES = 15\n\
\
#Plate reader inspection\n\
RDR_DDW_VOL = 55\n\
\
ELN_OLIGO_MIX_VOL = 4\n\
\
LIZ_500_VOL = 12\n\
CE_Sample = 1.5\n\
CE_SEQ_DIL_VOL = 28.5\n\
\
PCR_FLOUR_VOL = 4\n\
\
############# END OF DEFAULT HEADER #########################\n\
# anything from here on is produced by the automation,\n\
# parameters should be changed through AUT_CONFIG only\n'
    return res
