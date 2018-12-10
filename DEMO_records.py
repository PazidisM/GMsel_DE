# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 22:08:10 2018

@author: Marios D. Pazidis
"""
#!/usr/bin/env pypy
#%reset -f
import time
import Database_functions
import Split_functions
import folder_fmt_functions
import DE_functions
import DE_functions_cython_test
import UIn

start_time=time()

SaveFolder='C:\Marios\Research\Tests\GM_selection_DEMO\Python\GMsel_DE\TEST'
folders=folder_fmt_functions.folder_init(SaveFolder)
selectionParams,allowedRecs,DE_par,split_data=UIn.User_Input()
'''
#################################   Database screening / Combinations    #################################
'''

## screen database & define target spectrum
Ndatabase, Rec_db_metadata, Sa, Periods=Database_functions.screen_database(selectionParams,allowedRecs)
Sa_Tgt=Database_functions.EC8_Elastic_Spectrum_Type_1(selectionParams['soiltype'],selectionParams['a_g'],selectionParams['zeta'],Periods['T_all'])

## Define globals for simplicity / truncate arrays
T=Periods['T_match']
Sa_Tgt=Sa_Tgt[:,Periods['idx_T_match'][0]]
Sa=Sa[:,Periods['idx_T_match'][0]]

## Calculate max/min scaling factors and further screen the database
Sa,Rec_db_metadata,sf_ind, selectionParams, Ndatabase=Database_functions.Ind_sc_factor(Sa,Rec_db_metadata,Sa_Tgt,selectionParams)

## Calculate combinations / Create GM suites
Combs, Sa_unsc_ave, NSeed, split_data=Database_functions.Combinations(selectionParams,Rec_db_metadata,Ndatabase,Sa_Tgt,Sa,split_data,sf_ind)

## Printing formats
formats=folder_fmt_functions.fmt(SaveFolder,split_data,NSeed,Rec_db_metadata,DE_par)

## Save data to files
Split_functions.Pre_run_split(split_data,NSeed,Combs,Sa_unsc_ave,folders,formats)
del Combs, Sa_unsc_ave

'''
#################################   DE    #################################
#'''
#
## Initialize populations
DE_functions.Initialization(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)
end_time=time()
dur=end_time-start_time
print(dur)


start_time=time.clock()
## Differential Evolution

DE_functions.jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)
#DE_functions_cython_test.jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)

end_time=time.clock()
dur=end_time-start_time
print(dur)


start_time=time.clock()
## Differential Evolution

#DE_functions.jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)
DE_functions_cython_test.jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)

end_time=time.clock()
dur2=end_time-start_time
print(dur2)

print('cython is '+str((dur-dur2)/dur*100 )+'% faster')


