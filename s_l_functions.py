# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 13:41:45 2018

@author: Marios D. Pazidis

"""
import shelve

def save(mainFolder):
    filename=mainFolder+'\\vars.out'
    my_shelf = shelve.open(filename,'n') # 'n' for new
    
    for key in dir():
        try:
            my_shelf[key] = globals()[key]
        except TypeError:
            #
            # __builtins__, my_shelf, and imported modules can not be shelved.
            #
            print('ERROR shelving: {0}'.format(key))
    my_shelf.close()

def load(mainFolder):
    filename=mainFolder+'\\vars.out'
    my_shelf = shelve.open(filename)
    for key in my_shelf:
        globals()[key]=my_shelf[key]
    my_shelf.close()