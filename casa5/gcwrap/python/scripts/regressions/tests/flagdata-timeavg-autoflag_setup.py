
import sys
import os
import string
from locatescript import copydata
from locatescript import locatescript
from casa_stack_manip import stack_frame_find

gl = stack_frame_find()


def description():
    return "Regression for flagdata with time averaging and auto-flagging"


def data():
    return ['Four_ants_3C286.ms', '3ctst.ms']


def run(fetch=False):
    #Fetch Data

    if fetch:
        for f in data():
            copydata(f, os.getcwd())

    lepath = locatescript('flagdata-timeavg-autoflag_regression.py')
    print('Script used is ', lepath)
    gl['regstate'] = True
    execfile(lepath, gl)
    print('regstate = ', gl['regstate'])
    if not gl['regstate']:
        raise Exception('regstate = False')