from mad_core import *
from src.madx.objects import mad_gvar


def mad_init(gvar):
    madx_start(gvar)


# TODO: Redefine function to take in gvar and filepath argument instead of int
def mad_run(gvar, file=0):
    madx_input(file, gvar)


def mad_finish(gvar):
    madx_finish(gvar)


global_variable = mad_gvar.GlobalVariable()
mad_init(global_variable)
mad_run(global_variable)
mad_finish(global_variable)
