from mad_core import *
import mad_gvar


def mad_init():
    madx_start(gvar)

# TODO: Redefine function to take in gvar and filepath argument instead of int
def mad_run(gvar, file=0):
    madx_input(file)

def mad_finish(gvar):
    madx_finish()


gvar = mad_gvar.GlobalVariable()
mad_init()
mad_run()
mad_finish()