class Macro:
    def __init__(self, name, n_formal, dead, formal, tokens, body, stamp, original):
        self.name = name
        self.n_formal = n_formal
        self.dead = dead
        self.formal = formal
        self.tokens = tokens
        self.body = body
        self.stamp = stamp
        self.original = original


class MacroList:
    def __init__(self, max_length):
        self.lst = []
        self.max_length = max_length

    def add_var(self, var):
        self.lst += [var]



def make_macro(statemenet):
    return None

def new_macro(n_formal, length, p_length):
    return None

def new_macro_list(max_length):
    return MacroList(max_length)

def add_to_macro_list(macro, macro_list):
    return None

def disable_line(name, macro_list):
    return None

def replace_lines(macro, replace, reps):
    return None

def save_macros2file():
    return None

def exec_macro(in_cmd, pos):
    return None