# TODO: Implement on-demand
class Command:
    def __init__(self, name, module, group, stamp, link_type, mad8_type, beam_def, par_names, par):
        self.name = name
        self.module = module
        self.group = group
        self.stamp = stamp
        self.link_type = link_type
        self.mad8_type = mad8_type
        self.beam_def = beam_def
        self.par_names = par_names
        self.par = par


class CommandList:
    def __init__(self, name, max_length):
        self.name = name
        self.lst = []
        self.names = []
        self.max_length = max_length

    def add_var(self, var):
        self.lst += [var]


class CommandListList:
    def __init__(self, max_length):
        self.lst = []
        self.names = []
        self.max_length = max_length

    def add_var(self, var):
        self.lst += [var]


# Interface
def new_command(name, nl_length, pl_length, module, group, link, mad_8):
    return None

def delete_command(command):
    return None

def clone_command(command):
    return None

def clone_command_flat(command):
    return None

def find_command(name, command_list):
    return None


# List Commands
def new_command_list(name, max_length):
    return CommandList(name, max_length)

def delete_command_list(command_list):
    return None

def find_command_list(name, command_list_list):
    return None

def grow_command_list(command_list):
    return None

def add_to_command_list(label, command, command_list, flag):
    return None

def new_command_list_list(max_length):
    return CommandListList(max_length)

def add_to_command_list_list(label, command_list, command_list_list):
    return None


# Command Methods
def exec_command():
    return None

def decode_command():
    return None

def dump_command(command): # For Debugging
    return None

def print_command(command):
    return None

def store_command_def(): # Processes command definition
    return None

def make_line(statement):
    return None

def get_stmt(file, supp_flag):
    return None

def get_defined_commands():
    return None

def remove_from_command_list(label, command_list):
    return None

def exec_add_expression(in_cmd):
    return None
