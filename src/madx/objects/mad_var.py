# TODO: Implement on-demand
class Constant:
    def __init__(self, name, expr, value, stamp):
        self.name = name
        self.expr = expr
        self.value = value
        self.stamp = stamp


class Variable:
    def __init__(self, name, val, val_type, type, expression, string):
        self.name = name
        self.val = val
        self.val_type = val_type
        self.type = type
        self.expression = expression
        self.string = string


class VariableList:
    def __init__(self, max_length):
        self.lst = []
        self.max_length = max_length

    def add_var(self, var):
        self.lst += [var]


def enter_variable(in_cmd):
    return None


def new_variable(name, val, val_type, typ, expression, string):
    return Variable(name, val, val_type, typ, expression, string)


def variable_value(var):
    return None


def new_var_list(max_length):
    return VariableList(max_length)


def clone_var_list(var_list):
    return None


def delete_var_list(var_list):
    return None


def add_to_var_list(var, var_list, flag):
    var_list.add_var(var)


def find_variable(name, var_list):
    return None


def make_string_variable(string):
    return None


def write_vars(var_list, command_list, file, noexpr):
    return None


def write_vars_8(var_list, command_list, file):
    return None


def set_variable(name, value):
    return None


def set_stringvar(name, string):
    return None


def get_variable(name):
    return None


def get_varstring(name):
    return None


def get_defined_constants():
    return None
