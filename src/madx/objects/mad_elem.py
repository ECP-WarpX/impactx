# TODO: Implement on-demand
class Element:
    def __init__(self, name, def_type, bv, length, defi, parent, stamp, base_type, aperture, multipole):
        self.name = name
        self.def_type = def_type
        self.bv = bv
        self.length = length
        self.defi = defi
        self.parent = parent
        self.stamp = stamp
        self.base_type = base_type
        self.aperture = aperture
        self.multipole = multipole

class ElementList:
    def __init__(self, max_length):
        self.lst = []
        self.max_length = max_length

    def add_var(self, var):
        self.lst += [var]


def make_element(name, parent, command, flag):
    return None

def clone_element(element):
    return None

def delete_element(element):
    return None

def update_element(element, update_command):
    return None

def update_element_children(element, update_command):
    return None


def dump_element(element):
    return None

def export_el_def(element, string, noexpr):
    return None

def export_el_def_8(element, string, noexpr):
    return None


def new_el_list(max_length):
    return ElementList(max_length)

def delete_el_list(el_list):
    return None

def find_element(name, el_list):
    return None

def write_elems(el_list, command_list, file, noexpr):
    return None

def write_elems_8(el_list, command_list, file):
    return None


def new_elem_node(element, occ_cnt):
    return None

def make_elem_node(element, occ_cnt):
    return None

def compound(e_name, occ_cnt):
    return None


def enter_element(in_cmd):
    return None

def element_name(name, l):
    return None

def element_value(node, params):
    return None

def element_vector(element, params, vector):
    return None





