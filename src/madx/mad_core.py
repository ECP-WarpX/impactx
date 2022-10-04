# include <signal.h>
# include "madx.h"

from src.madx.objects.mad_var import *
from src.madx.objects.mad_elem import *
from src.madx.objects.mad_macro import *
from src.madx.objects.mad_cmd import *

from src.madx.utils.mad_str import *


def mad_init_c(gvar):
    # TODO: Implement subset of mad_init_c
    # allocate variable ptr
    ione = 1

    # make stdout unbuffered for proper I/O synchronization with fortran
    # setvbuf(stdout, 0, _IONBF, 0);
    #
    # init55(123456789);          /* random generator */
    # if (watch_flag == 1)  debug_file = fopen("madx.debug", "w");
    # else if (watch_flag == 2)  debug_file = stdout;
    # if (stamp_flag == 1)  stamp_file = fopen("madx.stamp", "w");
    # else if (stamp_flag == 2)  stamp_file = stdout;

    # in = new_in_buff_list(100); /* list of input buffers, dynamic */

    # // quick and dirty fix to accept jobs from filename
    # in->input_files[0] = mad_argc > 1 ? fopen(mad_argv[1], "r") : stdin;
    # if (!in->input_files[0]) {
    # mad_error("invalid input filename ", " %s", mad_argv[1]);
    # in->input_files[0] = stdin;
    # }
    # interactive = intrac();
    # prt_file = stdout;
    # pro = new_in_buff_list(100); /* list of process buffers, dynamic */
    # pro->buffers[0] = new_in_buffer(IN_BUFF_SIZE);
    # pro->curr = 1;
    # c_dum = new_char_array(AUX_LG);
    # c_join = new_char_array(AUX_LG);
    # work = new_char_array(AUX_LG);
    # l_wrk = new_char_array(AUX_LG);

    # Physical Constants:
    # deco_init();
    # get_defined_constants();
    # get_defined_commands(command_def);
    # get_defined_commands(element_def);
    # get_sxf_names();
    #
    # pi = get_variable("pi");
    # twopi = two * pi;
    # degrad = 180 / pi;
    # raddeg = pi / 180;
    # e = get_variable("e");
    # clight = get_variable("clight");
    # hbar = get_variable("hbar");

    # Global Variable Initialization:
    # aux_buff = new_char_array(AUX_LG);  /* dynamic temporary buffer */
    gvar.variable_list = new_var_list(2000)  # dynamic list of variables
    gvar.beam_list = new_command_list("beam_list", 10)  # dynamic beam list
    gvar.table_deselect = new_command_list_list(10)  # dynamic
    gvar.table_select = new_command_list_list(10)  # dynamic
    gvar.defined_commands = new_command_list("defined_commands", 100)  # dynamic
    gvar.stored_commands = new_command_list("stored_commands", 500)  # dynamic
    gvar.line_list = new_macro_list(100)  # dynamic
    gvar.macro_list = new_macro_list(100)  # dynamic
    gvar.base_type_list = new_el_list(60)  # dynamic
    gvar.element_list = new_el_list(20000)  # dynamic
    # buffered_cmds = new_in_cmd_list(10000); /* dynamic */
    # sequences = new_sequence_list(20); /* dynamic */
    # match_sequs = new_sequence_list(2);
    # selected_ranges = new_node_list(10000); /* dynamic */
    gvar.selected_elements = new_el_list(10000)  # dynamic
    # tmp_p_array = new_char_p_array(1000); /* dynamic */
    # tmp_l_array = new_char_p_array(1000); /* dynamic */
    # sxf_list = new_name_list("sxf_list", 50); /* dynamic */
    # all_table_lists = new_table_list_list(10); /* dynamic */

    # # Initialize needed work for commands
    var = new_variable("twiss_tol", 1.e-6, 1, 1, None, None)
    add_to_var_list(var, gvar.variable_list, 1)
    # title = permbuff("no-title");
    # set_defaults("option");
    # if (interactive) { int false = 0; set_option("echo", &false); }
    # set_defaults("beam");
    # add_to_command_list("default_beam", current_beam, beam_list, 0);
    # set_defaults("set");
    # set_defaults("setplot");
    # set_defaults("threader");
    # table_register = new_table_list(10); /* dynamic */
    # beta0_list = new_command_list("beta0_list", 10); /* dynamic */
    # savebeta_list = new_command_list("savebeta_list", 10); /* dynamic */
    # seqedit_select = /* dynamic - for "select seqedit" commands */
    #     new_command_list("seqedit_select", 10);
    # error_select = /* dynamic - for "select error" commands */
    #     new_command_list("error-select", 10);
    # save_select = /* dynamic - for "select save" commands */
    #     new_command_list("save_select", 10);
    # slice_select = /* dynamic - for "select makethin" commands */
    #     new_command_list("slice_select", 10);
    # sector_select = /* dynamic - for "select sectormap" commands */
    #     new_command_list("sector_select", 10);
    # interp_select = /* dynamic - for "select sectormap" commands */
    #     new_command_list("interp_select", 10);
    # s_range = new_int_array(10); /* dynamic */
    # e_range = new_int_array(10); /* dynamic */
    # sd_range = new_int_array(10); /* dynamic */
    # ed_range = new_int_array(10); /* dynamic */
    # zero_double(orbit0, 6);
    # zero_double(disp0, 6);
    # zero_double(guess_orbit,6);
    # zero_double(oneturnmat, 36);
    # set_option("twiss_print", &ione);
    return None


def madx_start(gvar):

    mad_init_c(gvar)

    return None


def madx_input(filename, gvar):

    # The original MadX file uses charbuffers, but in the Python version we want to keep
    # implementation simple with filenames and char lists, and whatever equivalent looping
    # logic. We assume that given ImpactX's use case, there is only one case for the input:
    # That is, anything and everything is in file_content, so we don't have to do any buffer logic.
    file_content = []
    with open(filename, 'r') as f:
        file_content += [f.read()]
    i = 0
    in_stop = False
    while not in_stop:
        # TODO: Implement dynamic metadata lists
        # We don't need line numbers per se, just the order of buffer consumption is enough for metadata.
        stolower_nq(file_content, i)

        # TODO: Figure out exactly where pro_input is defined. This is the entry point into the interpreter logic but isn't clearly defined in the source code.
        # pro_input(file_content, i)
        if gvar.stop_flag:
            return

        in_stop = i < len(file_content)
        i += 1  # whenever a buffer char is popped for interpreter

    return None


def madx_finish(gvar):
    for variable in gvar.variable_list:
        # TODO: Convert variable into ImpactX known type and add to other lists
        variable = variable

    for element in gvar.el_list:
        # TODO: Convert element into ImpactX known type and add to other lists
        element = element

    for line in gvar.line_list:
        # TODO: Convert line into ImpactX known type and add to other lists
        line = line

    return None
