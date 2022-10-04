# include <signal.h>
# include "madx.h"

from mad_var import *
from mad_elem import *
from mad_macro import *
from mad_cmd import *


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
    # TODO: Implement normal madx_start

    mad_init_c(gvar)

    time = None
    # /*  setbuf(stdout,(char *)0); */ /* no buffering - for debugging */
    # time(&start_time); /* initialize timing */
    # tm = localtime(&start_time); /* split system time */
    # last_time = start_time;

    # // compute padding of OSTYPE
    # const char *pad[] = { "", " ", "  ", "    " };
    # const int pad_sz = sizeof pad/sizeof *pad;
    # int pad_idx = strlen("Windows")-strlen(version_ostype);
    # if (pad_idx >= pad_sz) pad_idx = pad_sz-1;

    #  printf("\n  ++++++++++++++++++++++++++++++++++++++++++++\n");
    # printf("  +     MAD-X %s  (%s bit, %s) %s    +\n", version_name, version_arch, version_ostype, pad[pad_idx]);
    # printf("  + Support: mad@cern.ch, http://cern.ch/mad +\n");
    # printf("  + Release   date: %s               +\n", version_date);
    # printf("  + Execution date: %04d.%02d.%02d %02d:%02d:%02d      +\n",
    #        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday,
    #                                            tm->tm_hour, tm->tm_min, tm->tm_sec);
    # printf("  ++++++++++++++++++++++++++++++++++++++++++++\n");
    # if (interactive) putchar('\n');
    return None


def madx_input(top):
    # TODO: Implement normal madx_input
    # TODO: Implement dynamic metadata lists
    # {
    #     fputs("This is where STDOUT is", stdout);
    # if (interactive && in->curr == 0) fputs("X:> ", stdout);
    # if (return_flag || get_stmt(in->input_files[in->curr], 0) == 0)
    # {
    # if (in->input_files[in->curr] != stdin) {
    # fclose(in->input_files[in->curr]);
    # in->input_files[in->curr] = 0;
    # }
    # if (in->curr == 0)
    # return;
    # else
    #     in->curr -= 1;
    # return_flag = 0;
    # if (in->curr == top) return;
    # }
    # else
    # {
    # stolower_nq(in->buffers[in->curr]->c_a->c);
    # pro_input(in->buffers[in->curr]->c_a->c);
    # if (stop_flag)  return;
    # }
    # }
    return None


def madx_finish():
    # TODO: Implement normal madx_finish
    # TODO: Implement MadX-ImpactX Interface by dumping lists
    # int warn_numb, warn_numbf, nwarnings;
    #
    # /* should work with Lahey on windows 24.03.2004 */
    #
    # match2_delete_expressions();
    # match2_delete_arrays();
    #
    # if (final_message == 0)
    #     {
    #         final_message = 1;
    # if (plots_made)
    # {
    #     gxterm_();
    # }
    # mad_err_getwarn(&warn_numb, &warn_numbf);
    # nwarnings = warn_numb + warn_numbf;
    # printf("\n  Number of warnings: %d\n",nwarnings);
    # if (nwarnings > 0)
    # {
    #     printf("%d in C and %d in Fortran\n",warn_numb,warn_numbf);
    # }
    # if (get_option("trace")) time_stamp("end");
    #
    # // printf("poly-mul= %llu\n", poly_mul_count);
    #
    # printf("\n"
    # "  ++++++++++++++++++++++++++++++++++++++++++++\n");
    # printf("  +          MAD-X finished normally         +\n");
    # printf("  ++++++++++++++++++++++++++++++++++++++++++++\n");
    # }
    return None
