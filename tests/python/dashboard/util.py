import os
import subprocess


def wait_for_ready(sb, element_name, timeout):
    # Code adapted from https://github.com/Kitware/trame-client/blob/master/trame_client/utils/testing.py
    """
    Function waits until element_name is present.
    """
    for i in range(timeout):
        print(f"wait_for_ready {i}")
        if sb.is_element_present(element_name):
            sb.sleep(1)
        else:
            print("Ready")
            return


def wait_for_dashboard_ready(process, timeout=60):
    """
    Function waits until the dashboard server is ready by checking the process output.
    """
    for i in range(timeout):
        line = process.stdout.readline()
        if line:
            print(line, end="")
            if "App running at:" in line:
                print("Dashboard is ready!")
                return
    raise Exception("Dashboard did not start correctly.")


def set_input_value(sb, element_id, value):
    """
    Function to clear, update, and trigger a change event on an input field by ID.
    """
    selector = f"#{element_id}"
    sb.clear(selector)

    if not isinstance(value, str):
        value = str(value)

    sb.update_text(selector, value)
    sb.send_keys(selector, "\n")


def start_dashboard():
    """
    Function which starts up impactx-dashboard server.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    working_directory = os.path.join(script_dir, "../../../src/python/impactx")
    working_directory = os.path.normpath(working_directory)

    return subprocess.Popen(
        ["python", "-m", "dashboard", "--server"],
        cwd=working_directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
