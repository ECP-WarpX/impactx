import os
import subprocess
import sys
import time

from selenium.common.exceptions import (
    ElementNotInteractableException,
    JavascriptException,
)


def wait_for_ready(sb, element_name, timeout=10):
    """
    Waits until the specified element is present in the DOM.
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


def set_input_value(sb, element_id, value, timeout=60):
    """
    Function to clear, update, and trigger a change event on an input field by ID.
    Waits until the element is interactable before performing actions.
    """

    selector = f"#{element_id}"
    end_time = time.time() + timeout

    while True:
        try:
            sb.clear(selector)
            sb.update_text(selector, value)
            sb.send_keys(selector, "\n")
            break
        except (ElementNotInteractableException, JavascriptException):
            if time.time() > end_time:
                raise Exception(
                    f"Element {selector} not interactable after {timeout} seconds."
                )


def look_for_text(sb, element_id, text_to_look_for, timeout):
    """
    Function which repeatedly checks for specific text every second until it appears or the timeout is reached.
    """

    start_time = time.time()
    while time.time() - start_time < timeout:
        try:
            if sb.is_text_visible(text_to_look_for, element_id):
                return True
        except Exception as e:
            print(f"Error checking for text: {e}")
        time.sleep(1)

    return False


def start_dashboard():
    """
    Function which starts up impactx-dashboard server.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    working_directory = os.path.join(script_dir, "../../../src/python/impactx")
    working_directory = os.path.normpath(working_directory)

    return subprocess.Popen(
        [sys.executable, "-m", "dashboard", "--server"],
        cwd=working_directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
