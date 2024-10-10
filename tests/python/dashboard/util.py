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


def check_until_visible(sb, selector, timeout=10, interval=1):
    """
    Function which retries checking if an element is visible.
    """

    end_time = time.time() + timeout
    while time.time() < end_time:
        if sb.is_element_visible(selector):
            return True
        time.sleep(interval)
    return False


def find_repo_root():
    """
    Finds the root directory of the repository by looking for a .git directory.
    """
    current_dir = os.path.abspath(os.path.dirname(__file__))
    while True:
        if os.path.isdir(os.path.join(current_dir, ".git")):
            return current_dir
        parent_dir = os.path.dirname(current_dir)
        if parent_dir == current_dir:
            raise Exception("Repository root not found.")
        current_dir = parent_dir


def start_dashboard():
    """
    Function which starts up impactx-dashboard server.
    """
    repo_root = find_repo_root()
    working_directory = os.path.join(repo_root, "src", "python", "impactx", "dashboard")
    working_directory = os.path.normpath(working_directory)

    return subprocess.Popen(
        [sys.executable, "-m", "dashboard", "--server"],
        cwd=working_directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
