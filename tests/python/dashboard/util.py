def wait_for_ready(sb, timeout=60):
    for i in range(timeout):
        print(f"wait_for_ready {i}")
        if sb.is_element_present(".trame__loader"):
            sb.sleep(1)
        else:
            print("Ready")
            return


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
