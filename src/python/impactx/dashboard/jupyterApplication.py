from .__main__ import application


class JupyterMainApplication:
    """
    Creates specific components from
    Trame Application for Jupyter Notebook.
    """

    def __init__(self):
        self.ui = self.generate_ui()

    def generate_ui(self):
        return application()


"""
from dashboard.jupyterApplication import JupyterMainApplication

# Create new application instance
app = JupyterMainApplication()

# Let's start the server by waiting for its UI to be ready
await app.ui.ready

# Put the UI into the resulting cell
app.ui
"""
