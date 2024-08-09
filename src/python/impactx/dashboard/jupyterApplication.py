"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

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
