import warnings

# Ignore mpld3 warning
warnings.filterwarnings("ignore", message="Blended transforms not yet supported. Zoom behavior may not work as expected.", category=UserWarning, module="mpld3.mplexporter.exporter")

from trame.app import get_server
from trame.ui.vuetify import SinglePageWithDrawerLayout
from trame.widgets import vuetify, matplotlib

import matplotlib.pyplot as plt
import numpy as np

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------
state.alpha_x = 1
state.beta_x = 1
state.epsilon_x = 1

state.alpha_y = 1
state.beta_y = 1
state.epsilon_y = 1

state.alpha_t = 1
state.beta_t = 1
state.epsilon_t = 1

# Update plots on parameter change
@state.change("alpha_x", "beta_x", "epsilon_x")
def on_params_change_x(**kwargs):
    visualizeTwiss.update_plot_x()

@state.change("alpha_y", "beta_y", "epsilon_y")
def on_params_change_y(**kwargs):
    visualizeTwiss.update_plot_y()

@state.change("alpha_t", "beta_t", "epsilon_t")
def on_params_change_t(**kwargs):
    visualizeTwiss.update_plot_t()

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

class visualizeTwiss:
    def draw_phase_space_ellipse(alpha, beta, epsilon, title, xlabel, ylabel, n_points=100):
        """
        Draw a phase space ellipse using Twiss parameters and beam emittance.

        Parameters:
        alpha (float): Twiss parameter alpha
        beta (float): Twiss parameter beta
        epsilon (float): Beam emittance
        n_points (int): Number of points used to draw the ellipse
        """
        # Calculate gamma from alpha and beta
        gamma = (1 + alpha**2) / beta

        # Covariance matrix
        sigma = epsilon * np.array([[beta, -alpha], [-alpha, gamma]])

        # Eigenvalues and eigenvectors
        eigvals, eigvecs = np.linalg.eigh(sigma)

        # Semi-major and semi-minor axes
        a = np.sqrt(eigvals[1])
        b = np.sqrt(eigvals[0])

        # Rotation angle
        theta = np.arctan2(eigvecs[1, 1], eigvecs[0, 1])

        # Generate points for the ellipse
        t = np.linspace(0, 2 * np.pi, n_points)
        cos_t = np.cos(t)
        sin_t = np.sin(t)

        # Generate the ellipse points before rotation
        q0 = a * cos_t
        p0 = b * sin_t

        # Apply the rotation
        q = q0 * np.cos(theta) - p0 * np.sin(theta)
        p = q0 * np.sin(theta) + p0 * np.cos(theta)

        # Create the plot with a smaller size
        fig, ax = plt.subplots(figsize=(3.25, 2.25)) 
        
        ax.plot(q, p, label=f'ε={epsilon}, α={alpha}, β={beta}', linestyle='-', marker='o')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        legend = ax.legend(loc='upper right')
        legend.get_frame().set_facecolor('lightgray') 
        legend.get_frame().set_alpha(0.8)
        legend.get_frame().set_edgecolor('black')
        legend.get_frame().set_linewidth(1.5)
        ax.grid(True)
        ax.axis('equal')

        # Fix axis limits
        # ax.set_xlim([-10, 10])
        # ax.set_ylim([-10, 10])

        fig.subplots_adjust(left=0.2, bottom=0.2, right=.97, top=.95)

        return fig

    def update_plot_x(**kwargs):
        alpha = state.alpha_x
        beta = state.beta_x
        epsilon = state.epsilon_x
        fig = visualizeTwiss.draw_phase_space_ellipse(alpha, beta, epsilon, title='Phase Space Ellipse (x-px)', xlabel='x', ylabel='p_x')
        ctrl.matplotlib_figure_update_x(fig)
        plt.close(fig)

    def update_plot_y(**kwargs):
        alpha = state.alpha_y
        beta = state.beta_y
        epsilon = state.epsilon_y
        fig = visualizeTwiss.draw_phase_space_ellipse(alpha, beta, epsilon, title='Phase Space Ellipse (y-py)', xlabel='y', ylabel='p_y')
        ctrl.matplotlib_figure_update_y(fig)
        plt.close(fig)

    def update_plot_t(**kwargs):
        alpha = state.alpha_t
        beta = state.beta_t
        epsilon = state.epsilon_t
        fig = visualizeTwiss.draw_phase_space_ellipse(alpha, beta, epsilon, title='Phase Space Ellipse (t-pt)', xlabel='t', ylabel='p_t')
        ctrl.matplotlib_figure_update_t(fig)
        plt.close(fig)

    def card_x():
        with vuetify.VCard(style="width: 340px; height: 300px"):
            with vuetify.VCardText():
                with vuetify.VRow(classes="pl-1"):
                    vuetify.VIcon("mdi-chart-arc")
                    vuetify.VCardTitle("Phase Space Ellipse (x-px)", classes="text-subtitle-2", style="color: black;")
            vuetify.VDivider()
            matplotlib_figure = matplotlib.Figure()
            ctrl.matplotlib_figure_update_x = matplotlib_figure.update

    def card_y():
        with vuetify.VCard(style="width: 340px; height: 300px"):
            with vuetify.VCardText():
                with vuetify.VRow(classes="pl-1"):
                    vuetify.VIcon("mdi-chart-arc")
                    vuetify.VCardTitle("Phase Space Ellipse (y-py)", classes="text-subtitle-2", style="color: black;")
            vuetify.VDivider()
            matplotlib_figure = matplotlib.Figure()
            ctrl.matplotlib_figure_update_y = matplotlib_figure.update

    def card_t():
        with vuetify.VCard(style="width: 340px; height: 300px"):
            with vuetify.VCardText():
                with vuetify.VRow(classes="pl-1"):
                    vuetify.VIcon("mdi-chart-arc")
                    vuetify.VCardTitle("Phase Space Ellipse (t-pt)", classes="text-subtitle-2", style="color: black;")
            vuetify.VDivider()
            matplotlib_figure = matplotlib.Figure()
            ctrl.matplotlib_figure_update_t = matplotlib_figure.update

# -----------------------------------------------------------------------------
# Layout
# -----------------------------------------------------------------------------

# with SinglePageWithDrawerLayout(server) as layout:

#     with layout.content:
#         with vuetify.VContainer(fluid=True):
#             with vuetify.VRow():
#                 with vuetify.VCol(cols=4):
#                     visualizeTwiss.card_x()
#                 with vuetify.VCol(cols=4):
#                     visualizeTwiss.card_y()
#                 with vuetify.VCol(cols=4):
#                     visualizeTwiss.card_t()

