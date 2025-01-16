import numpy as np
import matplotlib.pyplot as plt
import os

def plot_pressure_data(pressure_data: np.ndarray, title: str, color: str,
                       size: int, file_path: str, file_name: str, show_plot: bool):
    """
    Plot the extracted pressure data (x, y, z, Pressure) in a 3D scatter plot with configurable color.

    :param pressure_data: Array of shape (n, 4), where each row contains the values [x, y, z, Pressure].
    :type pressure_data: np.ndarray
    :param title: The title to display on the plot.
    :type title: str
    :param color: The color or colormap to use for the scatter plot.
    :type color: str
    :param size: The size of the markers in the scatter plot.
    :type size: int
    :param file_path: Path where the plot will be saved.
    :type file_path: str
    :param file_name: Name of the image file (e.g., "image.png").
    :type file_name: str
    :param show_plot: Whether to display the plot or not (default is False).
    :type show_plot: bool
    """
    # Extract x, y, z, Pex values
    x_vals = pressure_data[:, 0]
    y_vals = pressure_data[:, 1]
    z_vals = pressure_data[:, 2]
    Pex_vals = pressure_data[:, 3]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(x_vals, y_vals, z_vals, c=Pex_vals, cmap=color, s=size)

    # Color bar for pressure
    cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.2, fraction=0.05)
    cbar.set_label(r'$P_\mathrm{ex}\ (\mathrm{kN/m}^2)$')

    # Set axis labels
    ax.set_xlabel(r'$x\ (\mathrm{m})$')
    ax.set_ylabel(r'$y\ (\mathrm{m})$')
    ax.set_zlabel(r'$z\ (\mathrm{m})$')

    # Set plot title
    plt.title(title)

    # Get the limits of each axis
    x_limits = ax.get_xlim()
    y_limits = ax.get_ylim()
    z_limits = ax.get_zlim()

    # Find the range of each axis
    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    # Find the largest range
    max_range = max(x_range, y_range, z_range)

    # Center the axes around the midpoint and set equal ranges
    x_mid = np.mean(x_limits)
    y_mid = np.mean(y_limits)
    z_mid = np.mean(z_limits)

    ax.set_xlim([x_mid - max_range / 2, x_mid + max_range / 2])
    ax.set_ylim([y_mid - max_range / 2, y_mid + max_range / 2])
    ax.set_zlim([z_mid - max_range / 2, z_mid + max_range / 2])

    if file_path and file_name:

        full_file_path = os.path.join(file_path.rstrip('/'), file_name) # Construct the full file path
        plt.savefig(full_file_path, dpi=300, bbox_inches='tight')

    if show_plot:
        plt.show()
    else:
        plt.close(fig)