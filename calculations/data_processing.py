import numpy as np

def interpolate(x: float, x_values: np.ndarray, y_values: np.ndarray):
    """
    Interpolates the value of y corresponding to x using linear interpolation.

    :param x: The value at which to interpolate.
    :type x: float
    :param x_values: The sorted array of x-values.
    :type x_values: np.ndarray
    :param y_values: The corresponding array of y-values.
    :type y_values: np.ndarray

    :return: The interpolated y-value at x.
    :rtype: float
    """
    if len(x_values) != len(y_values):
        raise ValueError("x_values and y_values must have the same length")

    if x < x_values[0] or x > x_values[-1]:  # Check if x is within the range of x_values
        raise ValueError(f"x ({x}) is out of bounds. Valid range is [{x_values[0]}, {x_values[-1]}].")

    y = np.interp(x, x_values, y_values)

    return y