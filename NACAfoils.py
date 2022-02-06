# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 07:51:31 2022

@author: yr3g17
"""
import numpy as np
import re


def parse_NACA4digit(foilcode: str, xres: int = None) -> tuple:
    """
    Return the coordinates of a NACA 4-digit foil, and a camber line function

    Parameters
    ----------
    foilcode : str
        The four-digit NACA aerofoil code, e.g. '0012', '2412'.
    xres : int, optional
        The number of points describing each of upper and lower surfaces. The
        default is 100.

    Raises
    ------
    ValueError
        Raises if an invalid four-digit NACA code was specified.

    Returns
    -------
    tuple
        A two-element tuple containing a tuple of surface coordinates, and a
        function that returns the position of camber at some non-dimensional
        position along the chord x/c.

    """
    if xres is None:
        xres = 100

    # Parse the 4 digit NACA code
    try:
        assert len(re.findall("[0-9]", foilcode)) == 4
        float(foilcode)
    except (AssertionError, ValueError):
        raise ValueError(f"Invalid 4-digit NACA code '{foilcode}'")
    else:
        digits = [int(x) for x in re.findall("[0-9]", foilcode)]

    # Decode the aerofoil code
    maxcamber = digits[0] / 100
    loccamber = digits[1] / 10 if digits[1] > 0 else 0.3  # avoid Zerodivision
    thickness = (10 * digits[2] + digits[3]) / 100

    # Find the camber line of the aerofoil
    def f_z(x):
        """Return camber line position for non-dimensional chord position"""
        commonterm = 2 * loccamber * x - x**2
        z = np.where(
            x <= loccamber,
            maxcamber * loccamber**-2 * commonterm,
            maxcamber * (1 - loccamber)**-2 * (1 - 2 * loccamber + commonterm))
        return z

    # Find the gradient of the camber line of the aerofoil
    def f_dzdx(x):
        """Return gradient of camber line for non-dimensional chord position"""
        dzdx = 2 * maxcamber * (loccamber - x) * np.where(
            x <= loccamber, loccamber**-2, (1-loccamber)**-2)
        return dzdx

    # Arctan of the camber gradient gives angle at which half-thickness applies
    def f_ztheta(x):
        """Return angle of half-thickness application to the camber"""
        theta_rad = np.arctan(f_dzdx(x))
        return theta_rad

    # The polynomial formula for half-thickness (surface to camber line)
    polyco = [0.2969, -0.1260, -0.3516, 0.2843, -0.1015]

    def f_yt(x):
        """Return the half-thickness of the aerofoil at some chord position"""
        y_halfthickness = 0
        for i in range(5):
            y_halfthickness += polyco[i]*x**0.5 if i == 0 else polyco[i]*x**i
        y_halfthickness *= 5 * thickness
        return y_halfthickness

    # Begin coordinate system from upper TE to LE to lower TE
    x_upper = np.linspace(1, 0, xres)
    x_lower = np.flip(x_upper)[1:]

    # Calculate x-coordinates for upper and lower surfaces
    yt_upper = f_yt(x_upper)
    yt_lower = f_yt(x_lower)
    xs_upper = x_upper - yt_upper * np.sin(f_ztheta(x_upper))
    xs_lower = x_lower + yt_lower * np.sin(f_ztheta(x_lower))

    # Upper and lower camber position (identical, but reversed orders)
    z_upper = f_z(x_upper)
    z_lower = f_z(x_lower)

    # Calculate y-coordinates for upper and lower surfaces
    ys_upper = z_upper + yt_upper * np.cos(f_ztheta(x_upper))
    ys_lower = z_lower - yt_lower * np.cos(f_ztheta(x_lower))

    xs = np.concatenate((xs_upper, xs_lower))
    ys = np.concatenate((ys_upper, ys_lower))

    return np.array([xs, ys]), f_z, f_dzdx


if __name__ == "__main__":

    xy, f_z, f_dz = parse_NACA4digit(foilcode="2412")
