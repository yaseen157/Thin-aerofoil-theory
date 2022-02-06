# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 19:57:51 2022

@author: yr3g17
"""
from matplotlib import pyplot as plt
import numpy as np
from scipy import integrate as sig

from NACAfoils import parse_NACA4digit


class NACA4digit:

    def __init__(self, foil: str, xres: int = None):
        """Instantiate a class object with a NACA foil code"""
        if xres is None:
            xres = 100

        coords, f_z, f_dzdx = parse_NACA4digit(foil, xres)

        self.foil = foil
        self.xres = xres
        self.coords = coords
        self.f_z = f_z
        self.f_dzdx = f_dzdx

        def f_0lift_integrand(x):
            theta0 = np.arccos(1 - 2 * x)
            result = f_dzdx(x) * (1 - np.cos(theta0)) * 2 / np.sin(theta0)
            return result
        self.zeroliftAoA_rad = 1 / np.pi * sig.quad(f_0lift_integrand, 0, 1)[0]

        return None

    def __repr__(self):
        foil = self.foil
        xres = self.xres
        return f"{type(self).__name__}({foil=}, {xres=})"

    def show(self) -> None:
        """Call to create a figure showing off the aerofoil's profile"""
        xs = np.linspace(0, 1, self.xres)

        fig, ax = plt.subplots(1, subplot_kw={"aspect": 1}, dpi=150)
        ax.set_title(f"NACA {self.foil}")
        ax.plot(*self.coords, "black")
        ax.plot(xs, self.f_z(xs), "r")
        ax.plot((0, 1), (0, 0), ls="-.", c="black")
        ax.set_xlabel("x/c")
        ax.set_ylabel("y/c")
        plt.show()
        return None

    def show2(self) -> None:
        """Create a tile of plots summarising the aerofoil's characteristics"""
        fig, axd = plt.subplot_mosaic(
            [['foil', 'foil'], ['lift', 'moment']],
            constrained_layout=True, dpi=100, figsize=(8*0.9, 5*0.9))
        fig.suptitle(f"NACA {self.foil} Aerofoil")

        # Plot the aerofoil profile
        xs = np.linspace(0, 1, self.xres)
        axd["foil"].plot(*self.coords, c="black")
        axd["foil"].fill_between(*self.coords, color="black", alpha=0.1)
        axd["foil"].plot(xs, self.f_z(xs), c="red", label="z=z(x) Camber line")
        axd["foil"].plot((0, 1), (0, 0), ls="-.", c="black")
        axd["foil"].set_xlabel("x/c")
        axd["foil"].set_ylabel("y/c")
        axd["foil"].set_aspect(1)
        hspace = xs[0::self.xres-1] - axd["foil"].get_xlim()
        vspace = np.array([min(self.coords[1]), max(self.coords[1])]) - hspace
        axd["foil"].set_ylim(*vspace)
        axd["foil"].legend()

        # Plot the lift-curve
        alphas = np.linspace(np.radians(-20), np.radians(20), 2)
        clift = [self.f_Cl(x) for x in alphas]
        zeroliftstr = "$\\alpha_{L=0}$="
        zeroliftstr += f"{np.degrees(self.zeroliftAoA_rad):.2f}"
        axd["lift"].plot(np.degrees(alphas), clift, "r")
        axd["lift"].plot(
            np.degrees(alphas), alphas * 2 * np.pi,
            ls="-.", c="black", alpha=0.8, label="d$C_l$/d$\\alpha$=$2\\pi$")
        axd["lift"].plot(
            np.degrees(self.zeroliftAoA_rad), 0,
            "x", c="blue", label=zeroliftstr)
        axd["lift"].set_title("$C_l$ vs $\\alpha$")
        axd["lift"].set_xlabel("Angle of Attack ($\\alpha$) [deg]")
        axd["lift"].set_ylabel("$C_l$")
        axd["lift"].grid()
        axd["lift"].legend()

        # Plot the moment-curve
        alphas = np.radians(1e-2 * np.array([-1, 1]))
        xs = np.linspace(0, 1, self.xres)
        cmoms0 = self.f_cm(alpha_rad=alphas[0], xs=xs)
        cmoms1 = self.f_cm(alpha_rad=alphas[-1], xs=xs)
        cmas = (cmoms1 - cmoms0) / (alphas[-1] - alphas[0])
        dcmasdx = (cmas[-1] - cmas[0]) / (xs[-1] - xs[0])
        zerocmax = (cmas[0] - dcmasdx * xs[0]) / (-dcmasdx)
        staticcmom = self.f_cm(alpha_rad=0, xs=zerocmax)
        staticcmomstr = f"$c_m={staticcmom:.6f}$"
        zerocmastr = "(x/c)$_{c_{m,\\alpha}=0}$=" + f"{zerocmax:.2f}"
        axd["moment"].plot(xs, cmas, "r")
        axd["moment"].plot(xs, np.zeros(len(xs)), ls="-.", c="black", alpha=.8)
        axd["moment"].plot(zerocmax, 0, "x", c="blue", label=zerocmastr)
        axd["moment"].set_title("$c_{m,\\alpha}$ vs x/c")
        axd["moment"].set_xlabel("x/c")
        axd["moment"].set_ylabel("d$c_m$ / d$\\alpha}$")
        axd["moment"].grid()
        axd["moment"].legend(title=staticcmomstr)

        plt.show()

        return None

    def f_An_fundamental(self, alpha_rad: float, nmax: int = None):
        """Return cos-series expansion of the fundamental eq. for thin foils"""

        if nmax is None:
            nmax = 100

        assert isinstance(nmax, int), "Argument nmax is of wrong type"
        assert nmax > 0, "nmax is too low"

        def f_A0_integrand(x):
            theta0 = np.arccos(1 - 2 * x)
            result = self.f_dzdx(x) * 2 / np.sin(theta0)
            return result

        def f_An_integrand(x, n):
            theta0 = np.arccos(1 - 2 * x)
            result = self.f_dzdx(x) * np.cos(n * theta0) * 2 / np.sin(theta0)
            return result

        An = np.zeros(nmax+1)
        An[0] = alpha_rad - 1 / np.pi * sig.quad(f_A0_integrand, 0, 1)[0]
        for i in range(1, len(An)):
            An[i] = 2 / np.pi
            An[i] *= sig.quad(f_An_integrand, 0, 1, args=i, limit=100)[0]
        return An

    def f_Cl(self, alpha_rad: float) -> float:
        """Return the Coefficient of lift at some Angle of Attack"""
        An = self.f_An_fundamental(alpha_rad, nmax=1)

        Cl = 2 * np.pi * (An[0] + An[1] / 2)
        return Cl

    def f_cm(self, alpha_rad: float, xs: float = 0.25) -> float:
        """Return the moment coefficient of the aerofoil at some position"""
        An = self.f_An_fundamental(alpha_rad, nmax=2)

        Cm_LE = -np.pi / 2 * (An[0] + An[1] - An[2] / 2)
        Cl = 2 * np.pi * (An[0] + An[1] / 2)
        Cm = Cm_LE + xs * Cl
        return Cm

    def f_gamma(self, xs, vfreestream_mps, alpha_rad):
        """Return vortex sheet strength along chordline, at position x"""
        An = self.f_An_fundamental(alpha_rad)

        theta = np.arccos(1 - 2 * xs)
        gamma = An[0] * (1 + np.cos(theta)) / np.sin(theta)
        gamma += np.sum([An[i] * np.sin(i * theta) for i in range(1, len(An))])
        gamma *= 2 * vfreestream_mps
        return gamma

    def f_Gamma(self, vfreestream_mps, alpha_rad):
        """Return source strength of vortex sheet"""
        An = self.f_An_fundamental(alpha_rad, nmax=1)

        Gamma = np.pi * vfreestream_mps * (An[0] + An[1] / 2)
        return Gamma


if __name__ == "__main__":

    myfoil = NACA4digit(foil="2412")
    myfoil.show2()
