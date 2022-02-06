# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 12:37:26 2022

@author: yr3g17
"""
from thinaerofoils.inviscidanalysis import NACA4digit


def demo_2412():
    myfoil = NACA4digit(foil="2412")
    myfoil.show2()


if __name__ == "__main__":
    demo_2412()
