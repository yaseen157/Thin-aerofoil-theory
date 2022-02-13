# Thin-aerofoil-theory Package

Tools for the inviscid analysis of aerofoils with "Thin aerofoil" theory.

## Windows Installation Instructions
> Using information from:<br>
> https://packaging.python.org/en/latest/tutorials/packaging-projects/<br>
> https://packaging.python.org/en/latest/tutorials/installing-packages/<br>
#### Preliminary setup
1. Download the package source files.
2. Open a command prompt inside the package folder.

#### Generating distribution archives
3. If no ".tar.gz" files exist in package/dist/, you need to type this in your command prompt.<br>
cmd > `py -m pip install --upgrade build`<br>
cmd > `py -m build`<br>

#### Installing a Local Package (without venv)
4. Follow these instructions to install the package, replacing path name and version number as appropriate.<br>
cmd > `py --version`<br>
cmd > `py -m pip --version`<br>
cmd > `py -m pip install --upgrade pip setuptools wheel`<br>
cmd > `py -m pip install ./downloads/thinaerofoils-yaseen157-x.x.x.tar.gz`<br>

## Windows Uninstallation Instructions
> Using information from:<br>
> https://pip.pypa.io/en/stable/cli/pip_uninstall/<br>
1. Open a command prompt in the environment thinaerofoils is installed in.<br>
2. Type the following into the command prompt
cmd > `py -m pip uninstall thinaerofoils-yaseen157`<br>

---

## Hello World.
This package currently only supports 4-digit NACA aerofoils:<br>
py > `from thinaerofoils.inviscidanalysis import NACA4digit`<br>
py > `myfoil = NACA4digit(foil="2412")`<br>
py > `myfoil.show2()`<br>
<br>
![NACA2412 demonstrative plot](./docs/images/demo_2412.png)
