comatMOR
========

comatMOR is a software enabling reduced basis computations in pyMOR (see http://pymor.org/) with input matrices in the .mat format originating by certain specific models. It is developed at the University of Muenster with the Python programming language.

NOTE comatMOR is in a very early development stage and therefore has only rudimentary as well as very specific applications. 

License
-------

Copyright (C) 2015 Falk Meyer

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

The following files contain modified source code originating from other open source software projects:

- src/IRT/timestepping.py (pyMOR)

Installation
------------

The comatMOR package is available on github via

https://github.com/fameyer/comatmor

Download it there or use git to keep on track with recent changes. You have to add comatMOR to your python path, for example by including a comatmor.pth file containing the path to the comatmor/src folder.

The comatMOR module depends on pyMOR and its dependencies, so be sure to have pyMOR installed correctly. For an instruction to install pyMOR see http://pymor.org/. If you intend to install pyMOR on Windows use the unofficial pyMOR for Windows version under https://github.com/fameyer/pymorWin.  

Be sure to state the correct python interpreter of your system in the .m Files you want to call later. An automatic installation script is so far not provided.  

Usage
-----

See the HOWTO.txt files in the problem specific subfolders of src/comatmor/ to get appropriate descriptions of how to use this software.

Contact
-------

Do you have any questions concerning comatMOR? Then feel free to write to:

falk.meyer@wwu.de
