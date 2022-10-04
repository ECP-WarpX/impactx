# MadX Python Interpreter (10/4/22)

Author: Joey Zhu (https://github.com/np-eazy)

This Python interpreter was written as a subset of MadX's C interpreter's full functionalities to suit the needs of ImpactX.
Implementation contains a lot of boilerplate, but this is to maintain a reliable, verbatim approach for transcribing a large
volume of code.

### Supported objects
- Variables
- Elements
- Macros (Lines)

### Removed/Changed Functionalities
- stdin stdout char buffers are replaced with Python file read
- Removed Fortran90 interface

### Other dev notes
- Execution begins at mad_main.py