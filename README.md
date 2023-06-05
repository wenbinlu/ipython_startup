# ipython_startup
pre-loaded startup files for ipython, including definitions for useful physical constants, and pre-defined functions

Step 1: place these files in ~/.ipython/profile_default/startup/ directory

Step 2 (optional): include the following line in .bashrc (or equivalent file)
alias ipython="ipython --InteractiveShellApp.exec_lines='%matplotlib Qt5'"