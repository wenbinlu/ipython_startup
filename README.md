# ipython_startup
pre-loaded startup files for ipython, including definitions for useful physical constants, and pre-defined functions

Step 1: place these files in ~/.ipython/profile_default/startup/ directory

Step 2 (optional): include the following line in .bashrc (or equivalent file) to allow for inline plots

alias ipython="ipython --InteractiveShellApp.exec_lines='%matplotlib Qt5'"



# ------ old README

This is the IPython startup directory

.py and .ipy files in this directory will be run *prior* to any code or files specified
via the exec_lines or exec_files configurables whenever you load this profile.

Files will be run in lexicographical order, so you can control the execution order of files
with a prefix, e.g.::

    00-first.py
    50-middle.py
    99-last.ipy