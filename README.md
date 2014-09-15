streamer_1d
===========

This is code for 1d streamer simulations with a particle of fluid code.

The code uses git submodules, so for the first compilation you should execute
these commands:
$ git submodule init
$ git submodule update
$ make

After that, you can compile with just
$ make

Running (sequential):
$ ./streamer_1d my_config_file.txt

You can also specify multiple configuration files, like:
$ ./streamer_1d cfg_base.txt cfg_1.txt

In each configuration file, you can specify a different "sim_name" variable.
These names will be appended to each other.
