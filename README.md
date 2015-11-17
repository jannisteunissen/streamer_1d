Particle & fluid simulations of 1D streamers
===========

This is code for 1d streamer simulations with a particle or fluid model.

Getting the code:
    $ git clone https://github.com/jannisteunissen/streamer_1d.git

Compiling the code:
   $ cd streamer_1d
   $ make

Running the code:
    $ ./streamer_1d cfg_example_1.txt
    $ ./streamer_1d cfg_example_2.txt

You can also specify multiple configuration files, like:
    $ ./streamer_1d cfg_example_1.txt my_other_options.txt

The later files will then override earlier options. In each configuration file,
you can specify a different "sim_name" variable. These names will be appended to
each other for the output.
