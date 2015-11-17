Streamer_1d
==

This is code for 1d streamer simulations with a particle-in-cell or fluid model.
Currently, fluid models can use the local field approximation or the local
energy approximation.

Getting the code:

    $ git clone https://github.com/jannisteunissen/streamer_1d.git

Compiling the code:

    $ cd streamer_1d
    $ make

Running the code:

    $ ./streamer_1d cfg_example_1.txt
    $ ./streamer_1d cfg_example_2.txt

Options are specified in the configuration files that are passed as argument.
You can also specify multiple configuration files, like:

    $ ./streamer_1d cfg_example_1.txt my_other_options.txt

The later files will then override earlier options. In each configuration file,
you can specify a different "sim_name" variable. These names will be appended to
each other for the output.

Program options
==

Please have a look at the provided example configuration files.

Program output
==

Output files are stored in the "output" directory. The filename consists of several parts:

    [sim_name]_[model_name]_[number].txt

So if sim_name is "test", and you use a particle model, then the first output will be

    test_p_1.txt

For the fluid model with the local field approximation "fl" is added, and for
the fluid model with the energy equation "fl_ee" is added. The columns in these files specify the:

1. position (m)
2. electric field (V/m)
3. electron density (1/m3)
4. positive ion density (1/m3)
5. negative ion density (1/m3)
6. energy density (eV/m3)
7. mean energy (eV)

(This information is also listed in the header of the output files.)

For example, to plot the electron density with **gnuplot**:

    gnuplot
    plot "test_part_1.txt" u 1:3


Todo
==

* Add example transport data for energy equation model
