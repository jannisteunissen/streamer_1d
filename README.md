streamer_1d
==

This is code for 1D discharge simulations with a particle-in-cell or fluid
model.

Getting the code:

    $ git clone https://github.com/jannisteunissen/streamer_1d.git

Compiling the code:

    $ cd streamer_1d
    $ make

Running the code:

    $ ./streamer_1d example_N2.cfg
    $ ./streamer_1d example_Ar.cfg

Options are specified in the configuration files that are passed as argument.
You can also specify multiple configuration files, like:

    $ ./streamer_1d a.cfg b.cfg

The later files will then override earlier options. Individual parameters can be
also be passed:

    $ ./streamer_1d example_N2.cfg -end%time=1e-9

Program options
==

Please have a look at the provided example configuration files for all the settings.

Program output
==

Output files are stored in the "output" directory. The file names consists of several parts:

    [output]_[model_name]_[number].txt

So if `output%filename` is "output/test", and you use a particle model, then the first output will be

    test_particle_000001.txt

For the fluid model with the local field approximation "fluid" is added. The files contain a header that describes what they store. For fluid simulations, the columns are:

1. position (m)
2. electric field (V/m)
3. electron density (1/m3)
4. positive ion density (1/m3)
5. potential (V)

To plot the electron density with **gnuplot**:

    gnuplot
    plot "output/example1_fluid_000001.txt" u 1:3

