# Variable-width-contouring

Code for computing beads (aka toolpaths, print-paths, print-trajectories) in a slice with varying width for minimizing underfill.
This is the code developed for the publication

> _Variable-width contouring for additive manufacturing_. Samuel Hornus, Tim Kuipers, Olivier Devillers, Monique Teillaud, Jonàs Martínez, Marc Glisse, Sylvain Lazard and Sylvain Lefebvre. ACM Transactions on Graphics 39(4) (Proc. SIGGRAPH 2020) [[link]](https://hal.inria.fr/hal-02568677)

#### Visualization and analysis

Code and scripts for visualizing and analyzing variable-width GCode can be
found [here](https://github.com/BagelOrb/ToolpathVisualizer).

#### Authors of the code

- [Samuel Hornus](https://members.loria.fr/samuel.hornus/) (Inria, France)
- [Tim Kuipers](https://kuipers.weblog.tudelft.nl/) (Ultimaker, TU Delft, Netherlands)

[Salim Perchy](https://members.loria.fr/SPerchy/) and [Pierre Bedell](https://fr.linkedin.com/in/pierre-bedell-667328a0) also contributed to Windows- and plugin-specific code.

## Compilation dependencies

### For compiling the `fill` binary:

- CMAKE
- TCLAP http://tclap.sourceforge.net
- Clipper http://angusj.com/delphi/clipper.php
- If the ``VWC_USE_CGAL_MEDIAL_AXIS`` CMake option is checked, then the code depends on CGAL https://www.cgal.org
- If **both** the ``VWC_USE_BOOST_MEDIAL_AXIS`` and ``VWC_USE_EXTERNAL_BOOST`` CMake options are checked, then the code depends on BOOST
- If the ``VWC_USE_CAIRO_PDF`` CMake option is checked, then the code depends on Cairo https://www.cairographics.org

### For compiling the [IceSL](https://icesl.loria.fr/) plugin:

If the ``VWC_MAKE_ICESL_PLUGIN`` CMake option is checked, then the plugin
`libVariableWidthContouring.so` is compiled.
This needs the `icesl-interface` library sources to be copied (typically as a git repoitory) alongside the `src` directory. (That library sources is **not yet** publicly available.)

## Compiling

It should work alright on macOS and Ubuntu.

	mkdir build
	cd build
	cmake ..
	make

You may have to adjust the CMakeLists.txt.
Please send us your ``pull`` requests!

There is the ``VWC_USE_CAIRO_PDF``CMake option to disable the use of Cairo. Then no PDF is ever output.

If you compile with both ``VWC_USE_BOOST_MEDIAL_AXIS`` and ``VWC_USE_CGAL_MEDIAL_AXIS`` then a command-line option ``-b / --boost`` is added to choose to use BOOST instead of the default CGAL for computing the medial axes.

## Usage of `fill`

``fill -h`` will give you a a lot of options.

The basic usage is ``fill -p input_file -o output_file``.
Without  ``-o`` just the PDF file ``input_file.pdf`` is generated (unless ``--no-pdf``).

### Input format

The ``input_file`` starts with the minimal and maximal allowed bead width.
The rest is a sequence of polygonal closed curves.
Each curve starts with the number of sample points (on one line).
Then there is one sample per line.
Each sample gives ``x, y``.

### Output format

The ``output_file`` is a (textual) sequence of print paths.
Each print path start with the number of samples.
Then there is one sample per line.
Each sample gives ``x, y, radius (half width), tangent_x, tangent_y``.

### Dataset

The ``input/`` directory has some sample input files.

The set of 300 input files in input/dataset/ is the one used in the paper and
originates from the following work:

> _A framework for adaptive width control of dense contour-parallel toolpaths in fused deposition modeling_. Tim Kuipers, Eugeni L. Doubrovski, Jun Wu, and Charlie C. L. Wang. 2020. _Computer Aided Design, vol. 128._ [arXiv](https://arxiv.org/abs/2004.13497), [DOI](http://dx.doi.org/10.1016/j.cad.2020.102907).

### Generating lots of output

You can use this simple shell command:

	for f in `ls ../input/*.txt ../input/dataset/*.txt`; do ./fill -p ${f}; done

This will generate one PDF file and store it alongside each input file. Each page shows the state of the medial axis and the new bead.

## Usage of the IceSL plugin

Copy the file `libVariableWidthContouring.so` into `icesl-plugins/` directory.
