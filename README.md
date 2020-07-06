# Variable-width-contouring

Code for computing beads (aka toolpaths, print-paths, print-trajectories) in a slice with varying width for minimizing underfill.
This is the code developed for the publication

> _Variable-width contouring for additive manufacturing_. Samuel Hornus, Tim Kuipers, Olivier Devillers, Monique Teillaud, Jonàs Martínez, Marc Glisse, Sylvain Lazard and Sylvain Lefebvre. ACM Transactions on Graphics 39(4) (Proc. SIGGRAPH 2020) [[link]](https://hal.inria.fr/hal-02568677)

## Authors of the code

- [Samuel Hornus](https://members.loria.fr/samuel.hornus/) (Inria, France)
- [Tim Kuipers](https://kuipers.weblog.tudelft.nl/) (Ultimaker, TU Delft, Netherlands)

## Visualization and analysis

Code and scripts for visualizing and analyzing variable-width GCode can be
found [here](https://github.com/BagelOrb/ToolpathVisualizer).

## Compilation dependencies

	- CMAKE
	- TCLAP http://tclap.sourceforge.net
	- Clipper http://angusj.com/delphi/clipper.php
	- CGAL https://www.cgal.org
	- Cairo https://www.cairographics.org

BOOST Voronoi support should be available in the future as an alternative to CGAL.

## Compiling

It should work alright on macOS and Ubuntu.

	mkdir build
	cd build
	cmake ..
	make

You may have to adjust the CMakeLists.txt.
Please send us your ``pull`` requests!

There is a CMake option to disable the use of Cairo. Then no PDF is ever output.

## Usage

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

## Dataset

The ``input/`` directory has some sample input files.

The set of 300 input files in input/dataset/ is the one used in the paper and
originates from the following work:

> _A framework for adaptive width control of dense contour-parallel toolpaths in fused deposition modeling_. Tim Kuipers, Eugeni L. Doubrovski, Jun Wu, and Charlie C. L. Wang. 2020. [arXiv:2004.13497 [cs.GR]](https://arxiv.org/abs/2004.13497) In submission.
