# distributed-anti-aliasing
Simple supersampling AA implementation designed to run on a linux cluster with MPICH and Open MP installed.
# How to run it from a linux terminal
Use the following command:

mpiexec –machinefile `machines` –np `num_of_cores` `executable` `input’s file name` `output’s file name` `input height` `input width` `superblock’s edge’s size` `num of iterations` `num of colors`

Where
- `machines` is the relative path to the file which contains the names of the nodes of the linux cluster
- `num_of_cores` is the number of the CPUs (cores) to be utilized on each node
- `executable` is the name of the executable produced by running the makefile
- `input’s file name` is a RAW image file either a grayscale or an RGB one
- `output’s file name` is as the name suggests the name of the output file
- `input height` is the input image's height in pixels
- `input width` the width in pixels, must be smaller than the height
- `superblock’s edge’s size` is the size of the edge of a superblock in pixels (if this argument is not given then the algorithm is ran on a single machine)
- `num of iterations` is the number of times the filter is going to be applied on the image
- `num of colors` 1 for grayscale and 3 for RGB input

# Further info
The input image is considered to be a grid of pixels which is subsequently divided into square superblocks. The superblock is itself divided into blocks that amount to the number of nodes in the cluster and to which are then assigned. The iterative process may stop in less time than defined by the `num of iterations` argument in the case where 2 succesive iterations have results that differ by less than a hardcoded lower bound (static variable `lower_bound_convergence`).
