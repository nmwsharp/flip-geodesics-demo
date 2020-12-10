C++ demo code and application for "[You Can Find Geodesic Paths in Triangle Meshes by Just Flipping Edges](https://nmwsharp.com/research/flip-geodesics/)", by [Nicholas Sharp](https://nmwsharp.com/) and [Keenan Crane](http://keenan.is/here) at SIGGRAPH Asia 2020.

- PDF: [link](https://nmwsharp.com/media/papers/flip-geodesics/flip_geodesics.pdf)
- Project: [link](https://nmwsharp.com/research/flip-geodesics/)
- Talk: TODO

The main algorithm is implemented in [geometry-central](http://geometry-central.net/). This repository contains a simple demo application including a GUI to invoke that implementation.

## Cloning and building

On unix-like environments, run:
```sh
git clone --recursive https://github.com/nmwsharp/flip-geodesics-demo.git
cd flip-geodesics-demo
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
./bin/flip_geodesics /path/to/your/mesh.obj
```

The provided `CMakeLists.txt` should also generate solutions which compile in Visual Studio (see many tutorials online).

## Basic usage

### Basic input

The simplest way to construct a path is to select two endpoints; the app will run Dijkstra's algorithm to generate an initial end path between the points. Click  <kbd>construct new Dijkstra path from endpoints</kbd> -- the app will then guide you to ctrl-click on two vertices (or instead enter vertex indices).

### Advanced input

The app also offers several methods to construct more interesting initial paths.

<details>
  <summary>Click to expand!</summary>

#### Fancy paths

This method allows you to manually construct more interesting paths along the surface beyond just Dijkstra paths between endpoints. Open the menu via the [construct fancy path] dropdown.

  You can input a path by selecting a sequential list of points on the surface. Once some sequence of points has been added, selecting [new path from these points] will run Dijkstra's algorithm between each consecutive pair of points in the list to create the initial path. The [push vertex] button adds a point to the sequence, while [pop vertex] removes the most recent point.

  Checking [created closed path] will connect the first and last points of the path to form a closed loop. Checking [mark interior vertices] will pin the curve to the selected vertex list during shortening.

#### Speciality loaders

Additionally, several loaders are included for other possible file formats. These interfaces are a bit ad-hoc, but are included to hopefully facilitate your own experiments and testing!

- [load edge set] Create a path by specifying a list of collection of edges which make up the path. Loads from a file in the current directory called `path_edges.txt`, where each line contains two, space-separated 0-indexed vertex indices which are the endpoints of some edge in the path.  Additionally, if `marked_vertices.txt` is present it should hold one vertex index per line, which will be pinned during straightening.
- [load line list obj]  Create a path network from [line elements](https://en.wikipedia.org/wiki/Wavefront_.obj_file#Line_elements) in an .obj file. Loads from the same file as the initial input to the program, which must be an .obj file. The line indices in this file must correspond to mesh vertex indices.
- [load Dijkstra list] Create a path network from one or more Dijkstra paths between vertices. Loads from a file in the current directory called `path_pairs.txt`, where each line contains two, space-separated 0-indexed vertex indices which are the endpoints of the path. If this file has many lines, a network will be created. 
- [load UV cut]  Create a path network from cuts (aka discontinuities aka island boundaries) in a UV map. Loads from the same file as the initial input to the program, which must be an .obj file with UVs specified.
- [load seg cut] Create a path network from the boundary of a per-face segmentation. Loads from a plaintext file in the current directory called `cut.seg`, where each line corresponds gives an integer segmentation ID for a face.

</details>

### FlipOut straightening

Once a path/loop/network has been loaded, the [make geodesics] button will straighten it to a geodesics. The optional checkboxes limit the number of `FlipOut()` iterations, or the limit the total length decrease. See the Visualization section 

To verify the resulting path is really an exact polyhedral geodesic, the [check path] button will measure the swept angles on either side of the path, and print the smallest such angle to the terminal. Mathematically, the FlipOut procedure is guaranteed to yield a geodesic; (very rare) failures in practice are due to the inaccuracies of floating point computation on degenerate meshes.

Expanding the [extras] dropdown gives additional options:

- **Bezier subdivision** iteratively constructs a smooth Bezier curve, treating the input path as control points. This option should be used when a single path between two endpoints is registered.
- **Mesh improvement** performs intrinsic refinment to improve the quality of the resulting triangulation

### Visualization

The app uses [polyscope](http://polyscope.run/) for visualization; see the documentation there for general details about the interface.

Once as path is loaded, it will be drawn with a red curve along the surface. Expanding the [path edges] dropdown on the leftmost menu allows modifying the color and curve size, etc.

By default, only the path itself is drawn, the [show intrinsic edges] checkbox draws _all_ edges in the underlying intrinsic triangulation, in yellow (which can again be tweaked via the options on the left).

The [export path lines] button writes a file called `lines_out.obj`, containing line entries for the path network. Note that you problably want to export _after_ straightening, to export the geodesic path network.


## Command line interface

The executable also supports scripted usage via a simple command line interface. See the `flip_geodesics --help` for additional documentation. This functionality essentially mimicks the GUI usage described above; see there for details.
