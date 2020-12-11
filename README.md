C++ demo code and application for "[You Can Find Geodesic Paths in Triangle Meshes by Just Flipping Edges](https://nmwsharp.com/research/flip-geodesics/)", by [Nicholas Sharp](https://nmwsharp.com/) and [Keenan Crane](http://keenan.is/here) at SIGGRAPH Asia 2020.

- PDF: [link](https://nmwsharp.com/media/papers/flip-geodesics/flip_geodesics.pdf)
- Project: [link](https://nmwsharp.com/research/flip-geodesics/)
- Talk: coming soon!

![shorten path to geodesic](https://raw.githubusercontent.com/nmwsharp/flip-geodesics-demo/master/media/make_geodesic.jpg)

This algorithm takes as input a path (or loop/network of paths) along the edges of a triangle mesh, and as output straightens that path to be a _geodesic_ (i.e. a straight line, or equivalently a locally-shortest path along a surface). The procedure runs in milliseconds, is quite robust, comes with a strong guarantee that no new crossings will be created in the path, and as an added benefit also generates a triangulation on the surface which conforms to the geodesic. Additionally, it even enables the construction of Bézier curves on a surface! 

The main algorithm is implemented in [geometry-central](http://geometry-central.net/surface/algorithms/flip_geodesics/). This repository contains a simple demo application including a GUI to invoke that implementation.

If this code contributes to academic work, please cite:

```
@article{sharp2020you,
  title={You can find geodesic paths in triangle meshes by just flipping edges},
  author={Sharp, Nicholas and Crane, Keenan},
  journal={ACM Transactions on Graphics (TOG)},
  volume={39},
  number={6},
  pages={1--15},
  year={2020},
  publisher={ACM New York, NY, USA}
}
```

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

## Usage

After building, run the demo application like `./bin/flip_geodesics /path/to/your/mesh.obj`. The accepted mesh types are [documented](http://geometry-central.net/surface/utilities/io/) in geometry-central (for now: obj, ply, off, stl). The input should be a manifold triangle mesh.

![gui screencap](https://raw.githubusercontent.com/nmwsharp/flip-geodesics-demo/master/media/gui_screencap.png)

### Basic input

The simplest way to construct a path is to select two endpoints; the app will run Dijkstra's algorithm to generate an initial end path between the points. Click  <kbd>construct new Dijkstra path from endpoints</kbd> -- the app will then guide you to ctrl-click on two vertices (or instead enter vertex indices).

### Advanced input

The app also offers several methods to construct more interesting initial paths.

<details>
  <summary>Click to expand!</summary>

#### Fancy paths

This method allows you to manually construct more interesting paths along the surface beyond just Dijkstra paths between endpoints. Open the menu via the  <kbd>construct fancy path</kbd> dropdown.

  You can input a path by selecting a sequential list of points on the surface. Once some sequence of points has been added, selecting  <kbd>new path from these points</kbd> will run Dijkstra's algorithm between each consecutive pair of points in the list to create the initial path. The  <kbd>push vertex</kbd> button adds a point to the sequence, while  <kbd>pop vertex</kbd> removes the most recent point.

  Checking  <kbd>created closed path</kbd> will connect the first and last points of the path to form a closed loop. Checking  <kbd>mark interior vertices</kbd> will pin the curve to the selected vertex list during shortening.

#### Speciality loaders

Additionally, several loaders are included for other possible file formats. These interfaces are a bit ad-hoc, but are included to hopefully facilitate your own experiments and testing!

-  <kbd>load edge set</kbd> Create a path by specifying a list of collection of edges which make up the path. Loads from a file in the current directory called `path_edges.txt`, where each line contains two, space-separated 0-indexed vertex indices which are the endpoints of some edge in the path.  Additionally, if `marked_vertices.txt` is present it should hold one vertex index per line, which will be pinned during straightening.
-  <kbd>load line list obj</kbd>  Create a path network from [line elements](https://en.wikipedia.org/wiki/Wavefront_.obj_file#Line_elements) in an .obj file. Loads from the same file as the initial input to the program, which must be an .obj file. The line indices in this file must correspond to mesh vertex indices.
-  <kbd>load Dijkstra list</kbd> Create a path network from one or more Dijkstra paths between vertices. Loads from a file in the current directory called `path_pairs.txt`, where each line contains two, space-separated 0-indexed vertex indices which are the endpoints of the path. If this file has many lines, a network will be created. 
-  <kbd>load UV cut</kbd>  Create a path network from cuts (aka discontinuities aka island boundaries) in a UV map. Loads from the same file as the initial input to the program, which must be an .obj file with UVs specified.
-  <kbd>load seg cut</kbd> Create a path network from the boundary of a per-face segmentation. Loads from a plaintext file in the current directory called `cut.seg`, where each line corresponds gives an integer segmentation ID for a face.

</details>

### FlipOut straightening

Once a path/loop/network has been loaded, the  <kbd>make geodesic</kbd> button will straighten it to a geodesics. The optional checkboxes limit the number of `FlipOut()` iterations, or the limit the total length decrease. See the Visualization section 

To verify the resulting path is really an exact polyhedral geodesic, the <kbd>check path</kbd> button will measure the swept angles on either side of the path, and print the smallest such angle to the terminal. Mathematically, the FlipOut procedure is guaranteed to yield a geodesic; (very rare) failures in practice are due to the inaccuracies of floating point computation on degenerate meshes.

Expanding the  <kbd>extras</kbd> dropdown gives additional options:

- **Bézier subdivision** iteratively constructs a smooth Bézier curve, treating the input path as control points. This option should be used when a single path between two endpoints is registered.
- **Mesh improvement** performs intrinsic refinement to improve the quality of the resulting triangulation.

### Visualization

The app uses [polyscope](http://polyscope.run/) for visualization; see the documentation there for general details about the interface.

Once as path is loaded, it will be drawn with a red curve along the surface. Expanding the <kbd>path edges</kbd> dropdown on the leftmost menu allows modifying the color and curve size, etc.

By default, only the path itself is drawn, the  <kbd>show intrinsic edges</kbd> checkbox draws _all_ edges in the underlying intrinsic triangulation, in yellow (which can again be tweaked via the options on the left).

The  <kbd>export path lines</kbd> button writes a file called `lines_out.obj`, containing line entries for the path network. Note that you probably want to export _after_ straightening, to export the geodesic path network.

## Command line interface

**coming soon :)**

The executable also supports scripted usage via a simple command line interface. See the `flip_geodesics --help` for additional documentation. This functionality essentially mimics the GUI usage described above; see there for details.
