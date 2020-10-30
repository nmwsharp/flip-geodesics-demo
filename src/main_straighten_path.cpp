#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/polygon_soup_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/timing.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

// Are we currently processing a path/loop/source?
enum class FlipMode { None = 0, Path, Loop, Source };
FlipMode currMode = FlipMode::None;

// Data that might be defined for modes above
std::unique_ptr<EdgeNetwork> edgeNetwork;

// UI parameters
bool iterativeShortenUseIterationCap = false;
int iterativeShortenIterationCap = 1;
bool straightenAtMarked = true;
bool useIterativeShortenLengthLim = false;
float iterativeShortenLengthLim = 0.5;

int nBezierIters = 3;

bool vizAllIntrinsicEdges = false;
bool maintainDelaunay = false;
bool makeMovie = false;
bool flipMovie = false;
bool ribbonPerEdge = true;
bool ribbonAllEdge = false;
float angleEPS = 1e-5;
float splitAngleDeg = 10;
float refineAreaThresh = std::numeric_limits<float>::infinity();
float refineAngleThresh = 25.;
int maxInsertions = -1;

float tubeWidth = -1.;

// ====== Path related stuff

void updatePathViz() {
  if (!edgeNetwork) {
    polyscope::error("tried to visualize path, but don't have");
    return;
  }

  // Copy settings
  edgeNetwork->vizAllIntrinsicEdges = vizAllIntrinsicEdges;

  edgeNetwork->visualizePath();
}

void createPathFromPoints() {

  long long int iVStart = psMesh->selectVertex();
  long long int iVEnd = psMesh->selectVertex();

  if (iVStart == -1 || iVEnd == -1)
    return;

  edgeNetwork = EdgeNetwork::constructFromDijkstraPath(
      *mesh, *geometry, mesh->vertex(iVStart), mesh->vertex(iVEnd));
  if (edgeNetwork == nullptr) {
    polyscope::warning("could not initialize edge path between vertices");
    return;
  }
  edgeNetwork->posGeom = geometry.get();
  currMode = FlipMode::Path;

  updatePathViz();
}

void createPathFromLast() {

  edgeNetwork =
      EdgeNetwork::constructFromSavedPath(*mesh, *geometry, "last_path.ply");
  if (edgeNetwork == nullptr) {
    polyscope::warning("could not initialize path from file");
    return;
  }
  edgeNetwork->posGeom = geometry.get();

  currMode = FlipMode::Path;
  updatePathViz();
}

void createRandomLoop() {

  int nTests = 1;
  int totalTime = 0;
  int maxTime = 0;
  size_t iterLim = iterativeShortenUseIterationCap
                       ? iterativeShortenIterationCap
                       : INVALID_IND;
  double lengthLim =
      useIterativeShortenLengthLim ? iterativeShortenLengthLim : 0.;

  for (int i = 0; i < nTests; i++) {
    std::cout << "=========== Test " << i + 1 << " of " << nTests << std::endl;
    EdgeData<char> edgeSet(*mesh, false);
    VertexData<char> extraMarkedVertices(*mesh, false);

    // pick a random plane through the origin
    Vector3 N{randomReal(-1., 1.), randomReal(-1., 1.), randomReal(-1., 1.)};
    N.normalize();

    // for each vertex, count the number of incident edges crossing the plane
    VertexData<int> endpointDegree(*mesh, 0);
    for (Vertex v : mesh->vertices()) {
      for (Edge e : v.adjacentEdges()) {
        Vertex v0 = e.halfedge().vertex();
        Vertex v1 = e.halfedge().twin().vertex();
        Vector3 p0 = geometry->inputVertexPositions[v0];
        Vector3 p1 = geometry->inputVertexPositions[v1];

        if (dot(N, p0) * dot(N, p1) < 0.) {
          endpointDegree[v]++;
        }
      }
    }

    // add edges with endpoints crossing the plane, but that
    // are not incident on a degree-1 endpoint, and do not
    // create a vertex of degree 3 or greater
    VertexData<int> endpointCount(*mesh, 0);
    for (Edge e : mesh->edges()) {
      Vertex v0 = e.halfedge().vertex();
      Vertex v1 = e.halfedge().twin().vertex();
      Vector3 p0 = geometry->inputVertexPositions[v0];
      Vector3 p1 = geometry->inputVertexPositions[v1];

      if (dot(N, p0) * dot(N, p1) < 0. && endpointDegree[v0] != 1 &&
          endpointDegree[v1] != 1 && endpointCount[v0] < 2 &&
          endpointCount[v1] < 2) {
        edgeSet[e] = true;
        endpointCount[v0]++;
        endpointCount[v1]++;
      }
    }

    edgeNetwork = EdgeNetwork::constructFromEdgeSet(*mesh, *geometry, edgeSet,
                                                    extraMarkedVertices);
    if (edgeNetwork == nullptr) {
      polyscope::warning("could not initialize path from random plane");
      return;
    }
    edgeNetwork->posGeom = geometry.get();

    // auto timer = NOW;
    // edgeNetwork->iterativeShorten(iterLim, lengthLim);
    // edgeNetwork->getPathPolyline();
    // auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(NOW
    // - timer); int time = elapsed.count(); totalTime += time; maxTime =
    // std::max( time, maxTime );
  }
  std::cout << "max time: " << maxTime / 1000. << "ms" << std::endl;
  std::cout << "avg time: " << totalTime / (double)nTests / 1000. << "ms"
            << std::endl;

  currMode = FlipMode::Loop;
  updatePathViz();
}

void createPathFromEdgeSet() {

  EdgeData<char> edgeSet(*mesh, false);
  VertexData<char> extraMarkedVertices(*mesh, false);

  { // Load the edge set from file
    std::ifstream inStream("path_edges.txt");
    if (!inStream) {
      polyscope::error("could not read path_edges.txt");
      return;
    }

    for (std::string line; std::getline(inStream, line);) {
      std::istringstream lineStream(line);

      size_t indA, indB;
      lineStream >> indA;
      lineStream >> indB;

      for (Halfedge he : mesh->vertex(indA).incomingHalfedges()) {
        if (he.vertex().getIndex() == indB) {
          edgeSet[he.edge()] = true;
        }
      }
    }
  }

  { // Load extra marked vertices from file
    std::ifstream inStream("marked_vertices.txt");
    if (!inStream) {
      polyscope::warning("could not read marked_vertices.txt");
    }

    for (std::string line; std::getline(inStream, line);) {
      std::istringstream lineStream(line);

      size_t ind;
      lineStream >> ind;
      extraMarkedVertices[mesh->vertex(ind)] = true;
    }
  }

  edgeNetwork = EdgeNetwork::constructFromEdgeSet(*mesh, *geometry, edgeSet,
                                                  extraMarkedVertices);
  if (edgeNetwork == nullptr) {
    polyscope::warning("could not initialize path from file");
    return;
  }
  edgeNetwork->posGeom = geometry.get();

  currMode = FlipMode::Path;
  updatePathViz();
}

void createPathFromObjLines() {

  auto findHalfedge = [&](Vertex vA, Vertex vB) {
    for (Halfedge he : vA.outgoingHalfedges()) {
      if (he.twin().vertex() == vB)
        return he;
    }
    return Halfedge();
  };

  std::vector<std::vector<Halfedge>> paths;

  { // Load the edge set from file
    std::ifstream inStream("path_list.obj");
    if (!inStream) {
      polyscope::error("could not read path_list.obj");
      return;
    }

    for (std::string line; std::getline(inStream, line);) {
      if (line.size() < 2 || line[0] != 'l' || line[1] != ' ')
        continue; // look for line ('l' first char)

      std::istringstream lineStream(
          line.substr(2)); // all but first char and space

      // parse out list of indices
      std::vector<size_t> inds;
      std::string token;
      while (std::getline(lineStream, token, ' ')) {
        std::stringstream tokenStream(token);
        size_t i;
        tokenStream >> i;
        inds.push_back(i - 1);
      }

      // build vertices
      paths.emplace_back();
      std::vector<Halfedge> &path = paths.back();
      for (size_t i = 1; i < inds.size(); i++) {
        Vertex vA = mesh->vertex(inds[i - 1]);
        Vertex vB = mesh->vertex(inds[i]);

        Halfedge he = findHalfedge(vA, vB);
        if (he == Halfedge()) {
          polyscope::warning("vertices " + std::to_string(inds[i - 1]) +
                             " and " + std::to_string(inds[i]) +
                             " are not connected");
          return;
        }
        path.push_back(he);
      }

      // Try to close a loop
      Halfedge lastHe =
          findHalfedge(mesh->vertex(inds.back()), mesh->vertex(inds.front()));
      if (lastHe != Halfedge()) {
        std::cout << "    closing loop with halfedge " << lastHe << std::endl;
        path.push_back(lastHe);
      }

      std::cout << "  ...found path with " << path.size() << " segments."
                << std::endl;
    }
  }

  std::cout << "Loaded line list with " << paths.size() << " paths."
            << std::endl;

  edgeNetwork.reset(new EdgeNetwork(*mesh, *geometry, paths));
  if (edgeNetwork == nullptr) {
    polyscope::warning("could not initialize path from file");
    return;
  }
  edgeNetwork->posGeom = geometry.get();

  currMode = FlipMode::Path;
  updatePathViz();
}

void createPathFromDijkstraList() {

  // Create an  (initially-empty) edge network
  edgeNetwork =
      std::unique_ptr<EdgeNetwork>(new EdgeNetwork(*mesh, *geometry, {}));
  edgeNetwork->posGeom = geometry.get();

  std::cout << "!!! NOTE: Constructing Dijkstra paths on Delaunay-fied mesh"
            << std::endl;
  edgeNetwork->makeDelaunay();

  { // Load the edge set from file
    std::ifstream inStream("path_pairs.txt");
    if (!inStream) {
      polyscope::error("could not read path_pairs.txt");
      return;
    }

    for (std::string line; std::getline(inStream, line);) {
      std::istringstream lineStream(line);

      size_t indA, indB;
      lineStream >> indA;
      lineStream >> indB;

      Vertex vA = edgeNetwork->tri->intrinsicMesh->vertex(indA);
      Vertex vB = edgeNetwork->tri->intrinsicMesh->vertex(indB);

      std::cout << "loading path from " << vA << " to " << vB << std::endl;

      std::vector<Halfedge> path = shortestEdgePath(*edgeNetwork->tri, vA, vB);
      edgeNetwork->addPath(path);
    }
  }

  currMode = FlipMode::Path;
  updatePathViz();
}

void createPathFromSeg() {

  // read face segmentation inds
  std::ifstream inStream("cut.seg");
  if (!inStream) {
    polyscope::warning("could not read cut.seg");
  }
  FaceData<int> segs(*mesh);
  size_t iF = 0;
  for (std::string line; std::getline(inStream, line);) {
    std::istringstream lineStream(line);
    int ind;
    lineStream >> ind;
    if (iF >= mesh->nFaces()) {
      polyscope::warning("segmentation file doesn't match number of faces");
      return;
    }
    segs[iF] = ind;
    iF++;
  }

  psMesh->addFaceScalarQuantity("segmentation", segs);

  // Make cut along boundary
  EdgeData<char> edgeSet(*mesh, false);
  for (Halfedge he : mesh->halfedges()) {
    if (he.edge().isBoundary())
      continue;
    if (segs[he.face()] != segs[he.twin().face()]) {
      edgeSet[he.edge()] = true;
    }
  }

  VertexData<char> extraMarkedVertices(*mesh, false); // none

  edgeNetwork = EdgeNetwork::constructFromEdgeSet(*mesh, *geometry, edgeSet,
                                                  extraMarkedVertices);
  if (edgeNetwork == nullptr) {
    polyscope::warning("could not initialize path from file");
    return;
  }
  edgeNetwork->posGeom = geometry.get();

  currMode = FlipMode::Path;
  updatePathViz();
}

void createPathFromUVCut() {

  // (re)-load the obj file, and pull out corner coords
  PolygonSoupMesh reloadMesh("uv_mesh.obj");
  if (reloadMesh.paramCoordinates.empty()) {
    polyscope::warning("could not load UVs from mesh file");
    return;
  }

  CornerData<Vector2> uvCoords(*mesh);
  for (Face f : mesh->faces()) {
    size_t iFace = f.getIndex();
    size_t iC = 0;
    for (Corner c : f.adjacentCorners()) {
      Vector2 uv = reloadMesh.paramCoordinates[iFace][iC];
      uvCoords[c] = uv;
      iC++;
    }
  }
  psMesh->addParameterizationQuantity("loaded UVs", uvCoords);

  // Detect cut as gap in UVs
  EdgeData<char> edgeSet(*mesh, false);
  for (Halfedge he : mesh->halfedges()) {
    if (he.edge().isBoundary())
      continue;
    if (uvCoords[he.corner()] != uvCoords[he.twin().next().corner()]) {
      edgeSet[he.edge()] = true;
    }
  }

  VertexData<char> extraMarkedVertices(*mesh, false); // none

  edgeNetwork = EdgeNetwork::constructFromEdgeSet(*mesh, *geometry, edgeSet,
                                                  extraMarkedVertices);
  if (edgeNetwork == nullptr) {
    polyscope::warning("could not initialize path from file");
    return;
  }
  edgeNetwork->posGeom = geometry.get();

  currMode = FlipMode::Path;
  updatePathViz();
}

void checkPath() {
  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }

  edgeNetwork->validate();

  double minAngle = edgeNetwork->minAngleIsotopy();
  std::cout << "min angle = " << minAngle << std::endl;
  if (minAngle < (M_PI - edgeNetwork->EPS_ANGLE)) {
    polyscope::warning("min angle in path is " + std::to_string(minAngle));
  }
}

void makeDelaunay() {
  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }
  edgeNetwork->makeDelaunay();
}

void locallyShorten() {
  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }

  // reset counters
  edgeNetwork->nFlips = 0;
  edgeNetwork->nShortenIters = 0;

  edgeNetwork->delaunayAwayFromPath = maintainDelaunay;
  edgeNetwork->EPS_ANGLE = angleEPS;
  edgeNetwork->straightenAroundMarkedVertices = straightenAtMarked;

  size_t iterLim = iterativeShortenUseIterationCap
                       ? iterativeShortenIterationCap
                       : INVALID_IND;
  double lengthLim =
      useIterativeShortenLengthLim ? iterativeShortenLengthLim : 0.;

  if (iterativeShortenUseIterationCap) {
    edgeNetwork->iterativeShorten(iterLim, lengthLim);
  } else {

    START_TIMING(shorten)
    edgeNetwork->iterativeShorten(iterLim, lengthLim);
    edgeNetwork->getPathPolyline();
    FINISH_TIMING_PRINT(shorten)

    checkPath();
  }

  std::cout << "shortening performed " << edgeNetwork->nShortenIters
            << " iterations, with a total of " << edgeNetwork->nFlips
            << " flips. " << std::endl;
}

void writeBoundaryAsSVG(ManifoldSurfaceMesh &mesh,
                        EmbeddedGeometryInterface &geom, std::string filename) {
  using namespace svg;

  geom.requireVertexPositions();

  // Setup
  double d = 1000;
  Dimensions dimensions(d, d);
  Document doc(filename, Layout(dimensions, Layout::BottomLeft));

  // Creates a single shape in output
  Polyline shape(Fill(Color::Silver), Stroke(.5, Color::Red));

  for (Vertex v : mesh.boundaryLoop(0).adjacentVertices()) {
    Vector3 p = geom.vertexPositions[v];
    p *= d;
    shape << Point(p.x, p.z);
  }
  doc << shape;

  doc.save();
}

void flattenCutSurface() {

  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }

  EdgeData<char> isOnCut(edgeNetwork->tri->mesh);
  for (Edge e : edgeNetwork->tri->mesh.edges()) {
    isOnCut[e] = edgeNetwork->edgeInPath(e);
  }

  // Cut the surface
  std::unique_ptr<ManifoldSurfaceMesh> cutMesh;
  HalfedgeData<Halfedge> parentHe;
  std::tie(cutMesh, parentHe) =
      cutAlongEdges(*edgeNetwork->tri->intrinsicMesh, isOnCut);

  /*
  // == Extrinsic
  // Copy geometry from the old mesh
  VertexData<Vector3> newPos(*cutMesh);
  for (Halfedge he : cutMesh->halfedges()) {
    Halfedge oldHe = parentHe[he];
    if (oldHe != Halfedge()) {
      newPos[he.vertex()] =
  geometry->inputVertexPositions[oldHe.vertex().getIndex()];
    }
  }
  VertexPositionGeometry cutGeom(*cutMesh, newPos);

  auto psCutMesh =
      polyscope::registerSurfaceMesh("cut mesh", cutGeom.inputVertexPositions,
  cutMesh->getFaceVertexList());

  // VertexData<Vector2> flatCoords = parameterizeDisk(*edgeNetwork->tri);
  VertexData<Vector2> flatCoords = parameterizeDisk(cutGeom);
  psCutMesh->addVertexParameterizationQuantity("flat coords", flatCoords);
  */

  // == Intrinsic
  EdgeData<double> newLen(*cutMesh);
  for (Halfedge he : cutMesh->halfedges()) {
    Halfedge oldHe = parentHe[he];
    if (oldHe != Halfedge()) {
      newLen[he.edge()] =
          edgeNetwork->tri->intrinsicEdgeLengths[oldHe.edge().getIndex()];
    }
  }
  EdgeLengthGeometry cutGeom(*cutMesh, newLen);

  // Refine
  SignpostIntrinsicTriangulation tri(*cutMesh, cutGeom);
  tri.delaunayRefine();

  // TODO compress here

  // VertexData<Vector2> flatCoords = parameterizeDisk(*edgeNetwork->tri);
  VertexData<Vector2> coords = parameterizeDisk(*(tri.intrinsicMesh), tri);

  // Copy coords to X-Y
  VertexData<Vector3> flatCoords(tri.mesh); // TODO use this once compress works
  std::vector<Vector3> flatCoordsVec;
  for (Vertex v : tri.mesh.vertices()) {
    Vector3 c{coords[v].x, 0., coords[v].y};
    flatCoordsVec.push_back(c);
    flatCoords[v] = c;
  }

  auto psCutMesh = polyscope::registerSurfaceMesh("cut mesh", flatCoordsVec,
                                                  tri.mesh.getFaceVertexList());

  // Write to file
  VertexPositionGeometry flatCoordsGeom(tri.mesh, flatCoords);
  WavefrontOBJ::write("flat_mesh.obj", flatCoordsGeom);
  writeBoundaryAsSVG(*tri.intrinsicMesh, flatCoordsGeom, "flat_mesh.svg");
}

void distanceToCurve() {

  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }
  ManifoldSurfaceMesh &intMesh = *(edgeNetwork->tri->intrinsicMesh);

  // == Refine to make sure we have a decent quality triangulation
  // edgeNetwork->delaunayRefine(areaThresh, maxInsertions);

  // == Compute distance from the curve
  HeatMethodDistanceSolver solver(*edgeNetwork->tri);

  // Build the boundary conditions
  VertexData<double> rhs(intMesh, 0.);
  std::vector<Vertex> sourceVerts;
  for (Vertex v : edgeNetwork->tri->mesh.vertices()) {
    for (Edge e : v.adjacentEdges()) {
      if (edgeNetwork->edgeInPath(e)) {
        sourceVerts.push_back(v);
        break;
      }
    }
  }
  for (Edge e : intMesh.edges()) {
    if (edgeNetwork->edgeInPath(e)) {
      double l = edgeNetwork->tri->intrinsicEdgeLengths[e];
      rhs[e.halfedge().vertex()] += l / 2;
      rhs[e.halfedge().twin().vertex()] += l / 2;
    }
  }

  Vector<double> distanceVec = solver.computeDistanceRHS(rhs.toVector());
  distanceVec =
      distanceVec.array() - (double)distanceVec.minCoeff(); // lazy shift
  VertexData<double> distance(intMesh, distanceVec);

  // == Make the triangulation explicit

  // Interpolate distance on each edge split
  auto updateDistOnSplit = [&](Edge oldE, Halfedge newHe1, Halfedge newHe2) {
    double lenHe1 = edgeNetwork->tri->intrinsicEdgeLengths[newHe1.edge()];
    double lenHe2 = edgeNetwork->tri->intrinsicEdgeLengths[newHe2.edge()];
    double totalLen = lenHe1 + lenHe2;

    double newDist = (lenHe2 / totalLen) * distance[newHe1.twin().vertex()] +
                     (lenHe1 / totalLen) * distance[newHe2.twin().vertex()];

    distance[newHe1.vertex()] = newDist;
  };
  auto callbackRef = edgeNetwork->tri->edgeSplitCallbackList.insert(
      std::end(edgeNetwork->tri->edgeSplitCallbackList), updateDistOnSplit);

  // Split bent edges until the triangulation is explicit
  edgeNetwork->splitBentEdges(
      splitAngleDeg, maxInsertions == -1 ? INVALID_IND : maxInsertions);

  edgeNetwork->tri->edgeSplitCallbackList.erase(
      callbackRef); // remove the callback we registered

  // TODO compress here

  // Copy geometry to a new extrinsic mesh
  std::vector<Vector3> pos3DVec;
  VertexData<Vector3> pos3D(intMesh);
  std::vector<double> denseDist;
  std::vector<Vector3> sourceVertPos;
  for (Vertex v : intMesh.vertices()) {
    Vector3 p = edgeNetwork->tri->vertexLocations[v].interpolate(
        geometry->inputVertexPositions);
    pos3DVec.push_back(p);
    pos3D[v] = p;
    denseDist.push_back(distance[v]);
  }
  for (Vertex v : sourceVerts) {
    sourceVertPos.push_back(edgeNetwork->tri->vertexLocations[v].interpolate(
        geometry->inputVertexPositions));
  }

  auto psIntMesh = polyscope::registerSurfaceMesh("int mesh", pos3DVec,
                                                  intMesh.getFaceVertexList());
  psIntMesh->addVertexDistanceQuantity("distance to curve", denseDist);
  polyscope::registerPointCloud("source verts", sourceVertPos);

  // == Output

  // Copy coords to u
  CornerData<Vector2> coords(intMesh);
  for (Vertex v : intMesh.vertices()) {
    Vector2 x{distance[v], 0.};
    for (Corner c : v.adjacentCorners()) {
      coords[c] = x;
    }
  }

  // Write to file
  VertexPositionGeometry outGeom(intMesh, pos3D);
  WavefrontOBJ::write("mesh_with_dist.obj", outGeom, coords);
}

void poissonFromCurve() {

  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }
  ManifoldSurfaceMesh &intMesh = *(edgeNetwork->tri->intrinsicMesh);
  IntrinsicGeometryInterface &intGeom = *(edgeNetwork->tri);

  // Operators
  intGeom.requireVertexIndices();
  intGeom.requireCotanLaplacian();
  intGeom.requireVertexLumpedMassMatrix();
  SparseMatrix<double> L = intGeom.cotanLaplacian;
  SparseMatrix<double> M = intGeom.vertexLumpedMassMatrix;
  size_t N = L.rows();

  // Set up boundary conditions
  std::vector<Vertex> pinnedVerts;
  for (Vertex v : intMesh.vertices()) {
    for (Edge e : v.adjacentEdges()) {
      if (edgeNetwork->edgeInPath(e)) {
        pinnedVerts.push_back(v);
        break;
      }
    }
  }

  // Dirichlet BCs
  size_t NBoundary = pinnedVerts.size();
  Vector<double> rhsVals = M * Vector<double>::Ones(N);
  Vector<double> bcVals = Vector<double>::Zero(NBoundary);
  Vector<bool> setAMembership = Vector<bool>::Ones(N);
  for (Vertex v : pinnedVerts) {
    size_t i = intGeom.vertexIndices[v];
    setAMembership(i) = false;
  }

  // Solve the poisson problem
  BlockDecompositionResult<double> decomp =
      blockDecomposeSquare(L, setAMembership, true);
  Vector<double> rhsValsA, rhsValsB;
  decomposeVector(decomp, rhsVals, rhsValsA, rhsValsB);
  Vector<double> combinedRHS = rhsValsA - decomp.AB * bcVals;
  Vector<double> Aresult = solve(decomp.AA, combinedRHS);
  Vector<double> resultVec = reassembleVector(decomp, Aresult, bcVals);
  VertexData<double> result(intMesh, resultVec);

  /*
  // Poisson problem with source along curve
  Vector<double> rhsVals = Vector<double>::Zero(N);
  for (Vertex v : pinnedVerts) {
    size_t i = intGeom.vertexIndices[v];
    rhsVals(i) = 1.;
  }
  rhsVals = M * rhsVals;
  double shift = rhsVals.sum() / N;
  rhsVals = rhsVals.array() - shift;
  Vector<double> resultVec = solve(L, rhsVals);
  VertexData<double> result(intMesh, resultVec);
  */

  // == Make the triangulation explicit

  // Interpolate result on each edge split
  auto updateDistOnSplit = [&](Edge oldE, Halfedge newHe1, Halfedge newHe2) {
    double lenHe1 = edgeNetwork->tri->intrinsicEdgeLengths[newHe1.edge()];
    double lenHe2 = edgeNetwork->tri->intrinsicEdgeLengths[newHe2.edge()];
    double totalLen = lenHe1 + lenHe2;

    double newDist = (lenHe2 / totalLen) * result[newHe1.twin().vertex()] +
                     (lenHe1 / totalLen) * result[newHe2.twin().vertex()];

    result[newHe1.vertex()] = newDist;
  };
  auto callbackRef = edgeNetwork->tri->edgeSplitCallbackList.insert(
      std::end(edgeNetwork->tri->edgeSplitCallbackList), updateDistOnSplit);

  // Split bent edges until the triangulation is explicit
  edgeNetwork->splitBentEdges(
      splitAngleDeg, maxInsertions == -1 ? INVALID_IND : maxInsertions);

  edgeNetwork->tri->edgeSplitCallbackList.erase(
      callbackRef); // remove the callback we registered

  // TODO compress here

  // Copy geometry to a new extrinsic mesh
  std::vector<Vector3> pos3DVec;
  VertexData<Vector3> pos3D(intMesh);
  std::vector<double> denseResult;
  std::vector<Vector3> sourceVertPos;
  for (Vertex v : intMesh.vertices()) {
    Vector3 p = edgeNetwork->tri->vertexLocations[v].interpolate(
        geometry->inputVertexPositions);
    pos3DVec.push_back(p);
    pos3D[v] = p;
    denseResult.push_back(result[v]);
  }
  for (Vertex v : pinnedVerts) {
    sourceVertPos.push_back(edgeNetwork->tri->vertexLocations[v].interpolate(
        geometry->inputVertexPositions));
  }

  auto psIntMesh = polyscope::registerSurfaceMesh("int mesh", pos3DVec,
                                                  intMesh.getFaceVertexList());
  psIntMesh->addVertexScalarQuantity("result", denseResult);
  polyscope::registerPointCloud("source verts", sourceVertPos);

  // == Output

  // Copy coords to u
  CornerData<Vector2> coords(intMesh);
  for (Vertex v : intMesh.vertices()) {
    Vector2 x{result[v], 0.};
    for (Corner c : v.adjacentCorners()) {
      coords[c] = x;
    }
  }

  // Write to file
  VertexPositionGeometry outGeom(intMesh, pos3D);
  WavefrontOBJ::write("mesh_with_poisson.obj", outGeom, coords);
}

void crossFieldFromCurve() {

  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }
  ManifoldSurfaceMesh &intMesh = *(edgeNetwork->tri->intrinsicMesh);

  // == Refine to make sure we have a decent quality triangulation
  // edgeNetwork->delaunayRefine(areaThresh, maxInsertions);

  // == Compute cross field from the curve

  edgeNetwork->tri->requireHalfedgeVectorsInVertex();

  // Build the boundary conditions
  std::cout << "building BCs" << std::endl;
  std::vector<std::pair<Vertex, Vector2>> BCs;
  for (Vertex v : edgeNetwork->tri->mesh.vertices()) {
    for (Halfedge he : v.outgoingHalfedges()) {
      if (edgeNetwork->edgeInPath(he.edge())) {

        Vector2 bcVal =
            unit(edgeNetwork->tri->halfedgeVectorsInVertex[he]).pow(4.);
        BCs.emplace_back(v, bcVal);
        std::cout << "  setting bc " << v << " = " << bcVal << std::endl;
        break;
      }
    }
  }

  // Solve
  std::cout << "solving" << std::endl;
  throw std::runtime_error("temporarily disabled");
  VertexData<Vector2> crossField;
  // VertexData<Vector2> crossField =
  // computeSmoothestVertexDirectionField(*edgeNetwork->tri, BCs, 4);
  std::cout << "done solving" << std::endl;

  // Pull back to original mesh
  VertexData<Vector2> crossFieldOnOrig(*mesh, Vector2::zero());
  for (Vertex v : edgeNetwork->tri->mesh.vertices()) {
    if (edgeNetwork->tri->vertexLocations[v].type == SurfacePointType::Vertex) {
      Vector2 val = crossField[v];
      Vertex origV = edgeNetwork->tri->vertexLocations[v].vertex;
      crossFieldOnOrig[origV] = val;
    }
  }
  std::cout << "pulling back" << std::endl;

  // Set vertex tangent spaces
  std::cout << "prepping polyscope" << std::endl;
  geometry->requireVertexTangentBasis();
  VertexData<Vector3> vBasisX(*mesh);
  for (Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
  }
  psMesh->setVertexTangentBasisX(vBasisX);

  // Add vectors
  std::cout << "register" << std::endl;
  psMesh->addVertexIntrinsicVectorQuantity("really great vectors",
                                           crossFieldOnOrig, 4);
  std::cout << "done register" << std::endl;

  // polyscope::registerPointCloud("source verts", sourceVertPos);
}

void bezierSubdivide() {
  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }

  START_TIMING(bezier)
  edgeNetwork->bezierSubdivide(nBezierIters);
  edgeNetwork->getPathPolyline();
  FINISH_TIMING_PRINT(bezier)

  updatePathViz();
}

void delaunayRefine() {

  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }

  edgeNetwork->delaunayRefine(refineAreaThresh,
                              maxInsertions == -1 ? INVALID_IND : maxInsertions,
                              refineAngleThresh);
  std::cout << "refined mesh has "
            << edgeNetwork->tri->intrinsicMesh->nVertices() << " verts\n";

  updatePathViz();
}

void splitBentEdges() {

  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }

  edgeNetwork->splitBentEdges(
      splitAngleDeg, maxInsertions == -1 ? INVALID_IND : maxInsertions);

  updatePathViz();
}

void bakePathToMesh() {
  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }

  std::cout << "BAKING PATH TO MESH" << std::endl;

  // Split edges in the original extrinsic mesh to include the path
  // WARNING: as implemeted, assumes the path crosses each face just once

  // Get polylines for the paths
  std::vector<std::vector<SurfacePoint>> polylines =
      edgeNetwork->getPathPolyline();

  // Populate with vertices on the new split mesh
  std::vector<std::vector<Vertex>> pathVerts;

  // First pass: do all edge splits
  for (auto &poly : polylines) {
    pathVerts.emplace_back();
    std::vector<Vertex> &thisPathVerts = pathVerts.back();
    for (SurfacePoint p : poly) {
      switch (p.type) {
      case SurfacePointType::Vertex: {
        thisPathVerts.push_back(p.vertex);
        break;
      }
      case SurfacePointType::Edge: {
        // Insert new edge

        double tSplit = p.tEdge;
        Halfedge heSplit = p.edge.halfedge();
        Vector3 newPos =
            (1. - tSplit) *
                geometry->inputVertexPositions[heSplit.tailVertex()] +
            tSplit * geometry->inputVertexPositions[heSplit.tipVertex()];

        std::cout << "inserting vertex at " << p << std::endl;

        Vertex newV = mesh->insertVertexAlongEdge(heSplit.edge()).vertex();
        geometry->inputVertexPositions[newV] = newPos;

        thisPathVerts.push_back(newV);

        break;
      }
      case SurfacePointType::Face: {
        throw std::runtime_error("unexpected face point");
        break;
      }
      }
    }
  }

  // Copy marked vertices
  VertexData<char> extraMarkedVertices(*mesh, false);
  for (Vertex v : edgeNetwork->tri->intrinsicMesh->vertices()) {
    if (edgeNetwork->isMarkedVertex[v] &&
        edgeNetwork->tri->vertexLocations[v].type == SurfacePointType::Vertex) {
      extraMarkedVertices[edgeNetwork->tri->vertexLocations[v].vertex] = true;
    }
  }

  // Second pass: connect up vertices with edges
  EdgeData<char> edgeInPath(*mesh, false);
  for (std::vector<Vertex> &thisPathVerts : pathVerts) {

    for (size_t iV = 1; iV < thisPathVerts.size(); iV++) {

      Vertex vA = thisPathVerts[iV - 1];
      Vertex vB = thisPathVerts[iV];

      // std::cout << "Finding edge between " << vA << " " << vB << std::endl;

      // Check if they're already connected
      Edge connectingEdge = mesh->connectingEdge(vA, vB);

      // std::cout << "  initial edge is " << connectingEdge << std::endl;

      if (connectingEdge == Edge()) {
        // Not connected

        // Find the halfedges in a shared face
        Halfedge heA, heB;
        for (Halfedge heACand : vA.outgoingHalfedges()) {
          for (Halfedge heBCand : vB.outgoingHalfedges()) {
            if (heACand.face() == heBCand.face()) {
              heA = heACand;
              heB = heBCand;
            }
          }
        }

        // Connect with a  new edge
        connectingEdge = mesh->connectVertices(heA, heB).edge();

        // std::cout << "  new edge is " << connectingEdge << std::endl;
      }

      edgeInPath[connectingEdge] = true;
    }
  }

  std::cout << "mesh has " << mesh->nVertices() << " verts\n";

  // Triangulate any remaining faces
  for (Face f : mesh->faces()) {
    mesh->triangulate(f);
  }

  geometry->refreshQuantities();

  { // TODO TEMP wiggle the vertices and output
    std::cout << "WIGGLING BAKED MESH\n";
    geometry->requireVertexNormals();
    geometry->requireShapeLengthScale();
    double wiggleScale = geometry->shapeLengthScale * 1e-4;
    for (Vertex v : mesh->vertices()) {
      geometry->inputVertexPositions[v] +=
          geometry->vertexNormals[v] * (unitRand() - 0.5) * wiggleScale;
    }
    geometry->refreshQuantities();

    writeSurfaceMesh(*mesh, *geometry, "baked_wiggled_mesh.obj");
  }

  // Update visualizations
  psMesh = polyscope::registerSurfaceMesh(
      "input mesh", geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  // rebuild a new matching path network
  edgeNetwork = EdgeNetwork::constructFromEdgeSet(*mesh, *geometry, edgeInPath,
                                                  extraMarkedVertices);
  if (edgeNetwork == nullptr) {
    polyscope::warning("could not initialize path from random plane");
    return;
  }
  edgeNetwork->posGeom = geometry.get();

  updatePathViz();
}

// ====== General viz

void clearData() {
  edgeNetwork.reset();
  currMode = FlipMode::None;
}

void makePathTubes() {
  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }

  edgeNetwork->savePathRibbons("", tubeWidth, ribbonPerEdge);

  if (ribbonAllEdge) {
    edgeNetwork->saveEdgeRibbons("", tubeWidth);
  }
}

void exportPathLines() {
  if (edgeNetwork == nullptr) {
    polyscope::warning("no path network");
    return;
  }

  edgeNetwork->savePathOBJLine("", ribbonAllEdge);
}

// Fancy path construction
bool fancyPathClosed = false;
std::vector<Vertex> fancyPathVerts;
std::vector<std::pair<size_t, int>> fancyPathVertsPs;
VertexData<double> fancyPathVertexNumbers;
bool fancyPathMarkVerts = false;
void buildFancyPathUI() {

  auto updateFancyPathViz = [&]() {
    psMesh->addVertexCountQuantity("fancy path vertices", fancyPathVertsPs);
  };

  if (ImGui::Button("Push Vertex")) {

    long long int iV = psMesh->selectVertex();
    if (iV != -1) {
      Vertex v = mesh->vertex(iV);
      fancyPathVerts.push_back(v);
      fancyPathVertsPs.emplace_back((size_t)iV, (int)fancyPathVertsPs.size());
      updateFancyPathViz();
    }
  }
  ImGui::SameLine();
  if (ImGui::Button("Pop Vertex")) {
    if (!fancyPathVerts.empty()) {
      fancyPathVerts.pop_back();
      fancyPathVertsPs.pop_back();
    }
    updateFancyPathViz();
  }

  if (ImGui::Button("New Path From These Points")) {
    edgeNetwork = EdgeNetwork::constructFromPiecewiseDijkstraPath(
        *mesh, *geometry, fancyPathVerts, fancyPathClosed, fancyPathMarkVerts);
    if (edgeNetwork == nullptr) {
      polyscope::warning(
          "could not initialize fancy edge path between vertices");
      return;
    }
    edgeNetwork->posGeom = geometry.get();
    currMode = FlipMode::Path;

    updatePathViz();
  }
  ImGui::Checkbox("Create Closed Path", &fancyPathClosed);
  ImGui::Checkbox("Mark interior vertices", &fancyPathMarkVerts);
}

// A user-defined callback, for creating control panels (etc)
void myCallback() {

  if (ImGui::Button("Load Last Path")) {
    clearData();
    createPathFromLast();
  }
  ImGui::SameLine();
  if (ImGui::Button("New Dijkstra Path")) {
    clearData();
    createPathFromPoints();
  }
  if (ImGui::Button("Random Loop")) {
    clearData();
    createRandomLoop();
  }

  // if (ImGui::Button("New Source Geodesics")) {
  // clearData();
  //}

  if (ImGui::Button("Clear Data")) {
    clearData();
  }

  if (ImGui::TreeNode("Construct Fancy Path")) {
    buildFancyPathUI();
    ImGui::TreePop();
  }

  if (ImGui::TreeNode("Specialty loaders")) {
    if (ImGui::Button("Load edge set")) {
      clearData();
      createPathFromEdgeSet();
    }
    ImGui::SameLine();
    if (ImGui::Button("Load line list obj")) {
      clearData();
      createPathFromObjLines();
    }
    if (ImGui::Button("Load Dijkstra list")) {
      clearData();
      createPathFromDijkstraList();
    }

    if (ImGui::Button("Load UV cut")) {
      clearData();
      createPathFromUVCut();
    }
    ImGui::SameLine();
    if (ImGui::Button("Load seg cut")) {
      clearData();
      createPathFromSeg();
    }
    ImGui::TreePop();
  }

  ImGui::PushItemWidth(150);

  ImGui::Separator();

  switch (currMode) {
  case FlipMode::None: {
    break;
  }
  case FlipMode::Path: {

    break;
  }
  case FlipMode::Loop: {
    break;
  }
  case FlipMode::Source: {
    break;
  }
  }

  // Straightening
  if (ImGui::Button("Locally Shorten")) {
    locallyShorten();
    updatePathViz();
  }
  ImGui::SameLine();
  ImGui::Checkbox("Limit Iteration Count", &iterativeShortenUseIterationCap);
  if (iterativeShortenUseIterationCap) {
    ImGui::SameLine();
    ImGui::InputInt("Iters: ", &iterativeShortenIterationCap);
  }

  ImGui::Checkbox("Limit Length Decrease", &useIterativeShortenLengthLim);
  if (useIterativeShortenLengthLim) {
    ImGui::SameLine();
    ImGui::InputFloat("Relative lim: ", &iterativeShortenLengthLim);
  }

  ImGui::Checkbox("Maintain Delaunay", &maintainDelaunay);
  ImGui::SameLine();
  ImGui::InputFloat("Angle EPS", &angleEPS);
  ImGui::Checkbox("Straighten at marked", &straightenAtMarked);

  if (ImGui::Button("Make Delaunay")) {
    makeDelaunay();
    updatePathViz();
  }
  ImGui::SameLine();
  if (ImGui::Button("Check Path")) {
    checkPath();
  }

  // ====  Extras
  if (ImGui::TreeNode("Extras")) {

    if (ImGui::Button("Flatten")) {
      flattenCutSurface();
    }

    if (ImGui::Button("Distance to curve")) {
      distanceToCurve();
    }

    if (ImGui::Button("Poisson from curve")) {
      poissonFromCurve();
    }

    ImGui::SameLine();
    if (ImGui::Button("Cross field from curve")) {
      crossFieldFromCurve();
    }

    if (ImGui::Button("Bezier subd")) {
      bezierSubdivide();
    }
    ImGui::SameLine();
    ImGui::InputInt("rounds", &nBezierIters);

    ImGui::InputInt("Max insert", &maxInsertions);
    ImGui::InputFloat("Area thresh", &refineAreaThresh);
    ImGui::InputFloat("Angle thresh", &refineAngleThresh);
    if (ImGui::Button("Delaunay refine")) {
      delaunayRefine();
    }
    if (ImGui::Button("Split bent edges")) {
      splitBentEdges();
    }
    ImGui::SameLine();
    ImGui::InputFloat("Bent thresh", &splitAngleDeg);

    if (ImGui::Button("Bake path")) {
      bakePathToMesh();
    }

    ImGui::TreePop();
  }

  // ==== Visualization
  ImGui::Separator();
  if (ImGui::TreeNode("Visualization")) {
    if (edgeNetwork) {
      if (ImGui::Button("Update")) {
        updatePathViz();
      }

      ImGui::Checkbox("Show Intrinsic Edges", &vizAllIntrinsicEdges);

      if (ImGui::Button("Export path lines")) {
        exportPathLines();
      }

      ImGui::SliderFloat("Tube Width", &tubeWidth, 0.,
                         0.1 * polyscope::state::lengthScale);
      ImGui::Checkbox("Ribbons per edge", &ribbonPerEdge);
      ImGui::Checkbox("Ribbons all edges", &ribbonAllEdge);
      if (ImGui::Button("Make Path Ribbons")) {
        makePathTubes();
      }

      ImGui::Checkbox("Make movie", &edgeNetwork->makeMovie);
      if (edgeNetwork->makeMovie) {
        ImGui::Checkbox("Include flips in movie", &edgeNetwork->flipMovie);
      }

      ImGui::InputInt("Movie frame every", &edgeNetwork->movieFrameEvery);

      if (ImGui::Button("Reset screenshot ind")) {
        polyscope::resetScreenshotIndex();
      }

      if (ImGui::Button("Save marked points")) {
        edgeNetwork->saveMarkedPointsAsCloud("marked_points.obj");
      }
    }
    ImGui::TreePop();
  }

  ImGui::PopItemWidth();
}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("Flip edges to find geodesic paths.");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      "input mesh", geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  tubeWidth = polyscope::state::lengthScale / 500;

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
