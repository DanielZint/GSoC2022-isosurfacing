#pragma once

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "Cartesian_grid_3.h"

typedef CGAL::Simple_cartesian<float> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Vector_3 Vector_3;
typedef typename Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor Vertex_descriptor;

typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point_3> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;