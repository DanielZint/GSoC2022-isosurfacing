#include "Cartesian_grid_3.h"
#include "Cartesian_grid_oracle.h"
#include "Function_oracle.h"
#include "Marching_cubes_3.h"
#include "types.h"

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <fstream>
#include <iostream>
#include <math.h>

int main() {
    auto sphere_function = []( const Point_3& point ) { return std::sqrt( point.x() * point.x() + point.y() * point.y() + point.z() * point.z() ); };
    CGAL::Function_oracle<Kernel, decltype( sphere_function )> sphere_oracle( sphere_function, { -1, -1, -1, 1, 1, 1 },
                                                                              Vector_3( 0.02f, 0.02f, 0.02f ) );

    Point_range points;
    Polygon_range polygons;

    CGAL::make_triangle_mesh_using_marching_cubes( sphere_oracle, 0.8f, points, polygons );

    // TODO: compare results with mesh_3

    Mesh mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points, polygons, mesh );

    CGAL::IO::write_OFF( "result.off", mesh );

    // write_off("result.off", result);
}
