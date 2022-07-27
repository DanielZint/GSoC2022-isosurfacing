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

Kernel::FT sphere_function( const Point_3& point ) { return std::sqrt( point.x() * point.x() + point.y() * point.y() + point.z() * point.z() ); };

int main() {
    Grid grid( 50, 50, 50, { -1, -1, -1, 1, 1, 1 } );

    CGAL::Cartesian_grid_oracle<Kernel> grid_oracle( grid );

    for( std::size_t x = 0; x < grid.xdim(); x++ ) {
        for( std::size_t y = 0; y < grid.ydim(); y++ ) {
            for( std::size_t z = 0; z < grid.zdim(); z++ ) {
                grid.value( x, y, z ) = sphere_function( grid_oracle.position( x, y, z ) );
            }
        }
    }

    Point_range points;
    Polygon_range polygons;

    CGAL::make_triangle_mesh_using_marching_cubes( grid_oracle, 0.8f, points, polygons );

    // TODO: compare results with mesh_3

    Mesh mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points, polygons, mesh );

    CGAL::IO::write_OFF( "result.off", mesh );

    // write_off("result.off", result);
}
