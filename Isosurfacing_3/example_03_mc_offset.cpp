#include "Cartesian_grid_3.h"
#include "Cartesian_grid_oracle.h"
#include "Function_oracle.h"
#include "Marching_cubes_3.h"
#include "types.h"

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <fstream>
#include <iostream>
#include <math.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

Kernel::FT sphere_function( const Point_3& point ) { return std::sqrt( point.x() * point.x() + point.y() * point.y() + point.z() * point.z() ); };

inline Kernel::FT distance_to_mesh( const Tree& tree, const Point_3& p ) {
    const Point_3& x = tree.closest_point( p );
    return std::sqrt( ( p - x ).squared_length() );
}

int main() {
    const std::string input_name = "D:/Documents/Projects/GSoC2022-isosurfacing-fork/bunny.off";
    const int n_voxels           = 50;
    const FT offset_value        = 0.02;

    Mesh mesh_input;
    if( !CGAL::IO::read_OFF( input_name, mesh_input ) ) {
        std::cout << "Could not read mesh" << std::endl;
        exit( -1 );
    }

    // compute bounding box
    CGAL::Bbox_3 aabb_grid     = PMP::bbox( mesh_input );
    Vector_3 aabb_increase_vec = Vector_3( offset_value * 1.1, offset_value * 1.1, offset_value * 1.1 );
    aabb_grid += ( Point_3( aabb_grid.xmax(), aabb_grid.ymax(), aabb_grid.zmax() ) + aabb_increase_vec ).bbox();
    aabb_grid += ( Point_3( aabb_grid.xmin(), aabb_grid.ymin(), aabb_grid.zmin() ) - aabb_increase_vec ).bbox();

    // construct AABB tree
    Tree tree( mesh_input.faces_begin(), mesh_input.faces_end(), mesh_input );

    CGAL::Side_of_triangle_mesh<Mesh, CGAL::GetGeomTraits<Mesh>::type> sotm( mesh_input );

    Grid grid( n_voxels, n_voxels, n_voxels, aabb_grid );

    CGAL::Cartesian_grid_oracle<Kernel> grid_oracle( grid );

    std::cout << "Init grid" << std::endl;

    const std::size_t size_k = grid.xdim();
    const std::size_t size_j = grid.ydim();
    const std::size_t size_i = grid.zdim();

    #pragma omp parallel for
    for( int z = 0; z < grid.zdim(); z++ ) {
        for( int y = 0; y < grid.ydim(); y++ ) {
            for( int x = 0; x < grid.xdim(); x++ ) {
                const auto& p         = grid_oracle.position( x, y, z );
                grid.value( x, y, z ) = distance_to_mesh( tree, p );
                // const bool is_inside  = ( sotm( p ) == CGAL::ON_BOUNDED_SIDE );
                // if( is_inside ) {
                //    grid.value( x, y, z ) *= -1;
                //}
            }
        }
    }

    Point_range points;
    Polygon_range polygons;

    std::cout << "Run MC" << std::endl;
    CGAL::make_triangle_mesh_using_marching_cubes( grid_oracle, offset_value, points, polygons );

    // TODO: compare results with mesh_3

    Mesh mesh_output;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points, polygons, mesh_output );

    CGAL::IO::write_OFF( "result.off", mesh_output );

    // write_off("result.off", result);
}
