#include "CLI11.hpp"
#include "Cartesian_grid_3.h"
#include "Cartesian_grid_oracle.h"
#include "Dual_contouring_3.h"
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

inline Kernel::FT distance_to_mesh( const Tree& tree, const Point_3& p ) {
    const Point_3& x = tree.closest_point( p );
    return std::sqrt( ( p - x ).squared_length() );
}

int main( int argc, char** argv ) {
    CLI::App app;
    app.description( "This program computes the offset to a mesh using Marching Cubes or Dual Contouring. To see the options use --help." );

    bool use_mc        = false;
    FT offset_value    = 0;
    int n_voxel_points = 3;
    std::string input_name;
    std::string output_name;
    int inside_outside = 0;
    bool use_bbox      = false;
    app.add_flag( "--mc,!--dc", use_mc, "Use Marching Cubes (--mc) or Dual Contouring (--dc) for offset computation" )->required();
    app.add_option( "--offset", offset_value, "Value by which the input mesh is offsetted" )->required()->check( CLI::NonNegativeNumber );
    app.add_option( "--np", n_voxel_points, "Number of sampling points in every dimension" )->required()->check( CLI::PositiveNumber );
    app.add_option( "--input,-i", input_name, "Input mesh in off file format" )->required()->check( CLI::ExistingFile );
    app.add_option( "--output,-o", output_name, "Output mesh in off file format" )->default_val( "offset.off" );
    app.add_flag( "--inside{1},--outside{2},--doublesided{0},--is{1},--os{2},--ds{0}", inside_outside,
                  "Compute the offset on the inside, outside or both sides of the mesh" )
        ->default_val( 0 );
    app.add_flag( "--bbox", use_bbox, "Make sure that vertices stay inside their voxel when using Dual Contouring" )->default_val( false );

    CLI11_PARSE( app, argc, argv );

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

    Grid grid( n_voxel_points, n_voxel_points, n_voxel_points, aabb_grid );

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
            }
        }
    }

    // add sign if output is not double sided
    if( inside_outside != 0 ) {
#pragma omp parallel for
        for( int z = 0; z < grid.zdim(); z++ ) {
            for( int y = 0; y < grid.ydim(); y++ ) {
                for( int x = 0; x < grid.xdim(); x++ ) {
                    const auto& p        = grid_oracle.position( x, y, z );
                    const bool is_inside = ( sotm( p ) == CGAL::ON_BOUNDED_SIDE );
                    if( is_inside && inside_outside == 2 ) {
                        grid.value( x, y, z ) *= -1;
                    } else if( !is_inside && inside_outside == 1 ) {
                        grid.value( x, y, z ) *= -1;
                    }
                }
            }
        }
    }

    Point_range points;
    Polygon_range polygons;

    if( use_mc ) {
        std::cout << "Run Marching Cubes" << std::endl;
        CGAL::make_triangle_mesh_using_marching_cubes( grid_oracle, offset_value, points, polygons );
    } else {
        std::cout << "Run Dual Contouring" << std::endl;
        CGAL::make_quad_mesh_using_dual_contouring( grid_oracle, offset_value, points, polygons, use_bbox );
    }

    CGAL::IO::write_OFF( output_name, points, polygons );

}
