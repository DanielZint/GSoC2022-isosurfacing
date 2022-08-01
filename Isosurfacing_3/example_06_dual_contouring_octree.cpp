#include "Cartesian_grid_3.h"
#include "Cartesian_grid_oracle.h"
#include "Dual_contouring_3.h"
#include "Octree_oracle.h"
#include "Octree_wrapper.h"
#include "types.h"

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <fstream>
#include <iostream>
#include <math.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

Kernel::FT sphere_function( const Point_3& point ) { return std::sqrt( point.x() * point.x() + point.y() * point.y() + point.z() * point.z() ); };

int main() {
    const int n_voxels = 11;

    Grid grid( n_voxels, n_voxels, n_voxels, { -1, -1, -1, 1, 1, 1 } );

    OctreeWrapper octree_wrap( 2, { -1, -1, -1, 1, 1, 1 } );
    octree_wrap.print( "../octree.off" );

    // CGAL::Cartesian_grid_oracle<Kernel> grid_oracle( grid );
    CGAL::Octree_oracle octree_oracle( octree_wrap );

    std::cout << "Init grid" << std::endl;

    const std::size_t size_k = grid.xdim();
    const std::size_t size_j = grid.ydim();
    const std::size_t size_i = grid.zdim();

    //#pragma omp parallel for
    for( int z = 0; z < octree_wrap.dim() + 1; z++ ) {
        for( int y = 0; y < octree_wrap.dim() + 1; y++ ) {
            for( int x = 0; x < octree_wrap.dim() + 1; x++ ) {
                const auto& p                  = octree_oracle.position( x, y, z );
                octree_wrap.value( x, y, z ) = sphere_function( p );
            }
        }
    }

    Point_range points;
    Polygon_range polygons;

    std::cout << "Run DC" << std::endl;
    CGAL::make_quad_mesh_using_dual_contouring( octree_oracle, 0.8, points, polygons, true );

    // Mesh mesh_output;
    // CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points, polygons, mesh_output );

    CGAL::IO::write_OFF( "../result.off", points, polygons );

    // write_off("result.off", result);
}
