#include "Dual_contouring_octree_3.h"
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
    OctreeWrapper octree_wrap( 4, { -1, -1, -1, 1, 1, 1 } );
    octree_wrap.print( "../octree.off" );

    CGAL::Octree_oracle octree_oracle( octree_wrap );

    std::cout << "Init grid" << std::endl;

    const std::size_t n_vertices = octree_oracle.n_vertices();

    //#pragma omp parallel for
    for( int i = 0; i < n_vertices; ++i ) {
        const auto& v = octree_oracle.vertices()[i];
        Point_3 p     = octree_oracle.position( v );
        const auto val         = sphere_function( p );
        Vector_3 gradient = p - CGAL::ORIGIN;
        gradient               = gradient / std::sqrt( gradient.squared_length() );
        octree_wrap.value( v ) = val;
        octree_wrap.gradient( v ) = gradient;
    }

    Point_range points;
    Polygon_range polygons;

    std::cout << "Run DC" << std::endl;
    CGAL::make_quad_mesh_using_dual_contouring( octree_oracle, 0.8, points, polygons, false );

    // Mesh mesh_output;
    // CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points, polygons, mesh_output );

    CGAL::IO::write_OFF( "../result.off", points, polygons );

    // write_off("result.off", result);
}
