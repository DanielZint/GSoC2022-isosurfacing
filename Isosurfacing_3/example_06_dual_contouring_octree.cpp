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

struct Refine_one_eighth {
    std::size_t min_depth_;
    std::size_t max_depth_;

    std::size_t octree_dim_;

    Octree_wrapper::Uniform_coords uniform_coordinates( const Octree_wrapper::Octree::Node& node ) const {
        auto coords                    = node.global_coordinates();
        const std::size_t depth_factor = std::size_t( 1 ) << ( max_depth_ - node.depth() );
        for( int i = 0; i < Octree_wrapper::Octree::Node::Dimension::value; ++i ) {
            coords[i] *= depth_factor;
        }

        return coords;
    }

    Refine_one_eighth( std::size_t min_depth, std::size_t max_depth ) : min_depth_( min_depth ), max_depth_( max_depth ) {
        octree_dim_ = std::size_t( 1 ) << max_depth_;
    }

    bool operator()( const Octree_wrapper::Octree::Node& n ) const {
        // n.depth()
        if( n.depth() < min_depth_ ) {
            return true;
        }
        if( n.depth() == max_depth_ ) {
            return false;
        }

        auto leaf_coords = uniform_coordinates( n );

        if( leaf_coords[0] >= octree_dim_ / 2 ) {
            return false;
        }
        if( leaf_coords[1] >= octree_dim_ / 2 ) {
            return false;
        }
        if( leaf_coords[2] >= octree_dim_ / 2 ) {
            return false;
        }
        return true;
    }
};

int main() {
    Octree_wrapper octree_wrap( { -1, -1, -1, 1, 1, 1 } );

    Refine_one_eighth split_predicate( 2, 4 );
    octree_wrap.refine( split_predicate );

    octree_wrap.print( "../octree.off" );

    CGAL::Octree_oracle octree_oracle( octree_wrap );

    std::cout << "Init grid" << std::endl;

    const std::size_t n_vertices = octree_oracle.n_vertices();

    //#pragma omp parallel for
    for( int i = 0; i < n_vertices; ++i ) {
        const auto& v             = octree_oracle.vertices( i );
        Point_3 p                 = octree_oracle.position( v );
        const auto val            = sphere_function( p );
        Vector_3 gradient         = p - CGAL::ORIGIN;
        gradient                  = gradient / std::sqrt( gradient.squared_length() );
        octree_wrap.value( v )    = val;
        octree_wrap.gradient( v ) = gradient;
    }

    Point_range points;
    Polygon_range polygons;

    CGAL::Dual_contouring_3::Positioning::QEM_SVD dc_positioning;

    std::cout << "Run DC" << std::endl;
    CGAL::make_quad_mesh_using_dual_contouring( octree_oracle, 0.8, points, polygons, dc_positioning );

    // Mesh mesh_output;
    // CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points, polygons, mesh_output );

    CGAL::IO::write_OFF( "../result.off", points, polygons );

    // write_off("result.off", result);
}
