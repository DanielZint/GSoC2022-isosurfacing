#pragma once

#include "types.h"

#include <CGAL/Octree.h>
#include <CGAL/Orthtree/Split_predicates.h>    // TODO just for implementation. Delete include later.

typedef CGAL::Octree<Kernel, std::vector<Point_3>> Octree;

struct Split_by_ratio {
    std::size_t ratio_;
    std::size_t min_depth_;
    std::size_t max_depth_;

    std::size_t octree_dim_;

    Split_by_ratio( std::size_t min_depth, std::size_t max_depth ) : min_depth_( min_depth ), max_depth_( max_depth ) {
        octree_dim_ = std::size_t( 1 ) << max_depth_;
    }
    /*!
    \brief returns `true` if `n` should be split, `false` otherwise.
   */
    // template<class Node>
    bool operator()( const Octree::Node& n ) const {
        // n.depth()
        if( n.depth() < min_depth_ ) {
            return true;
        }
        if( n.depth() == max_depth_ ) {
            return false;
        }

        auto coords = n.global_coordinates();

        if( coords[0] >= octree_dim_ / 4 ) {
            return false;
        }
        if( coords[1] >= octree_dim_ / 4 ) {
            return false;
        }
        if( coords[2] >= octree_dim_ / 4 ) {
            return false;
        }
        return true;
    }
};

class OctreeWrapper {
    std::size_t max_depth_ = 0;

    FT offset_x_ = 0;
    FT offset_y_ = 0;
    FT offset_z_ = 0;

    FT hx_ = 0;

    Octree octree_;

  public:
    OctreeWrapper( const std::size_t max_depth, const CGAL::Bbox_3& bbox )
        : max_depth_( max_depth ), offset_x_( bbox.xmin() ), offset_y_( bbox.ymin() ), offset_z_( bbox.zmin() ),
          octree_( std::vector<Point_3> { { bbox.xmin(), bbox.ymin(), bbox.zmin() }, { bbox.xmax(), bbox.ymax(), bbox.zmax() } } ) {
        if( max_depth_ < 1 ) {
            std::cout << "tree depth must be larger than 0" << std::endl;
            exit( -1 );
        }
        hx_ = bbox.x_span() / ( 1 << max_depth_ );

        // octree_.refine( CGAL::Orthtree:: );
        octree_.refine( Split_by_ratio( 1, max_depth_ ) );

        Mesh m;

        // Print out the octree
        for( Octree::Node node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            std::cout << node << std::endl;
            auto p = node_points( node );

            std::array<Vertex_descriptor, 8> v;
            std::transform( p.begin(), p.end(), v.begin(), [&m]( const Point_3& e ) { return m.add_vertex( e ); } );

            m.add_face( v[0], v[2], v[3], v[1] );
            m.add_face( v[0], v[4], v[6], v[2] );
            m.add_face( v[2], v[6], v[7], v[3] );
            m.add_face( v[1], v[3], v[7], v[5] );
            m.add_face( v[0], v[1], v[5], v[4] );
            m.add_face( v[7], v[6], v[4], v[5] );
        }

        CGAL::write_off( m, "../octree.off" );
    }

    std::array<Point_3, 8> node_points( const Octree::Node& node ) {
        auto coords = node.global_coordinates();
        const std::size_t depth_factor = std::size_t(1) << ( max_depth_ - node.depth() );

        // x0
        FT x0 = offset_x_ + coords[0] * depth_factor * hx_;
        FT y0 = offset_y_ + coords[1] * depth_factor * hx_;
        FT z0 = offset_z_ + coords[2] * depth_factor * hx_;
        FT x1 = offset_x_ + ( coords[0] + 1 ) * depth_factor * hx_;
        FT y1 = offset_y_ + ( coords[1] + 1 ) * depth_factor * hx_;
        FT z1 = offset_z_ + ( coords[2] + 1 ) * depth_factor * hx_;

        std::array<Point_3, 8> points;
        points[0] = { x0, y0, z0 };
        points[1] = { x1, y0, z0 };
        points[2] = { x0, y1, z0 };
        points[3] = { x1, y1, z0 };

        points[4] = { x0, y0, z1 };
        points[5] = { x1, y0, z1 };
        points[6] = { x0, y1, z1 };
        points[7] = { x1, y1, z1 };

        return points;
    }

  private:
};