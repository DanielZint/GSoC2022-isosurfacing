#pragma once

#include "Octree_edge_index.h"
#include "types.h"

#include <CGAL/Octree.h>
#include <CGAL/Orthtree/Split_predicates.h>    // TODO just for implementation. Delete include later.

typedef CGAL::Octree<Kernel, std::vector<Point_3>> Octree;
typedef Octree::Node::Global_coordinates Octree_vertex_id;

struct Split_by_ratio {
    std::size_t min_depth_;
    std::size_t max_depth_;

    std::size_t octree_dim_;

    Octree::Node::Global_coordinates global_leaf_coordinates( const Octree::Node& node ) const {
        auto coords                    = node.global_coordinates();
        const std::size_t depth_factor = std::size_t( 1 ) << ( max_depth_ - node.depth() );
        for( int i = 0; i < Octree::Node::Dimension::value; ++i ) {
            coords[i] *= depth_factor;
        }

        return coords;
    }

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

        // auto leaf_coords = global_leaf_coordinates( n );
        //
        // if( leaf_coords[0] >= octree_dim_ / 2 ) {
        //    return false;
        //}
        // if( leaf_coords[1] >= octree_dim_ / 2 ) {
        //    return false;
        //}
        // if( leaf_coords[2] >= octree_dim_ / 2 ) {
        //    return false;
        //}
        return true;
    }
};

class OctreeWrapper {
    std::size_t max_depth_ = 0;

    FT offset_x_ = 0;
    FT offset_y_ = 0;
    FT offset_z_ = 0;

    std::size_t dim_;

    FT hx_ = 0;

    Octree octree_;

    // TODO add node values
    std::map<Octree_vertex_id, FT> node_values_;

  public:
    OctreeWrapper( const std::size_t max_depth, const CGAL::Bbox_3& bbox )
        : max_depth_( max_depth ), offset_x_( bbox.xmin() ), offset_y_( bbox.ymin() ), offset_z_( bbox.zmin() ),
          octree_( std::vector<Point_3> { { bbox.xmin(), bbox.ymin(), bbox.zmin() }, { bbox.xmax(), bbox.ymax(), bbox.zmax() } } ) {
        if( max_depth_ < 1 ) {
            std::cout << "tree depth must be larger than 0" << std::endl;
            exit( -1 );
        }
        dim_ = std::size_t( 1 ) << max_depth_;
        hx_  = bbox.x_span() / ( 1 << max_depth_ );

        // octree_.refine( CGAL::Orthtree:: );
        octree_.refine( Split_by_ratio( 1, max_depth_ ) );

        // init node values
        for( Octree::Node node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            node_values_[global_vertex_coordinates( node, 0b000 )] = 0;
            node_values_[global_vertex_coordinates( node, 0b001 )] = 0;
            node_values_[global_vertex_coordinates( node, 0b010 )] = 0;
            node_values_[global_vertex_coordinates( node, 0b011 )] = 0;
            node_values_[global_vertex_coordinates( node, 0b100 )] = 0;
            node_values_[global_vertex_coordinates( node, 0b101 )] = 0;
            node_values_[global_vertex_coordinates( node, 0b110 )] = 0;
            node_values_[global_vertex_coordinates( node, 0b111 )] = 0;
        }
    }

    std::size_t dim() const { return dim_; }
    FT hx() const { return hx_; }
    FT offset_x() const { return offset_x_; }
    FT offset_y() const { return offset_y_; }
    FT offset_z() const { return offset_z_; }

    std::size_t depth_factor( const Octree::Node& node ) const { return std::size_t( 1 ) << ( max_depth_ - node.depth() ); }

    Octree::Node::Global_coordinates global_leaf_coordinates( const Octree::Node& node ) {
        auto coords          = node.global_coordinates();
        const std::size_t df = depth_factor( node );
        for( int i = 0; i < Octree::Node::Dimension::value; ++i ) {
            coords[i] *= df;
        }

        return coords;
    }

    std::array<Point_3, 8> node_points( const Octree::Node& node ) {
        auto coords          = node.global_coordinates();
        const std::size_t df = depth_factor( node );

        const FT x0 = offset_x_ + coords[0] * df * hx_;
        const FT y0 = offset_y_ + coords[1] * df * hx_;
        const FT z0 = offset_z_ + coords[2] * df * hx_;
        const FT x1 = offset_x_ + ( coords[0] + 1 ) * df * hx_;
        const FT y1 = offset_y_ + ( coords[1] + 1 ) * df * hx_;
        const FT z1 = offset_z_ + ( coords[2] + 1 ) * df * hx_;

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

    Point_3 point( const Octree_vertex_id& vertex_coordinates ) const {
        const FT x0 = offset_x_ + vertex_coordinates[0] * hx_;
        const FT y0 = offset_y_ + vertex_coordinates[1] * hx_;
        const FT z0 = offset_z_ + vertex_coordinates[2] * hx_;
        return { x0, y0, z0 };
    }

    Octree_vertex_id global_vertex_coordinates( const Octree::Node& node, const Octree::Node::Local_coordinates local_coords ) const {
        const auto node_coords = node.global_coordinates();
        auto v_coords          = node_coords;
        for( int i = 0; i < Octree::Node::Dimension::value; ++i ) {
            v_coords[i] += std::size_t( local_coords[i] );
        }

        const auto df = depth_factor( node );
        for( int i = 0; i < Octree::Node::Dimension::value; ++i ) {
            v_coords[i] *= df;
        }

        return v_coords;
    }

    std::size_t lex_index( const std::size_t& i, const std::size_t& j, const std::size_t& k ) const { return k * dim_ * dim_ + j * dim_ + i; }

    /// <summary>
    /// compute unique edge global index.
    /// </summary>
    /// <param name="e">local edge index</param>
    /// <param name="i_idx">i-index of cell</param>
    /// <param name="j_idx">j-index of cell</param>
    /// <param name="k_idx">k-index of cell</param>
    /// <returns></returns>
    std::size_t e_glIndex( const std::size_t& e, const std::size_t& i_idx, const std::size_t& j_idx, const std::size_t& k_idx ) const {
        const unsigned long long gei_pattern_ = 670526590282893600ull;
        const size_t i                        = i_idx + (size_t)( ( gei_pattern_ >> 5 * e ) & 1 );            // global_edge_id[eg][0];
        const size_t j                        = j_idx + (size_t)( ( gei_pattern_ >> ( 5 * e + 1 ) ) & 1 );    // global_edge_id[eg][1];
        const size_t k                        = k_idx + (size_t)( ( gei_pattern_ >> ( 5 * e + 2 ) ) & 1 );    // global_edge_id[eg][2];
        const size_t offs                     = (size_t)( ( gei_pattern_ >> ( 5 * e + 3 ) ) & 3 );
        return ( 3 * lex_index( i, j, k ) + offs );
    }

    FT value( const std::size_t i, const std::size_t j, const std::size_t k ) const {
        Octree_vertex_id v_id = { (Octree_vertex_id::value_type)i, (Octree_vertex_id::value_type)j, (Octree_vertex_id::value_type)k };
        return node_values_.at( v_id );
    }

    FT& value( const std::size_t i, const std::size_t j, const std::size_t k ) {
        Octree_vertex_id v_id = { (Octree_vertex_id::value_type)i, (Octree_vertex_id::value_type)j, (Octree_vertex_id::value_type)k };
        return node_values_[v_id];
    }

    void print( const std::string& filename ) const {
        Mesh m;
        for( Octree::Node node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            auto coords = node.global_coordinates();
            // if( coords[1] == 0 && coords[2] == 0 ) {
            //    auto leaf_coords = global_leaf_coordinates( node );
            //    std::cout << node << std::endl;
            //    std::cout << leaf_coords[0] << " " << leaf_coords[1] << " " << leaf_coords[2] << std::endl;
            //}
            std::array<Point_3, 8> p;
            p[0] = point( global_vertex_coordinates( node, 0b000 ) );
            p[1] = point( global_vertex_coordinates( node, 0b001 ) );
            p[2] = point( global_vertex_coordinates( node, 0b010 ) );
            p[3] = point( global_vertex_coordinates( node, 0b011 ) );
            p[4] = point( global_vertex_coordinates( node, 0b100 ) );
            p[5] = point( global_vertex_coordinates( node, 0b101 ) );
            p[6] = point( global_vertex_coordinates( node, 0b110 ) );
            p[7] = point( global_vertex_coordinates( node, 0b111 ) );

            std::array<Vertex_descriptor, 8> v;
            std::transform( p.begin(), p.end(), v.begin(), [&m]( const Point_3& e ) { return m.add_vertex( e ); } );

            m.add_face( v[0], v[2], v[3], v[1] );
            m.add_face( v[0], v[4], v[6], v[2] );
            m.add_face( v[2], v[6], v[7], v[3] );
            m.add_face( v[1], v[3], v[7], v[5] );
            m.add_face( v[0], v[1], v[5], v[4] );
            m.add_face( v[7], v[6], v[4], v[5] );
        }

        CGAL::write_off( m, filename );
    }

  private:
};