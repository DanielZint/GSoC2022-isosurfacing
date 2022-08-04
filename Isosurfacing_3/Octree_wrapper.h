#pragma once

#include "types.h"
#include "Octree_tables.h"

#include <CGAL/Octree.h>
#include <CGAL/Orthtree/Split_predicates.h>    // TODO just for implementation. Delete include later.

typedef CGAL::Octree<Kernel, std::vector<Point_3>> Octree;
typedef Octree::Node::Global_coordinates Octree_uniform_coords;                               // coordinates on max depth level
typedef std::tuple<Octree::Node::Global_coordinates, std::size_t> Octree_coords_and_depth;    // global coordinates + depth
typedef std::tuple<std::size_t, std::size_t> Octree_edge_index;                               // edge uniform index + depth

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

        auto leaf_coords = global_leaf_coordinates( n );

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

class OctreeWrapper {
    /*
     * Naming convention from "A parallel dual marching cubes approach to quad only surface reconstruction - Grosso & Zint"
     *
     *        ^ y
     *        |
     *       v2------e2------v3
     *       /|             /|
     *     e11|           e10|
     *     /  e3          /  e1
     *   v6------e6------v7  |
     *    |   |          |   |
     *    |  v0------e0--|---v1 --> x
     *    e7 /           e5 /
     *    | e8           | e9
     *    |/             |/
     *   v4------e4------v5
     *   /
     *  < z
     */
    std::size_t max_depth_ = 0;

    FT offset_x_ = 0;
    FT offset_y_ = 0;
    FT offset_z_ = 0;

    std::size_t dim_;

    FT hx_ = 0;

    Octree octree_;

    std::set<Octree_uniform_coords> leaf_node_max_depth_coordinates_;
    std::set<Octree_edge_index> leaf_edges_;
    std::map<Octree_uniform_coords, FT> vertex_values_;

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
        octree_.refine( Split_by_ratio( 2, max_depth_ ) );

        for( Octree::Node node: octree_.traverse<CGAL::Orthtrees::Leaves_traversal>() ) {
            // write all leaf nodes in a set
            leaf_node_max_depth_coordinates_.insert( uniform_coordinates( node ) );

            // init node values
            vertex_values_[vertex_uniform_coordinates( node, 0b000 )] = 0;
            vertex_values_[vertex_uniform_coordinates( node, 0b001 )] = 0;
            vertex_values_[vertex_uniform_coordinates( node, 0b010 )] = 0;
            vertex_values_[vertex_uniform_coordinates( node, 0b011 )] = 0;
            vertex_values_[vertex_uniform_coordinates( node, 0b100 )] = 0;
            vertex_values_[vertex_uniform_coordinates( node, 0b101 )] = 0;
            vertex_values_[vertex_uniform_coordinates( node, 0b110 )] = 0;
            vertex_values_[vertex_uniform_coordinates( node, 0b111 )] = 0;

            // write all leaf edges in a set
            const auto& coords_global  = node.global_coordinates();
            const auto& coords_uniform = uniform_coordinates( node );
            const auto& depth          = node.depth();
            const auto& df             = depth_factor( node );
            for( const auto& edge_voxels: Tables::edge_to_voxel_neighbor ) {
                bool are_all_nodes_leafs = true;
                for( const auto& node_ijk: edge_voxels ) {
                    const std::size_t x = coords_uniform[0] + df * node_ijk[0];
                    const std::size_t y = coords_uniform[1] + df * node_ijk[1];
                    const std::size_t z = coords_uniform[2] + df * node_ijk[2];
                    // check for overflow / ignore edges on boundary
                    if( x >= dim_ || y >= dim_ || z >= dim_ ) {
                        are_all_nodes_leafs = false;
                        break;
                    }

                    const Octree::Node n = get_node( x, y, z );
                    if( n.depth() > depth ) {
                        are_all_nodes_leafs = false;
                        break;
                    }
                }
                if( are_all_nodes_leafs ) {
                    // add to leaf edge set
                    std::size_t e_gl = e_glIndex( edge_voxels[0][3], coords_global[0], coords_global[1], coords_global[2], depth );
                    if( e_gl == 17 ) {
                        std::cout << "DEBUG " << __LINE__ << std::endl;
                    }
                    leaf_edges_.insert( { e_gl, depth } );
                }
            }
        }
    }

    std::size_t dim() const { return dim_; }
    FT hx() const { return hx_; }
    FT offset_x() const { return offset_x_; }
    FT offset_y() const { return offset_y_; }
    FT offset_z() const { return offset_z_; }

    std::size_t depth_factor( const std::size_t& depth ) const { return std::size_t( 1 ) << ( max_depth_ - depth ); }

    std::size_t depth_factor( const Octree::Node& node ) const { return std::size_t( 1 ) << ( max_depth_ - node.depth() ); }

    Octree_uniform_coords uniform_coordinates( const Octree::Node& node ) const  {
        auto coords          = node.global_coordinates();
        const std::size_t df = depth_factor( node );
        for( int i = 0; i < Octree::Node::Dimension::value; ++i ) {
            coords[i] *= df;
        }

        return coords;
    }

    std::array<Point_3, 8> node_points( const Octree::Node& node ) const {
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

    Point_3 point( const Octree_uniform_coords& vertex_coordinates ) const {
        const FT x0 = offset_x_ + vertex_coordinates[0] * hx_;
        const FT y0 = offset_y_ + vertex_coordinates[1] * hx_;
        const FT z0 = offset_z_ + vertex_coordinates[2] * hx_;
        return { x0, y0, z0 };
    }

    Octree_uniform_coords vertex_uniform_coordinates( const Octree::Node& node, const Octree::Node::Local_coordinates local_coords ) const {
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

    Octree::Node get_node( const std::size_t& i, const std::size_t& j, const std::size_t& k ) const {
        Octree::Node node    = octree_.root();
        const std::size_t& x = i;
        const std::size_t& y = j;
        const std::size_t& z = k;
        while( !node.is_leaf() ) {
            std::size_t dist_to_max = max_depth_ - node.depth() - 1;
            Octree::Node::Local_coordinates loc;
            if( x & ( std::size_t( 1 ) << dist_to_max ) ) {
                loc[0] = true;
            }
            if( y & ( std::size_t( 1 ) << dist_to_max ) ) {
                loc[1] = true;
            }
            if( z & ( std::size_t( 1 ) << dist_to_max ) ) {
                loc[2] = true;
            }
            node = node[loc.to_ulong()];
        }
        return node;
    }

    Octree::Node get_node( const Octree_coords_and_depth& cad ) const {
        Octree::Node node           = octree_.root();
        const auto& [coords, depth] = cad;
        const std::size_t& x        = coords[0];
        const std::size_t& y        = coords[1];
        const std::size_t& z        = coords[2];
        while( node.depth() != depth ) {
            if( node.is_leaf() ) {
                return node;
            }
            const std::size_t dist_to_max       = depth - node.depth() - 1;
            const std::size_t mask              = std::size_t( 1 ) << dist_to_max;
            Octree::Node::Local_coordinates loc = 0;
            if( x & mask ) {
                loc[0] = true;
            }
            if( y & mask ) {
                loc[1] = true;
            }
            if( z & mask ) {
                loc[2] = true;
            }
            node = node[loc.to_ulong()];
        }
        return node;
    }

    /// <summary>
    /// Check if the cell (i,j,k) exists in the leaf node set
    /// </summary>
    /// <param name="i">i-index of cell</param>
    /// <param name="j">j-index of cell</param>
    /// <param name="k">k-index of cell</param>
    /// <returns></returns>
    bool exists( const std::size_t& i, const std::size_t& j, const std::size_t& k ) const {
        Octree::Node::Global_coordinates coords = { (Octree_uniform_coords::value_type)i, (Octree_uniform_coords::value_type)j,
                                                    (Octree_uniform_coords::value_type)k };
        if( leaf_node_max_depth_coordinates_.find( coords ) != leaf_node_max_depth_coordinates_.end() ) {
            return true;
        } else {
            return false;
        }
    }

    std::size_t lex_index( const std::size_t& i, const std::size_t& j, const std::size_t& k, const std::size_t& depth ) const {
        std::size_t dim = std::size_t( 1 ) << depth;
        return k * dim * dim + j * dim + i;
    }

    std::size_t i_index( const std::size_t& lex_index, const std::size_t& depth ) const {
        std::size_t dim = std::size_t( 1 ) << depth;
        return lex_index % dim;
    }

    std::size_t j_index( const std::size_t& lex_index, const std::size_t& depth ) const {
        std::size_t dim = std::size_t( 1 ) << depth;
        return ( ( lex_index / dim ) % dim );
    }

    std::size_t k_index( const std::size_t& lex_index, const std::size_t& depth ) const {
        std::size_t dim = std::size_t( 1 ) << depth;
        return ( lex_index / ( dim * dim ) );
    }

    std::tuple<std::size_t, std::size_t, std::size_t> ijk_index( const std::size_t& lex_index, const std::size_t& depth ) const {
        return std::make_tuple( i_index( lex_index, depth ), j_index( lex_index, depth ), k_index( lex_index, depth ) );
    }

    /// <summary>
    /// compute unique edge global index.
    /// </summary>
    /// <param name="e">local edge index</param>
    /// <param name="i_idx">i-index of cell</param>
    /// <param name="j_idx">j-index of cell</param>
    /// <param name="k_idx">k-index of cell</param>
    /// <returns></returns>
    std::size_t e_glIndex( const std::size_t& e, const std::size_t& i_idx, const std::size_t& j_idx, const std::size_t& k_idx,
                           const std::size_t& depth ) const {
        const unsigned long long gei_pattern_ = 670526590282893600ull;
        const size_t i                        = i_idx + (size_t)( ( gei_pattern_ >> 5 * e ) & 1 );            // global_edge_id[eg][0];
        const size_t j                        = j_idx + (size_t)( ( gei_pattern_ >> ( 5 * e + 1 ) ) & 1 );    // global_edge_id[eg][1];
        const size_t k                        = k_idx + (size_t)( ( gei_pattern_ >> ( 5 * e + 2 ) ) & 1 );    // global_edge_id[eg][2];
        const size_t offs                     = (size_t)( ( gei_pattern_ >> ( 5 * e + 3 ) ) & 3 );
        return ( 3 * lex_index( i, j, k, depth ) + offs );
    }

    bool is_leaf_edge( const std::size_t& e_global_index, const std::size_t& depth ) {
        if( leaf_edges_.find( { e_global_index, depth } ) != leaf_edges_.end() ) {
            return true;
        } else {
            return false;
        }
    }

    const std::set<Octree_edge_index>& leaf_edges() const { return leaf_edges_; }

    FT value( const std::size_t i, const std::size_t j, const std::size_t k ) const {
        Octree_uniform_coords v_id = { (Octree_uniform_coords::value_type)i, (Octree_uniform_coords::value_type)j,
                                       (Octree_uniform_coords::value_type)k };
        return vertex_values_.at( v_id );
    }

    FT& value( const std::size_t i, const std::size_t j, const std::size_t k ) {
        Octree_uniform_coords v_id = { (Octree_uniform_coords::value_type)i, (Octree_uniform_coords::value_type)j,
                                       (Octree_uniform_coords::value_type)k };
        return vertex_values_[v_id];
    }

    std::array<FT, 8> voxel_values( const std::size_t i, const std::size_t j, const std::size_t k ) const {
        Octree::Node node = get_node( i, j, k );

        std::array<Octree_uniform_coords, 8> v;
        v[0] = vertex_uniform_coordinates( node, 0b000 );
        v[1] = vertex_uniform_coordinates( node, 0b001 );
        v[2] = vertex_uniform_coordinates( node, 0b010 );
        v[3] = vertex_uniform_coordinates( node, 0b011 );
        v[4] = vertex_uniform_coordinates( node, 0b100 );
        v[5] = vertex_uniform_coordinates( node, 0b101 );
        v[6] = vertex_uniform_coordinates( node, 0b110 );
        v[7] = vertex_uniform_coordinates( node, 0b111 );

        std::array<FT, 8> s;
        s[0] = vertex_values_.at( v[0] );
        s[1] = vertex_values_.at( v[1] );
        s[2] = vertex_values_.at( v[2] );
        s[3] = vertex_values_.at( v[3] );
        s[4] = vertex_values_.at( v[4] );
        s[5] = vertex_values_.at( v[5] );
        s[6] = vertex_values_.at( v[6] );
        s[7] = vertex_values_.at( v[7] );

        return s;
    }

    std::array<FT, 2> edge_values( const Octree_edge_index& e_id ) const {
        const auto& [e_global_id, depth] = e_id;
        const auto df                    = depth_factor( depth );

        const size_t v0_lex_index = e_global_id / 3;
        const auto [i0, j0, k0]   = ijk_index( v0_lex_index, depth );

        // v1
        const std::size_t e_local_index = Tables::edge_store_index[e_global_id % 3];
        const auto& v1_local            = Tables::local_vertex_position[Tables::edge_to_vertex[e_local_index][1]];

        const std::size_t i1 = i0 + v1_local[0];
        const std::size_t j1 = j0 + v1_local[1];
        const std::size_t k1 = k0 + v1_local[2];

        return { value( df * i0, df * j0, df * k0 ), value( df * i1, df * j1, df * k1 ) };
    }

    /// <summary>
    /// Get the 4 voxels incident to an edge. If an edge has only three incident voxels, one will appear twice. The voxels are given with the uniform lexicographical index.
    /// </summary>
    /// <param name="e_id"></param>
    /// <returns></returns>
    std::array<std::size_t, 4> voxels_incident_to_edge( const Octree_edge_index& e_id ) const {
        const auto& [e_global_id, depth] = e_id;
        const std::size_t e_local_index  = Tables::edge_store_index[e_global_id % 3];

        const auto df = depth_factor( depth );

        const size_t v0_lex_index = e_global_id / 3;
        auto [i, j, k]   = ijk_index( v0_lex_index, depth );
        i *= df;
        j *= df;
        k *= df;

        const auto& voxel_neighbors = Tables::edge_to_voxel_neighbor[e_local_index];
        Octree::Node n0 = get_node(i + voxel_neighbors[0][0],j + voxel_neighbors[0][1],k + voxel_neighbors[0][2]);
        Octree::Node n1 = get_node(i + voxel_neighbors[1][0],j + voxel_neighbors[1][1],k + voxel_neighbors[1][2]);
        Octree::Node n2 = get_node(i + voxel_neighbors[2][0],j + voxel_neighbors[2][1],k + voxel_neighbors[2][2]);
        Octree::Node n3 = get_node(i + voxel_neighbors[3][0],j + voxel_neighbors[3][1],k + voxel_neighbors[3][2]);

        const Octree_uniform_coords n0_uniform_coords = uniform_coordinates(n0);
        const Octree_uniform_coords n1_uniform_coords = uniform_coordinates(n1);
        const Octree_uniform_coords n2_uniform_coords = uniform_coordinates(n2);
        const Octree_uniform_coords n3_uniform_coords = uniform_coordinates(n3);

        std::size_t n0_lex = lex_index(n0_uniform_coords[0],n0_uniform_coords[1],n0_uniform_coords[2],max_depth_);
        std::size_t n1_lex = lex_index(n1_uniform_coords[0],n1_uniform_coords[1],n1_uniform_coords[2],max_depth_);
        std::size_t n2_lex = lex_index(n2_uniform_coords[0],n2_uniform_coords[1],n2_uniform_coords[2],max_depth_);
        std::size_t n3_lex = lex_index(n3_uniform_coords[0],n3_uniform_coords[1],n3_uniform_coords[2],max_depth_);

        return {n0_lex, n1_lex, n2_lex, n3_lex};

        //return { value( i0, j0, k0 ), value( i1, j1, k1 ) };
    }

    std::array<Point_3, 8> voxel_vertex_positions( const std::size_t i, const std::size_t j, const std::size_t k ) const {
        Octree::Node node = get_node( i, j, k );
        return node_points( node );
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
            p[0] = point( vertex_uniform_coordinates( node, 0b000 ) );
            p[1] = point( vertex_uniform_coordinates( node, 0b001 ) );
            p[2] = point( vertex_uniform_coordinates( node, 0b010 ) );
            p[3] = point( vertex_uniform_coordinates( node, 0b011 ) );
            p[4] = point( vertex_uniform_coordinates( node, 0b100 ) );
            p[5] = point( vertex_uniform_coordinates( node, 0b101 ) );
            p[6] = point( vertex_uniform_coordinates( node, 0b110 ) );
            p[7] = point( vertex_uniform_coordinates( node, 0b111 ) );

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