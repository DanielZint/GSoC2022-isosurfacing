#ifndef CGAL_OCTREE_GRID_ORACLE_H
#define CGAL_OCTREE_GRID_ORACLE_H

#include "Octree_wrapper.h"
#include "types.h"

#include <array>

namespace CGAL {

    class Octree_oracle {
      public:
        typedef Kernel Geom_traits;
        typedef typename Geom_traits::FT FT;
        typedef typename Geom_traits::Point_3 Point_3;
        typedef typename Geom_traits::Vector_3 Vector_3;

      public:
        Octree_oracle( const OctreeWrapper& octree ) : octree_( &octree ), octree_edges_( octree.leaf_edges().begin(), octree.leaf_edges().end() ) {}

        std::size_t size_x() const { return octree_->dim(); }
        std::size_t size_y() const { return octree_->dim(); }
        std::size_t size_z() const { return octree_->dim(); }

        std::size_t lex_index( const std::size_t& i, const std::size_t& j, const std::size_t& k ) const {
            return k * size_z() * size_y() + j * size_x() + i;
        }

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

        /// compute the gradient of the scalar field with central difference
        void gradient( std::array<Vector_3, 8>& n, const std::array<float, 8>& s, const std::size_t& i, const std::size_t& j,
                       const std::size_t& k ) const {
            auto index = []( const int dim, const int ii ) { return ( ii < 0 ? 0 : ii >= dim ? ( dim - 1 ) : ii ); };

            const FT& dx = octree_->hx();
            const FT& dy = octree_->hx();
            const FT& dz = octree_->hx();

            const std::size_t dim = octree_->dim() + 1;

            // vertex 0
            n[0] = { 0.5f * ( s[1] - value( index( dim, i - 1 ), j, k ) ) / dx, 0.5f * ( s[2] - value( i, index( dim, j - 1 ), k ) ) / dy,
                     0.5f * ( s[4] - value( i, j, index( dim, k - 1 ) ) ) / dz };

            // vertex 1
            n[1] = { 0.5f * ( value( index( dim, i + 2 ), j, k ) - s[0] ) / dx, 0.5f * ( s[3] - value( i + 1, index( dim, j - 1 ), k ) ) / dy,
                     0.5f * ( s[5] - value( i + 1, j, index( dim, k - 1 ) ) ) / dz };

            // vertex 2
            n[2] = { 0.5f * ( s[3] - value( index( dim, i - 1 ), j + 1, k ) ) / dx, 0.5f * ( value( i, index( dim, j + 2 ), k ) - s[0] ) / dy,
                     0.5f * ( s[6] - value( i, j + 1, index( dim, k - 1 ) ) ) / dz };

            // vertex 3
            n[3] = { 0.5f * ( value( index( dim, i + 2 ), j + 1, k ) - s[2] ) / dx, 0.5f * ( value( i + 1, index( dim, j + 2 ), k ) - s[1] ) / dy,
                     0.5f * ( s[7] - value( i + 1, j + 1, index( dim, k - 1 ) ) ) / dz };

            // vertex 4
            n[4] = { 0.5f * ( s[5] - value( index( dim, i - 1 ), j, k + 1 ) ) / dx, 0.5f * ( s[6] - value( i, index( dim, j - 1 ), k + 1 ) ) / dy,
                     0.5f * ( value( i, j, index( dim, k + 2 ) ) - s[0] ) / dz };

            // vertex 5
            n[5] = { 0.5f * ( value( index( dim, i + 2 ), j, k + 1 ) - s[4] ) / dx, 0.5f * ( s[7] - value( i + 1, index( dim, j - 1 ), k + 1 ) ) / dy,
                     0.5f * ( value( i + 1, j, index( dim, k + 2 ) ) - s[1] ) / dz };

            // vertex 6
            n[6] = { 0.5f * ( s[7] - value( index( dim, i - 1 ), j + 1, k + 1 ) ) / dx, 0.5f * ( value( i, index( dim, j + 2 ), k + 1 ) - s[4] ) / dy,
                     0.5f * ( value( i, j + 1, index( dim, k + 2 ) ) - s[2] ) / dz };

            // vertex 7
            n[7] = { 0.5f * ( value( index( dim, i + 2 ), j + 1, k + 1 ) - s[6] ) / dx,
                     0.5f * ( value( i + 1, index( dim, j + 2 ), k + 1 ) - s[5] ) / dy,
                     0.5f * ( value( i + 1, j + 1, index( dim, k + 2 ) ) - s[3] ) / dz };
        }

        Point_3 position( const std::size_t x, const std::size_t y, const std::size_t z ) const {
            const FT vx = octree_->hx();
            const FT vy = octree_->hx();
            const FT vz = octree_->hx();

            return Point_3( x * vx + octree_->offset_x(), y * vy + octree_->offset_y(), z * vz + octree_->offset_z() );
        }

        FT value( const std::size_t x, const std::size_t y, const std::size_t z ) const { return octree_->value( x, y, z ); }

        std::array<FT, 8> voxel_values( const std::size_t i, const std::size_t j, const std::size_t k ) const {
            return octree_->voxel_values( i, j, k );
        }

        std::array<Point_3, 8> voxel_vertex_positions( const std::size_t i, const std::size_t j, const std::size_t k ) const {
            return octree_->voxel_vertex_positions( i, j, k );
        }

        bool exists( const std::size_t& i, const std::size_t& j, const std::size_t& k ) const { return octree_->exists( i, j, k ); }

        std::size_t n_edges() const { return octree_->leaf_edges().size(); }

        std::array<FT, 2> edge_values( const Octree_edge_index& e_id ) const { return octree_->edge_values( e_id ); }
        std::array<std::size_t, 4> voxels_incident_to_edge( const Octree_edge_index& e_id ) const { return octree_->voxels_incident_to_edge( e_id ); }

        const std::vector<Octree_edge_index>& edges() const { return octree_edges_; }

      private:
        const OctreeWrapper* octree_;
        std::vector<Octree_edge_index> octree_edges_;
    };

}    // namespace CGAL

#endif    // CGAL_OCTREE_GRID_ORACLE_H
