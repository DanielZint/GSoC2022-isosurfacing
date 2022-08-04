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

        typedef typename OctreeWrapper::Vertex_handle Vertex_handle;
        typedef typename OctreeWrapper::Edge_handle Edge_handle;
        typedef typename OctreeWrapper::Voxel_handle Voxel_handle;

      public:
        Octree_oracle( const OctreeWrapper& octree )
            : octree_( &octree ), octree_edges_( octree.leaf_edges().begin(), octree.leaf_edges().end() ),
              octree_vertices_( octree.leaf_vertices().begin(), octree.leaf_vertices().end() ),
              octree_voxels_( octree.leaf_voxels().begin(), octree.leaf_voxels().end() ) {}

        /// compute the gradient of the scalar field with central difference
        void gradient( std::array<Vector_3, 8>& n, const std::array<float, 8>& s, const std::size_t& i, const std::size_t& j,
                       const std::size_t& k ) const {
            const auto node       = octree_->get_node( i, j, k );
            const auto df         = octree_->depth_factor( node );
            const std::size_t dim = octree_->dim() + 1;
            auto index            = [&df, &dim]( const int ii ) { return ( ii < 0 ? 0 : ii >= dim ? ( dim - 1 ) : ii ); };

            const FT& dx = octree_->hx() * df;
            const FT& dy = octree_->hx() * df;
            const FT& dz = octree_->hx() * df;

            //n[0] = { ( s[1] - s[0] ) / dx, ( s[2] - s[0] ) / dy, ( s[4] - s[0] ) / dz };
            //n[1] = { ( s[1] - s[0] ) / dx, ( s[3] - s[1] ) / dy, ( s[5] - s[1] ) / dz };
            //n[2] = { ( s[3] - s[2] ) / dx, ( s[2] - s[0] ) / dy, ( s[6] - s[2] ) / dz };
            //n[3] = { ( s[3] - s[2] ) / dx, ( s[3] - s[1] ) / dy, ( s[7] - s[3] ) / dz };
            //n[4] = { ( s[5] - s[4] ) / dx, ( s[6] - s[4] ) / dy, ( s[4] - s[0] ) / dz };
            //n[6] = { ( s[7] - s[6] ) / dx, ( s[6] - s[4] ) / dy, ( s[6] - s[2] ) / dz };
            //n[7] = { ( s[7] - s[6] ) / dx, ( s[7] - s[5] ) / dy, ( s[7] - s[3] ) / dz };

            // vertex 0
            n[0] = { 0.5f * ( s[1] - value( index( i - df ), j, k ) ) / dx, 0.5f * ( s[2] - value( i, index( j - df ), k ) ) / dy,
                     0.5f * ( s[4] - value( i, j, index( k - df ) ) ) / dz };
            // vertex 1
            n[1] = { 0.5f * ( value( index( i + 2 * df ), j, k ) - s[0] ) / dx, 0.5f * ( s[3] - value( i + df, index( j - df ), k ) ) / dy,
                     0.5f * ( s[5] - value( i + df, j, index( k - df ) ) ) / dz };
            // vertex 2
            n[2] = { 0.5f * ( s[3] - value( index( i - df ), j + df, k ) ) / dx, 0.5f * ( value( i, index( j + 2 * df ), k ) - s[0] ) / dy,
                     0.5f * ( s[6] - value( i, j + df, index( k - df ) ) ) / dz };
            // vertex 3
            n[3] = { 0.5f * ( value( index( i + 2 * df ), j + df, k ) - s[2] ) / dx, 0.5f * ( value( i + df, index( j + 2 * df ), k ) - s[1] ) / dy,
                     0.5f * ( s[7] - value( i + df, j + df, index( k - df ) ) ) / dz };
            // vertex 4
            n[4] = { 0.5f * ( s[5] - value( index( i - df ), j, k + df ) ) / dx, 0.5f * ( s[6] - value( i, index( j - df ), k + df ) ) / dy,
                     0.5f * ( value( i, j, index( k + 2 * df ) ) - s[0] ) / dz };
            // vertex 5
            n[5] = { 0.5f * ( value( index( i + 2 * df ), j, k + df ) - s[4] ) / dx, 0.5f * ( s[7] - value( i + df, index( j - df ), k + df ) ) / dy,
                     0.5f * ( value( i + df, j, index( k + 2 * df ) ) - s[1] ) / dz };
            // vertex 6
            n[6] = { 0.5f * ( s[7] - value( index( i - df ), j + df, k + df ) ) / dx, 0.5f * ( value( i, index( j + 2 * df ), k + df ) - s[4] ) / dy,
                     0.5f * ( value( i, j + df, index( k + 2 * df ) ) - s[2] ) / dz };
            // vertex 7
            n[7] = { 0.5f * ( value( index( i + 2 * df ), j + df, k + df ) - s[6] ) / dx,
                     0.5f * ( value( i + df, index( j + 2 * df ), k + df ) - s[5] ) / dy,
                     0.5f * ( value( i + df, j + df, index( k + 2 * df ) ) - s[3] ) / dz };
        }

        void gradient( std::array<Vector_3, 8>& n, const std::array<float, 8>& s, const Voxel_handle& vh ) const {
            n = octree_->voxel_gradients( vh );
        }

        Point_3 position( const std::size_t x, const std::size_t y, const std::size_t z ) const {
            const FT vx = octree_->hx();
            const FT vy = octree_->hx();
            const FT vz = octree_->hx();

            return Point_3( x * vx + octree_->offset_x(), y * vy + octree_->offset_y(), z * vz + octree_->offset_z() );
        }

        Point_3 position( const Vertex_handle& v ) const { return octree_->point( v ); }

        FT value( const std::size_t x, const std::size_t y, const std::size_t z ) const { return octree_->value( x, y, z ); }

        std::array<FT, 8> voxel_values( const Vertex_handle& vh ) const { return octree_->voxel_values( vh ); }

        std::array<Point_3, 8> voxel_vertex_positions( const Voxel_handle& vh ) const { return octree_->voxel_vertex_positions( vh ); }

        std::size_t n_edges() const { return octree_edges_.size(); }

        const std::vector<Octree_edge_index>& edges() const { return octree_edges_; }

        std::size_t n_vertices() const { return octree_vertices_.size(); }

        const std::vector<OctreeWrapper::Vertex_handle>& vertices() const { return octree_vertices_; }

        std::size_t n_voxels() const { return octree_voxels_.size(); }

        const std::vector<OctreeWrapper::Vertex_handle>& voxels() const { return octree_voxels_; }

        std::array<FT, 2> edge_values( const Octree_edge_index& e_id ) const { return octree_->edge_values( e_id ); }
        std::array<std::size_t, 4> voxels_incident_to_edge( const Octree_edge_index& e_id ) const { return octree_->voxels_incident_to_edge( e_id ); }

      private:
        const OctreeWrapper* octree_;
        std::vector<Octree_edge_index> octree_edges_;
        std::vector<Vertex_handle> octree_vertices_;
        std::vector<Voxel_handle> octree_voxels_;
    };

}    // namespace CGAL

#endif    // CGAL_OCTREE_GRID_ORACLE_H
