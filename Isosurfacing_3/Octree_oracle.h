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

        std::array<Vector_3, 8> gradient( const Voxel_handle& vh ) const { return octree_->voxel_gradients( vh ); }

        Point_3 position( const Vertex_handle& v ) const { return octree_->point( v ); }

        std::array<FT, 8> voxel_values( const Vertex_handle& vh ) const { return octree_->voxel_values( vh ); }

        std::array<Point_3, 8> voxel_vertex_positions( const Voxel_handle& vh ) const { return octree_->voxel_vertex_positions( vh ); }

        std::size_t n_edges() const { return octree_edges_.size(); }
        std::size_t n_vertices() const { return octree_vertices_.size(); }
        std::size_t n_voxels() const { return octree_voxels_.size(); }

        const Octree_edge_index& edges( const std::size_t& i ) const { return octree_edges_[i]; }
        const OctreeWrapper::Vertex_handle& vertices(const std::size_t& i) const { return octree_vertices_[i]; }
        const OctreeWrapper::Vertex_handle& voxels( const std::size_t& i ) const { return octree_voxels_[i]; }

        std::array<FT, 2> edge_values( const Octree_edge_index& e_id ) const { return octree_->edge_values( e_id ); }
        std::array<std::size_t, 4> voxels_incident_to_edge( const Octree_edge_index& e_id ) const { return octree_->edge_voxels( e_id ); }

      private:
        const OctreeWrapper* octree_;
        std::vector<Octree_edge_index> octree_edges_;
        std::vector<Vertex_handle> octree_vertices_;
        std::vector<Voxel_handle> octree_voxels_;
    };

}    // namespace CGAL

#endif    // CGAL_OCTREE_GRID_ORACLE_H
