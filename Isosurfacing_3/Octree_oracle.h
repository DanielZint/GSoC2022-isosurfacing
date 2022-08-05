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

        typedef typename Octree_wrapper::Vertex_handle Vertex_handle;
        typedef typename Octree_wrapper::Edge_handle Edge_handle;
        typedef typename Octree_wrapper::Voxel_handle Voxel_handle;

      public:
        Octree_oracle( const Octree_wrapper& octree ) : octree_( &octree ) {}

        std::array<Vector_3, 8> gradient( const Voxel_handle& vh ) const { return octree_->voxel_gradients( vh ); }

        Point_3 position( const Vertex_handle& v ) const { return octree_->point( v ); }

        std::array<FT, 8> voxel_values( const Vertex_handle& vh ) const { return octree_->voxel_values( vh ); }

        std::array<Point_3, 8> voxel_vertex_positions( const Voxel_handle& vh ) const { return octree_->voxel_vertex_positions( vh ); }

        std::size_t n_edges() const { return octree_->leaf_edges().size(); }
        std::size_t n_vertices() const { return octree_->leaf_vertices().size(); }
        std::size_t n_voxels() const { return octree_->leaf_voxels().size(); }

        const Edge_handle& edges( const std::size_t& i ) const { return octree_->leaf_edges()[i]; }
        const Octree_wrapper::Vertex_handle& vertices( const std::size_t& i ) const { return octree_->leaf_vertices()[i]; }
        const Octree_wrapper::Vertex_handle& voxels( const std::size_t& i ) const { return octree_->leaf_voxels()[i]; }

        std::array<FT, 2> edge_values( const Edge_handle& e_id ) const { return octree_->edge_values( e_id ); }
        std::array<std::size_t, 4> voxels_incident_to_edge( const Edge_handle& e_id ) const { return octree_->edge_voxels( e_id ); }

      private:
        const Octree_wrapper* octree_;
    };

}    // namespace CGAL

#endif    // CGAL_OCTREE_GRID_ORACLE_H
