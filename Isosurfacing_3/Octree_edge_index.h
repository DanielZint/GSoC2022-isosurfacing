#pragma once

#include "types.h"

#include <CGAL/Octree.h>
#include <iostream>

class Octree_edge_index {
  public:
    typedef CGAL::Octree<Kernel, std::vector<Point_3>> Octree;
    typedef Octree::Node::Global_coordinates Octree_vertex_id;
    typedef std::tuple<Octree_vertex_id, std::size_t> type;

  private:
    type id_;

  public:
    Octree_edge_index( const Octree_vertex_id& v_id, const std::size_t& e_local_id ) : id_ { v_id, e_local_id } {}

    operator type() { return id_; }

    Octree_vertex_id v_id() const { return std::get<0>( id_ ); }
    std::size_t e_local_id() const { return std::get<1>( id_ ); }
};

std::ostream& operator<<( std::ostream& os, const Octree_edge_index& oei ) {
    auto v_id = oei.v_id();

    for( int i = 0; i < Octree_edge_index::Octree::Dimension::value; ++i ) {
        os << v_id[i] << ' ';
    }

    os << "e_local: " << oei.e_local_id();
    return os;
}