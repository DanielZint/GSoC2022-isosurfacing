#ifndef CGAL_DUAL_CONTOURING_3_H
#define CGAL_DUAL_CONTOURING_3_H

#include "Isosurfacing_3/internal/Marching_cubes_3_internal.h"
#include "util.h"

#include <Eigen/SVD>
#include <array>
#include <map>
#include <vector>

namespace CGAL {

    /// <summary>
    /// Compute vertex position for Dual Contouring
    /// </summary>
    /// <typeparam name="Domain_"></typeparam>
    /// <param name="domain"></param>
    /// <param name="iso_value"></param>
    /// <param name="i"></param>
    /// <param name="j"></param>
    /// <param name="k"></param>
    /// <returns> true, if there is a point in the cell</returns>
    template<class Domain_, bool use_bbox = false>
    bool get_vertex_position( const Domain_& domain, const typename Domain_::FT iso_value, const typename Domain_::Voxel_handle& vh,
                              typename Domain_::Point_3& point ) {
        typedef typename Domain_::Point_3 Point_3;
        typedef typename Domain_::Vector_3 Vector_3;

        std::array<Domain_::FT, Tables::N_VERTICES> s = domain.voxel_values( vh );

        std::array<bool, Tables::N_VERTICES> b;
        std::transform( s.begin(), s.end(), b.begin(), [iso_value]( const auto& e ) { return e <= iso_value; } );

        unsigned int cubeindex = 0;
        // set bit if corresponding corner is below iso
        for( int i = 0; i < Tables::N_VERTICES; ++i ) {
            cubeindex |= b[i] << i;
        }

        if( cubeindex == 0 || cubeindex == 255 ) {
            return false;
        }

        std::array<Point_3, Tables::N_VERTICES> p = domain.voxel_vertex_positions( vh );
        std::array<Vector_3, Tables::N_VERTICES> pos;
        std::transform( p.begin(), p.end(), pos.begin(), []( const auto& e ) { return e - CGAL::ORIGIN; } );

        point = CGAL::ORIGIN + ( pos[0] + 0.5 * ( pos[7] - pos[0] ) );    // set point to voxel center

        std::array<Vector_3, Tables::N_VERTICES> normals = domain.gradient( vh );

        // compute edge intersections
        std::vector<Point_3> edge_intersections;
        std::vector<Vector_3> edge_intersection_normals;

        for( int i = 0; i < Tables::N_EDGES; ++i ) {
            const auto& v0 = Tables::edge_to_vertex[i][0];
            const auto& v1 = Tables::edge_to_vertex[i][1];

            if( b[v0] != b[v1] ) {    // e0
                const FT u           = ( s[v0] - iso_value ) / ( s[v0] - s[v1] );
                const Point_3 p_lerp = CGAL::ORIGIN + ( ( 1 - u ) * pos[v0] + u * pos[v1] );
                edge_intersections.push_back( p_lerp );
                const Vector_3 n_lerp = ( 1 - u ) * normals[v0] + u * normals[v1];
                edge_intersection_normals.push_back( n_lerp );
            }
        }

        // MC Polygon Center of Mass
        if( false ) {
            Vector_3 com_vec( 0, 0, 0 );

            for( int i = 0; i < edge_intersections.size(); ++i ) {
                com_vec += edge_intersections[i] - CGAL::ORIGIN;
            }

            Point_3 p = CGAL::ORIGIN + com_vec / edge_intersections.size();
            point     = p;
        }

        // SVD QEM
        if( true ) {
            Eigen::Matrix3d A;
            A.setZero();
            Eigen::Vector3d b;
            b.setZero();
            for( int i = 0; i < edge_intersections.size(); ++i ) {
                Eigen::Vector3d n_k = { edge_intersection_normals[i].x(), edge_intersection_normals[i].y(), edge_intersection_normals[i].z() };
                Eigen::Vector3d p_k = { edge_intersections[i].x(), edge_intersections[i].y(), edge_intersections[i].z() };
                double d_k          = n_k.transpose() * p_k;

                Eigen::Matrix3d A_k = n_k * n_k.transpose();
                Eigen::Vector3d b_k = d_k * n_k;
                A += A_k;
                b += b_k;
            }

            Eigen::JacobiSVD<Eigen::MatrixXd> svd( A, Eigen::ComputeThinU | Eigen::ComputeThinV );
            // set threshold as in Peter Lindstrom's paper, "Out-of-Core
            // Simplification of Large Polygonal Models"
            svd.setThreshold( 1e-3 );

            // Init x hat
            Eigen::Vector3d x_hat;
            x_hat << point.x(), point.y(), point.z();

            // Lindstrom formula for QEM new position for singular matrices
            Eigen::VectorXd v_svd = x_hat + svd.solve( b - A * x_hat );

            point = Point_3( v_svd[0], v_svd[1], v_svd[2] );
        }

        // bbox
        if constexpr( use_bbox ) {
            CGAL::Bbox_3 bbox = ( CGAL::ORIGIN + pos[0] ).bbox() + ( CGAL::ORIGIN + pos[7] ).bbox();

            FT x  = std::min<FT>( std::max<FT>( point.x(), bbox.xmin() ), bbox.xmax() );
            FT y  = std::min<FT>( std::max<FT>( point.y(), bbox.ymin() ), bbox.ymax() );
            FT z  = std::min<FT>( std::max<FT>( point.z(), bbox.zmin() ), bbox.zmax() );
            point = Point_3( x, y, z );
        }

        return true;
    }

    template<class Domain_, class PointRange, class PolygonRange>
    void make_quad_mesh_using_dual_contouring( const Domain_& domain, const typename Domain_::FT iso_value, PointRange& points,
                                               PolygonRange& polygons, const bool use_bbox = false ) {
        typedef typename Domain_::FT FT;
        typedef typename Domain_::Point_3 Point_3;

        ScopeTimer timer;

        // compute dc-vertices
        std::map<size_t, size_t> map_voxel_to_point_id;
        std::map<size_t, Point_3> map_voxel_to_point;
        size_t points_counter = 0;
        std::map<Octree_edge_index, std::array<size_t, 4>> quads;

        const std::size_t n_voxel = domain.n_voxels();

        // save all points
        if( use_bbox ) {
            for( int i = 0; i < n_voxel; ++i ) {
                const auto& vh = domain.voxels(i);
                Point_3 p;
                if( get_vertex_position<Domain_, true>( domain, iso_value, vh, p ) ) {
                    map_voxel_to_point[vh]    = p;
                    map_voxel_to_point_id[vh] = points_counter++;
                }
            }
        } else {
            for( int i = 0; i < n_voxel; ++i ) {
                const auto& vh = domain.voxels(i);
                Point_3 p;
                if( get_vertex_position( domain, iso_value, vh, p ) ) {
                    map_voxel_to_point[vh]    = p;
                    map_voxel_to_point_id[vh] = points_counter++;
                }
            }
        }

        const std::size_t n_edges = domain.n_edges();
        
        // save all quads
        for( int i = 0; i < n_edges; ++i ) {
            const auto& e        = domain.edges( i );
            const auto& [s0, s1] = domain.edge_values( e );

            if( s0 <= iso_value && s1 > iso_value ) {
                const auto voxels = domain.voxels_incident_to_edge( e );
                quads[e]          = voxels;
            } else if( s1 <= iso_value && s0 > iso_value ) {
                const auto voxels = domain.voxels_incident_to_edge( e );
                quads[e]          = { voxels[0], voxels[3], voxels[2], voxels[1] };
            }
        }

        // write points and quads in ranges
        points.resize( points_counter );
        for( const auto& vtop: map_voxel_to_point ) {
            points[map_voxel_to_point_id[vtop.first]] = vtop.second;
        }

        polygons.reserve( quads.size() );
        for( const auto& q: quads ) {
            std::vector<std::size_t> vertex_ids;
            for( const auto& v_id: q.second ) {
                vertex_ids.push_back( map_voxel_to_point_id[v_id] );
            }
            polygons.push_back( vertex_ids );
        }
    }

}    // namespace CGAL

#endif    // CGAL_DUAL_CONTOURING_3_H
