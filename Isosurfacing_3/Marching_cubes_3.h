#ifndef CGAL_MARCHING_CUBES_3_H
#define CGAL_MARCHING_CUBES_3_H

#include "Isosurfacing_3/internal/Marching_cubes_3_internal.h"
#include "util.h"

#include <array>
#include <mutex>

namespace CGAL {

    template<class Domain_, class PointRange, class PolygonRange>
    void make_triangle_mesh_using_marching_cubes( const Domain_& domain, const typename Domain_::FT iso_value, PointRange& points,
                                                  PolygonRange& polygons ) {
        ScopeTimer timer;

        std::mutex mutex;

        const std::size_t size_k = domain.size_x();
        const std::size_t size_j = domain.size_y();
        const std::size_t size_i = domain.size_z();

        const std::size_t blocking_size = 100;

        // TODO: look at polygon mesh processing for tbb (also linking)

#pragma omp parallel for
        for( int bj = 1; bj < size_j - 1; bj += blocking_size ) {
            for( int k = 1; k < size_k - 1; k++ ) {
                const std::size_t j_start = bj;
                const std::size_t j_end   = std::min( size_j - 1, bj + blocking_size );

                for( std::size_t j = j_start; j < j_end; j++ ) {
                    for( std::size_t i = 1; i < size_i - 1; i++ ) {
                        internal::Marching_cubes_3::marching_cubes_cell( i, j, k, domain, iso_value, points, polygons, mutex );
                    }
                }
            }
        }
    }

}    // namespace CGAL

#endif    // CGAL_MARCHING_CUBES_3_H
