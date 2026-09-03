#include "utils_H.h"

//! std library
#include <cstdlib>

#include <cmath>
#include <iosfwd>
#include <limits>
#include <numeric>
#include <string>

#include <chrono>
#include <ctime>
#include <functional>

//==============================================================================

double compute_dt_sediment(const double alpha, const double beta,
                           const double S_x, const double S_y,
                           const std::vector<double> &u,
                           const std::vector<double> &v,
                           const double pixel_size, const double dt_DSV,
                           unsigned int *numberOfSteps) {
  static constexpr double max_courant_number = .95;
  // +-----------------------------------------------+
  // |      Estimate vertical max Courant number     |
  // +-----------------------------------------------+

  const double dt_y =
      max_courant_number * pixel_size /
      (alpha * std::pow(S_y, beta) *
       std::max(*std::max_element(v.begin(), v.end()),
                std::abs(*std::min_element(v.begin(), v.end()))));

  // +-----------------------------------------------+
  // |    Estimate horizontal max Courant number     |
  // +-----------------------------------------------+

  const double dt_x =
      max_courant_number * pixel_size /
      (alpha * std::pow(S_x, beta) *
       std::max(*std::max_element(u.begin(), u.end()),
                std::abs(*std::min_element(u.begin(), u.end()))));

  double dt_sed = std::min(dt_y, dt_x);

  dt_sed = std::min(dt_DSV / std::floor(dt_DSV / dt_sed), dt_DSV);

  *numberOfSteps = std::floor(dt_DSV / dt_sed);

  return dt_sed;
}

//==============================================================================

int current_start_chunk(const int &rank,
                        const std::vector<int> &chunk_length_vec) {
  int result;
  if (rank - 1 < 0) {
    result = 0;
  } else {
    result = chunk_length_vec[rank - 1] +
             current_start_chunk(rank - 1, chunk_length_vec);
  }
  return result;
}

//==============================================================================

#ifdef ENABLE_CUDA
void make_sparsity_pattern(thrust::host_vector<unsigned int>& idBasinVect,
		std::vector<unsigned int>& basin_mask,
		thrust::host_vector<unsigned int>& idBasinVectReIndex,
#else
void make_sparsity_pattern(std::vector<unsigned int>& idBasinVect,
		std::vector<unsigned int>& basin_mask,
		std::vector<unsigned int>& idBasinVectReIndex,
#endif
		int *row_offsets, int* columns,
		unsigned int const N_rows, unsigned int const N_cols) {

    int it = 0; // next unused index into `columns`/`values`

    // Solve-cell bitmap taken straight from idBasinVect, so the pattern is always
    // consistent with idBasinVectReIndex and with the number of unknowns
    // m = idBasinVect.size().  The `basin_mask` argument is the *non-excluded*
    // mask and must NOT be used here: it would couple solve rows to excluded
    // cells, write more entries than the CSR was sized for and corrupt it
    // (this is what crashed the 10 m case in the IC factorization).
    (void)basin_mask;
    std::vector<unsigned char> is_solve((size_t)N_rows * N_cols, 0);
    for (const auto Id : idBasinVect)
        is_solve[Id] = 1;

#define INSERT(u,v)                                             \
    if(0<=(u) && (u)<(int)N_rows &&                             \
       0<=(v) && (v)<(int)N_cols &&                             \
	is_solve[((u) * N_cols + (v))]==1)                      \
    {                                                           \
        columns[it] = idBasinVectReIndex[((u) * N_cols + (v))]; \
        ++it;                                                   \
    }

    int row = 0;
    row_offsets[row] = 0;
    for (const auto Id : idBasinVect) {
	const int i = (int)(Id / N_cols);
	const int j = (int)(Id % N_cols);
 
	// Our sparsity pattern is a star-stencil
	INSERT(i - 1, j    );
        INSERT(i    , j - 1);
        INSERT(i    , j    );
        INSERT(i    , j + 1);
        INSERT(i + 1, j    );
        row_offsets[++row] = it;
    }
#undef INSERT
}

//==============================================================================
