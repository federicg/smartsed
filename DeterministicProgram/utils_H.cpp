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

void compute_dt_adaptive(const std::vector<double> &H,
                         const std::vector<double> &H_old,
                         const std::vector<double> &H_oldold,
                         const std::vector<unsigned int> &idBasinVect,
                         double &dt,
                         const double &local_estimator_time_tolerance,
                         const double &time, const double &timed,
                         const double &timedd) {

  double dh_t, h1, h2, h3, a_coeff, b_coeff, Nu_hmean_cell = 0., nu_htot = 0.;
  for (const auto &Id : idBasinVect) {
    const double &hcell = H[Id], &hcell_old = H_old[Id],
                 &hcell_oldd = H_oldold[Id];

    dh_t = (hcell - hcell_old) / (time - timed);

    h1 = hcell_oldd / ((timedd - timed) * (timedd - time));
    h2 = hcell_old / ((timed - timedd) * (timed - time));
    h3 = hcell / ((time - timedd) * (time - timed));

    a_coeff = h1 + h2 + h3;
    b_coeff =
        -(h1 * (time + timed) + h2 * (time + timedd) + h3 * (timed + timedd));

    Nu_hmean_cell += (1. / 3. * a_coeff * a_coeff *
                          (time * time + time * timed + timed * timed) +
                      a_coeff * (b_coeff - dh_t) * (time + timed) +
                      (b_coeff - dh_t) * (b_coeff - dh_t));
    nu_htot += Nu_hmean_cell * (time - timed) * (time - timed);
  }

  double dt_candidate =
      local_estimator_time_tolerance / std::sqrt(nu_htot) * (time - timed);

  dt_candidate = (nu_htot > 0 && dt_candidate < dt) ? dt_candidate : dt;

  dt = dt_candidate;
}

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
