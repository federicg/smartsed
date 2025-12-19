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

void bilinearInterpolation(
    const std::vector<double> &u, const std::vector<double> &v,
    std::vector<double> &u_star, std::vector<double> &v_star,
    const unsigned int &nrows, const unsigned int &ncols, const double &dt_DSV,
    const double &pixel_size,
    const std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    const std::vector<unsigned int> &idStaggeredInternalVectVertical,
    const std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    const std::vector<unsigned int> &idStaggeredBoundaryVectEast,
    const std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    const std::vector<unsigned int> &idStaggeredBoundaryVectSouth) {

  // +-----------------------------------------------+
  // |              Horizontal Velocity              |
  // +-----------------------------------------------+

  for (const auto &Id : idStaggeredInternalVectHorizontal) {

    const unsigned int i = Id / (ncols + 1), j = Id % (ncols + 1),
                       ID_NE = Id - i, // v
        ID_NW = ID_NE - 1,             // v
        ID_SE = ID_NE + ncols,         // v
        ID_SW = ID_NW + ncols;         // v

    Vector2D vel(std::array<double, 2>{
        {u[Id], (v[ID_NE] + v[ID_NW] + v[ID_SE] + v[ID_SW]) / 4.}});

    const auto Dx = vel / pixel_size * dt_DSV;

    Vector2D xx(std::array<double, 2>{
        {double(j), double(i)}}); // Top-Left reference frame

    xx = xx - Dx;

    auto x = xx(0), y = xx(1);

    auto x_1 = std::floor(x), //  11    21     ----> x
        y_1 = std::floor(y),  //
        x_2 = x_1 + 1,        //
        y_2 = y_1 + 1;        //  12    22
                              //
                              // |
                              // | y

    if (x_1 == ncols) {
      x_2 -= 1;
      x_1 -= 1;
      x -= 1;
    }
    if (x_2 == 0) {
      x_2 += 1;
      x_1 += 1;
      x += 1;
    }

    if (y_1 == (nrows - 1)) {
      y_2 -= 1;
      y_1 -= 1;
      y -= 1;
    }
    if (y_2 == 0) {
      y_2 += 1;
      y_1 += 1;
      y += 1;
    }

    // return 0 for target values that are out of bounds
    if (x_2 < 0 || x_1 > ncols || y_2 < 0 || y_1 > (nrows - 1)) {
      u_star[Id] = 0.;
    } else {
      const auto Id_11 = x_1 + y_1 * (ncols + 1), // u
          Id_12 = x_1 + y_2 * (ncols + 1),        // u
          Id_21 = x_2 + y_1 * (ncols + 1),        // u
          Id_22 = x_2 + y_2 * (ncols + 1);        // u

      // compute weights
      const double w_x2 = x_2 - x, w_x1 = x - x_1, w_y2 = y_2 - y,
                   w_y1 = y - y_1;

      const auto a = u[Id_11] * w_x2 + u[Id_21] * w_x1,
                 b = u[Id_12] * w_x2 + u[Id_22] * w_x1;

      u_star[Id] = a * w_y2 + b * w_y1;
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectWest) {
    const auto Idd = Id + 1;
    u_star[Id] = u_star[Idd] * (u[Idd] < 0.);
  }

  for (const auto &Id : idStaggeredBoundaryVectEast) {
    const auto Idd = Id - 1;
    u_star[Id] = u_star[Idd] * (u[Idd] > 0.);
  }

  // +-----------------------------------------------+
  // |              Vertical Velocity                |
  // +-----------------------------------------------+

  for (const auto &Id : idStaggeredInternalVectVertical) {

    const auto i = Id / ncols, j = Id % ncols,

               ID_SW = Id + i,       // u
        ID_SE = ID_SW + 1,           // u
        ID_NW = ID_SW - (ncols + 1), // u
        ID_NE = ID_NW + 1;           // u

    Vector2D vel(std::array<double, 2>{
        {(u[ID_SW] + u[ID_SE] + u[ID_NW] + u[ID_NE]) / 4., v[Id]}});

    const auto Dx = vel / pixel_size * dt_DSV;

    Vector2D xx(std::array<double, 2>{
        {double(j), double(i)}}); // Top-Left reference frame

    xx = xx - Dx;

    auto x = xx(0), y = xx(1);

    auto x_1 = std::floor(x), //  11    21     ----> x
        y_1 = std::floor(y),  //
        x_2 = x_1 + 1,        //
        y_2 = y_1 + 1;        //  12    22
                              //
                              // |
                              // | y

    if (x_1 == (ncols - 1)) {
      x_2 -= 1;
      x_1 -= 1;
      x -= 1;
    }
    if (x_2 == 0) {
      x_2 += 1;
      x_1 += 1;
      x += 1;
    }

    if (y_1 == nrows) {
      y_2 -= 1;
      y_1 -= 1;
      y -= 1;
    }
    if (y_2 == 0) {
      y_2 += 1;
      x_2 += 1;
      y += 1;
    }

    // return 0 for target values that are out of bounds
    if (x_2 < 0 || x_1 > (ncols - 1) || y_2 < 0 || y_1 > nrows) {
      v_star[Id] = 0.;
    } else {
      const auto Id_11 = x_1 + y_1 * ncols, Id_12 = x_1 + y_2 * ncols,
                 Id_21 = x_2 + y_1 * ncols, Id_22 = x_2 + y_2 * ncols;

      // compute weights
      const double w_x2 = x_2 - x, w_x1 = x - x_1, w_y2 = y_2 - y,
                   w_y1 = y - y_1;

      const auto a = v[Id_11] * w_x2 + v[Id_21] * w_x1,
                 b = v[Id_12] * w_x2 + v[Id_22] * w_x1;

      v_star[Id] = a * w_y2 + b * w_y1;
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectNorth) {
    const auto Idd = Id + ncols;
    v_star[Id] = v_star[Idd] * (v[Idd] < 0.);
  }

  for (const auto &Id : idStaggeredBoundaryVectSouth) {
    const auto Idd = Id - ncols;
    v_star[Id] = v_star[Idd] * (v[Idd] > 0.);
  }
}

//==============================================================================

void buildMatrix(
    const std::vector<double> &H_int_x, const std::vector<double> &H_int_y,
    const std::vector<double> &orography, const std::vector<double> &u_star,
    const std::vector<double> &v_star, const std::vector<double> &u,
    const std::vector<double> &v, const std::vector<double> &H,
    const unsigned int &N_cols, const unsigned int &N_rows,
    const unsigned int &N, const double &c1, const double &c3,
    const double &H_min, const std::vector<double> &precipitation,
    const double &dt_DSV, const std::vector<double> &alfa_x,
    const std::vector<double> &alfa_y,
    const std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    const std::vector<unsigned int> &idStaggeredInternalVectVertical,
    const std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    const std::vector<unsigned int> &idStaggeredBoundaryVectEast,
    const std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    const std::vector<unsigned int> &idStaggeredBoundaryVectSouth,
    const std::vector<unsigned int> &idBasinVect,
    const std::vector<unsigned int> &idBasinVect_not_excluded,
    const std::vector<unsigned int>
        &idStaggeredInternalVectHorizontal_not_excluded,
    const std::vector<unsigned int>
        &idStaggeredInternalVectVertical_not_excluded,
    const std::vector<unsigned int> &idBasinVectReIndex,
    const bool &isNonReflectingBC, const bool &isH,

    const std::vector<std::tuple<bool, int>>
        &excluded_ids, // bool, value 1 in excluded ids  unsigned int, Id del
                       // pour point
    std::vector<double> &additional_source_term,

    std::vector<Eigen::Triplet<double>> &coefficients, Eigen::VectorXd &rhs) {

  // Be careful to the mass conservation
  // cycle over boundary interfaces and check if one of the left, right cells
  // are excluded ones

  for (const auto &k_ex : idBasinVect_not_excluded) {
    const auto &current_tuple = excluded_ids[k_ex];
    if (std::get<0>(current_tuple)) {
      const auto &k_pour = std::get<1>(current_tuple);
      if (k_pour >= 0) {
        additional_source_term[k_pour] += precipitation[k_ex] * dt_DSV;
      }
    }
  }

  for (const auto &Id : idBasinVect) {
    const auto IDreIndex = idBasinVectReIndex[Id];

    coefficients.push_back(Eigen::Triplet<double>(IDreIndex, IDreIndex, 1.));
    rhs(IDreIndex) = H[Id] + precipitation[Id] * dt_DSV;
  }

  for (const auto &Id : idStaggeredInternalVectHorizontal) {
    const unsigned int i = Id / (N_cols + 1),
                       IDleft = Id - i - 1,           // H
        IDright = Id - i,                             // H
        IDleftReIndex = idBasinVectReIndex[IDleft],   // H
        IDrightReIndex = idBasinVectReIndex[IDright]; // H

    // define H at interfaces
    const auto H_interface = H_int_x[Id];
    const double coeff_m = H_interface * alfa_x[Id];

    if (H_interface > H_min) {
      coefficients.push_back(
          Eigen::Triplet<double>(IDleftReIndex, IDrightReIndex, -c3 * coeff_m));
      coefficients.push_back(
          Eigen::Triplet<double>(IDleftReIndex, IDleftReIndex, c3 * coeff_m));

      rhs(IDleftReIndex) +=
          -c1 * (+coeff_m * u_star[Id]) -
          (orography[IDleft] - orography[IDright]) * c3 * coeff_m * isH;

      coefficients.push_back(
          Eigen::Triplet<double>(IDrightReIndex, IDleftReIndex, -c3 * coeff_m));
      coefficients.push_back(
          Eigen::Triplet<double>(IDrightReIndex, IDrightReIndex, c3 * coeff_m));

      rhs(IDrightReIndex) +=
          -c1 * (-coeff_m * u_star[Id]) -
          (orography[IDright] - orography[IDleft]) * c3 * coeff_m * isH;
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectWest) {
    const unsigned int i = Id / (N_cols + 1), IDright = Id - i,
                       IDrightright = IDright + 1,    // H
        IDrightReIndex = idBasinVectReIndex[IDright]; // H

    // define H at interfaces
    const auto H_interface = H_int_x[Id];

    const double coeff_m = H_interface * alfa_x[Id];

    if (H_interface > H_min) {
      rhs(IDrightReIndex) +=
          isNonReflectingBC * (-c1 * (-coeff_m * u_star[Id]));
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectEast) {
    const unsigned int i = Id / (N_cols + 1), IDleft = Id - i - 1,
                       IDleftleft = IDleft - 1,     // H
        IDleftReIndex = idBasinVectReIndex[IDleft]; // H

    // define H at interfaces
    const auto H_interface = H_int_x[Id];

    const double coeff_m = H_interface * alfa_x[Id];

    if (H_interface > H_min) {
      rhs(IDleftReIndex) += isNonReflectingBC * (-c1 * (+coeff_m * u_star[Id]));
    }
  }

  for (const auto &Id : idStaggeredInternalVectVertical) {
    const unsigned int IDleft = Id - N_cols,          // H
        IDright = Id,                                 // H
        IDleftReIndex = idBasinVectReIndex[IDleft],   // H
        IDrightReIndex = idBasinVectReIndex[IDright]; // H

    // define H at interfaces
    const auto H_interface = H_int_y[Id];

    const double coeff_m = H_interface * alfa_y[Id];

    if (H_interface > H_min) {
      coefficients.push_back(
          Eigen::Triplet<double>(IDleftReIndex, IDrightReIndex, -c3 * coeff_m));
      coefficients.push_back(
          Eigen::Triplet<double>(IDleftReIndex, IDleftReIndex, c3 * coeff_m));

      rhs(IDleftReIndex) +=
          -c1 * (+coeff_m * v_star[Id]) -
          (orography[IDleft] - orography[IDright]) * c3 * coeff_m * isH;

      coefficients.push_back(
          Eigen::Triplet<double>(IDrightReIndex, IDleftReIndex, -c3 * coeff_m));
      coefficients.push_back(
          Eigen::Triplet<double>(IDrightReIndex, IDrightReIndex, c3 * coeff_m));

      rhs(IDrightReIndex) +=
          -c1 * (-coeff_m * v_star[Id]) -
          (orography[IDright] - orography[IDleft]) * c3 * coeff_m * isH;
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectNorth) {
    const unsigned int IDright = Id,
                       IDrightright = Id + N_cols,    //
        IDrightReIndex = idBasinVectReIndex[IDright]; // H

    // define H at interfaces
    const auto H_interface = H_int_y[Id];

    const double coeff_m = H_interface * alfa_y[Id];

    if (H_interface > H_min) {
      rhs(IDrightReIndex) +=
          isNonReflectingBC * (-c1 * (-coeff_m * v_star[Id]));
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectSouth) {
    const unsigned int IDleft = Id - N_cols,        // H
        IDleftleft = IDleft - N_cols,               // H
        IDleftReIndex = idBasinVectReIndex[IDleft]; // H

    // define H at interfaces
    const auto H_interface = H_int_y[Id];

    const double coeff_m = H_interface * alfa_y[Id];

    if (H_interface > H_min) {
      rhs(IDleftReIndex) += isNonReflectingBC * (-c1 * (+coeff_m * v_star[Id]));
    }
  }

  for (const auto &Id : idStaggeredInternalVectHorizontal_not_excluded) {

    const unsigned int i = Id / (N_cols + 1),
                       IDleft = Id - i - 1, // H
        IDright = Id - i;                   // H

    // define H at interfaces
    const auto H_interface = H_int_x[Id];

    if (H_interface > H_min) {
      if (std::get<0>(excluded_ids[IDleft])) {
        const auto &k_pour = std::get<1>(excluded_ids[IDleft]);
        if (k_pour >= 0) {
          additional_source_term[k_pour] += H_interface * std::abs(u[Id]) * c1;
        }
      }

      if (std::get<0>(excluded_ids[IDright])) {
        const auto &k_pour = std::get<1>(excluded_ids[IDright]);
        if (k_pour >= 0) {
          additional_source_term[k_pour] += H_interface * std::abs(u[Id]) * c1;
        }
      }
    }
  }

  for (const auto &Id : idStaggeredInternalVectVertical_not_excluded) {
    const unsigned int IDleft = Id - N_cols, // H
        IDright = Id;                        // H

    // define H at interfaces
    const auto H_interface = H_int_y[Id];

    if (H_interface > H_min) {
      if (std::get<0>(excluded_ids[IDleft])) {
        const auto &k_pour = std::get<1>(excluded_ids[IDleft]);
        if (k_pour >= 0) {
          additional_source_term[k_pour] += H_interface * std::abs(v[Id]) * c1;
        }
      }

      if (std::get<0>(excluded_ids[IDright])) {
        const auto &k_pour = std::get<1>(excluded_ids[IDright]);
        if (k_pour >= 0) {
          additional_source_term[k_pour] += H_interface * std::abs(v[Id]) * c1;
        }
      }
    }
  }

  for (const auto &Id : idBasinVect) {
    const auto IDreIndex = idBasinVectReIndex[Id];
    if (!std::get<0>(excluded_ids[Id])) {
      rhs(IDreIndex) += additional_source_term[Id];
    }
  }
}

//==============================================================================

void updateVel(
    std::vector<double> &u, std::vector<double> &v,
    const std::vector<double> &u_star, const std::vector<double> &v_star,
    const std::vector<double> &alfa_x, const std::vector<double> &alfa_y,
    const double &N_rows, const double &N_cols, const double &c2,
    const double &H_min, const std::vector<double> &eta,
    const std::vector<double> &H, const std::vector<double> &orography,
    const std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    const std::vector<unsigned int> &idStaggeredInternalVectVertical,
    const std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    const std::vector<unsigned int> &idStaggeredBoundaryVectEast,
    const std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    const std::vector<unsigned int> &idStaggeredBoundaryVectSouth,
    const bool &isNonReflectingBC) {

  // +-----------------------------------------------+
  // |              Update Vertical Velocity         |
  // +-----------------------------------------------+

  for (const auto &Id : idStaggeredInternalVectVertical) {
    const unsigned int IDsouth = Id, IDnorth = Id - N_cols;

    const auto &H_interface =
        (H[IDsouth] + H[IDnorth]) * .5 +
        signum(-eta[IDsouth] + eta[IDnorth]) * (-H[IDsouth] + H[IDnorth]) * .5;
    if (H_interface > H_min) {
      v[Id] = alfa_y[Id] * (v_star[Id] - c2 * (eta[IDsouth] - eta[IDnorth]));
    } else {
      v[Id] = 0.;
    }
  }

  // first row
  for (const auto &Id : idStaggeredBoundaryVectNorth) {
    const unsigned int IDsouth = Id + N_cols, IDnorth = Id;

    const auto &H_interface =
        (H[IDnorth]) * .5 +
        signum(-eta[IDsouth] + eta[IDnorth]) * (-H[IDnorth]) * .5;
    if (H_interface > H_min) {
      v[Id] = isNonReflectingBC * alfa_y[Id] * (v_star[Id]);
    } else {
      v[Id] = 0.;
    }
  }

  // last row
  for (const auto &Id : idStaggeredBoundaryVectSouth) {

    const unsigned int IDsouth = Id - N_cols, IDnorth = Id - 2 * N_cols;

    const auto &H_interface =
        (H[IDsouth]) * .5 +
        signum(-eta[IDsouth] + eta[IDnorth]) * (H[IDsouth]) * .5;
    if (H_interface > H_min) {
      v[Id] = isNonReflectingBC * alfa_y[Id] * (v_star[Id]);
    } else {
      v[Id] = 0.;
    }
  }

  // +-----------------------------------------------+
  // |              Update Horizontal Velocity       |
  // +-----------------------------------------------+

  for (const auto &Id : idStaggeredInternalVectHorizontal) {

    const unsigned int i = Id / (N_cols + 1), IDeast = Id - i,
                       IDwest = Id - i - 1;

    const auto &H_interface =
        (H[IDeast] + H[IDwest]) * .5 +
        signum(-eta[IDeast] + eta[IDwest]) * (-H[IDeast] + H[IDwest]) * .5;
    if (H_interface > H_min) {
      u[Id] = alfa_x[Id] * (u_star[Id] - c2 * (eta[IDeast] - eta[IDwest]));
    } else {
      u[Id] = 0.;
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectWest) {

    const unsigned int i = Id / (N_cols + 1), IDeast = Id - i + 1,
                       IDwest = Id - i;

    const auto &H_interface =
        (H[IDwest]) * .5 +
        signum(-eta[IDeast] + eta[IDwest]) * (-H[IDwest]) * .5;
    if (H_interface > H_min) {
      u[Id] = isNonReflectingBC * alfa_x[Id] * (u_star[Id]);
    } else {
      u[Id] = 0.;
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectEast) {

    const unsigned int i = Id / (N_cols + 1), IDeast = Id - i - 1,
                       IDwest = Id - i - 2;

    const auto &H_interface =
        (H[IDeast]) * .5 +
        signum(-eta[IDeast] + eta[IDwest]) * (H[IDeast]) * .5;
    if (H_interface > H_min) {
      u[Id] = isNonReflectingBC * alfa_x[Id] * (u_star[Id]);
    } else {
      u[Id] = 0.;
    }
  }
}

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

double maxdt(const std::vector<double> &u, const std::vector<double> &v,
             const double &Hmax, const double &pixel_size) {

  static constexpr double gravity = 9.81;
  // +-----------------------------------------------+
  // |      Estimate vertical max Courant number     |
  // +-----------------------------------------------+

  const double vel_max_y =
      std::max(*std::max_element(v.begin(), v.end()),
               std::abs(*std::min_element(v.begin(), v.end())));

  // +-----------------------------------------------+
  // |    Estimate horizontal max Courant number     |
  // +-----------------------------------------------+

  const double vel_max_x =
      std::max(*std::max_element(u.begin(), u.end()),
               std::abs(*std::min_element(u.begin(), u.end())));

  const double Co = 2.;      // 0.9; // 0.3
  const double Co_cel = 1e4; // 10

  const double cel = std::sqrt(Hmax * gravity);

  double dt_candidate =
      Co * pixel_size /
      (std::max(vel_max_x, vel_max_y) + std::numeric_limits<double>::epsilon());
  dt_candidate = std::min(dt_candidate,
                          Co_cel * pixel_size /
                              (cel + std::numeric_limits<double>::epsilon()));

  return (dt_candidate);
}

//==============================================================================

double maxCourant(const std::vector<double> &u, const std::vector<double> &v,
                  const double &c1) {

  // +-----------------------------------------------+
  // |      Estimate vertical max Courant number     |
  // +-----------------------------------------------+

  const double Courant_y =
      std::max(*std::max_element(v.begin(), v.end()),
               std::abs(*std::min_element(v.begin(), v.end())));

  // +-----------------------------------------------+
  // |    Estimate horizontal max Courant number     |
  // +-----------------------------------------------+

  const double Courant_x =
      std::max(*std::max_element(u.begin(), u.end()),
               std::abs(*std::min_element(u.begin(), u.end())));

  // output the result
  return (std::max(Courant_y, Courant_x) * c1);
}

//==============================================================================

double maxCourant(const std::vector<double> &H, const double &c1) {
  static constexpr double gravity = 9.81;
  const double Courant_cel =
      std::sqrt(*std::max_element(H.begin(), H.end()) * gravity);
  return (Courant_cel * c1);
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

  return (dt_sed);
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
  return (result);
}

//==============================================================================
// For gravitational layer
void computeResiduals(
    const std::vector<double> &n_x, const std::vector<double> &n_y,
    const unsigned int &N_cols, const unsigned int &N_rows,
    const std::vector<double> &h,
    const std::vector<double> &coeff, // hydraulic conductivity
    const std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    const std::vector<unsigned int> &idStaggeredInternalVectVertical,
    const std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    const std::vector<unsigned int> &idStaggeredBoundaryVectEast,
    const std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    const std::vector<unsigned int> &idStaggeredBoundaryVectSouth,
    const std::vector<unsigned int> &idBasinVect,
    std::vector<double> &h_interface_x, std::vector<double> &h_interface_y,
    std::vector<double> &Res_x, std::vector<double> &Res_y) {

  // +-----------------------------------------------+
  // |                  Horizontal                   |
  // +-----------------------------------------------+

  for (const auto &Id : idStaggeredInternalVectHorizontal) {
    const unsigned int i = Id / (N_cols + 1), // u
        IDeast = Id - i,                      // H
        IDwest = Id - i - 1;                  // H

    const double &h_left = h[IDwest], &h_right = h[IDeast];

    const double k_c_left = coeff[IDwest], k_c_right = coeff[IDeast];

    h_interface_x[Id] = n_x[Id] *
                        ((k_c_left * h_left + k_c_right * h_right) +
                         n_x[Id] * (k_c_left * h_left - k_c_right * h_right)) *
                        .5;
  }

  for (const auto &Id : idStaggeredBoundaryVectWest) {
    const unsigned int i = Id / (N_cols + 1);
    const double h_left = 0, h_right = h[Id - i];
    const double k_c_left = 0., k_c_right = coeff[Id - i];
    h_interface_x[Id] = n_x[Id] *
                        ((k_c_left * h_left + k_c_right * h_right) +
                         n_x[Id] * (k_c_left * h_left - k_c_right * h_right)) *
                        .5;
  }

  for (const auto &Id : idStaggeredBoundaryVectEast) {
    const unsigned int i = Id / (N_cols + 1);
    const double h_left = h[Id - i - 1], h_right = 0;
    const double k_c_left = coeff[Id - i - 1], k_c_right = 0.;
    h_interface_x[Id] = n_x[Id] *
                        ((k_c_left * h_left + k_c_right * h_right) +
                         n_x[Id] * (k_c_left * h_left - k_c_right * h_right)) *
                        .5;
  }

  for (const auto &Id : idBasinVect) {
    const unsigned int i = Id / N_cols;
    Res_x[Id] = h_interface_x[Id + 1 + i] - h_interface_x[Id + i];
  }

  // +-----------------------------------------------+
  // |                   Vertical                    |
  // +-----------------------------------------------+

  for (const auto &Id : idStaggeredInternalVectVertical) {
    const unsigned int IDsouth = Id, // H
        IDnorth = Id - N_cols;       // H

    const double h_left = h[IDnorth], h_right = h[IDsouth];

    const double k_c_left = coeff[IDnorth], k_c_right = coeff[IDsouth];

    h_interface_y[Id] = n_y[Id] *
                        ((k_c_left * h_left + k_c_right * h_right) +
                         n_y[Id] * (k_c_left * h_left - k_c_right * h_right)) *
                        .5;
  }

  for (const auto &Id : idStaggeredBoundaryVectNorth) {
    const double h_left = 0, h_right = h[Id];

    const double k_c_left = 0., k_c_right = coeff[Id];

    h_interface_y[Id] = n_y[Id] *
                        ((k_c_left * h_left + k_c_right * h_right) +
                         n_y[Id] * (k_c_left * h_left - k_c_right * h_right)) *
                        .5;
  }

  for (const auto &Id : idStaggeredBoundaryVectSouth) {
    const double h_left = h[Id - N_cols], h_right = 0;

    const double k_c_left = coeff[Id - N_cols], k_c_right = 0.;

    h_interface_y[Id] = n_y[Id] *
                        ((k_c_left * h_left + k_c_right * h_right) +
                         n_y[Id] * (k_c_left * h_left - k_c_right * h_right)) *
                        .5;
  }

  for (const auto &Id : idBasinVect) {
    Res_y[Id] = h_interface_y[Id + N_cols] - h_interface_y[Id];
  }
}

//==============================================================================
// For sediment transport
void computeResidualsTruncated(
    const std::vector<double> &u, const std::vector<double> &v,
    const unsigned int &N_cols, const unsigned int &N_rows,
    const unsigned int &N, const double &c1, const std::vector<double> &S_x,
    const std::vector<double> &S_y, const double &alpha, const double &beta,
    const double &gamma,
    const std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    const std::vector<unsigned int> &idStaggeredInternalVectVertical,
    const std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    const std::vector<unsigned int> &idStaggeredBoundaryVectEast,
    const std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    const std::vector<unsigned int> &idStaggeredBoundaryVectSouth,
    std::vector<std::array<double, 2>> &Gamma_x,
    std::vector<std::array<double, 2>> &Gamma_y) {

  for (unsigned int ii = 0; ii < idStaggeredInternalVectHorizontal.size(); ii++) {
    const auto &Id = idStaggeredInternalVectHorizontal[ii];

    const double coeff_right = c1 * alpha * std::pow(std::abs(S_x[Id]), beta) *
                               u[Id] * (.5 - .5 * signum(u[Id])),

                 coeff_left = c1 * alpha * std::pow(std::abs(S_x[Id]), beta) *
                              u[Id] * (.5 + .5 * signum(u[Id]));

    Gamma_x[Id][0] = coeff_right;
    Gamma_x[Id][1] = coeff_left;
  }

  for (unsigned int ii = 0; ii < idStaggeredBoundaryVectWest.size(); ii++) {
    const auto &Id = idStaggeredBoundaryVectWest[ii];

    const double coeff_right = c1 * alpha * std::pow(std::abs(S_x[Id]), beta) *
                               u[Id] * (.5 - .5 * signum(u[Id])),

                 coeff_left = c1 * alpha * std::pow(std::abs(S_x[Id]), beta) *
                              u[Id] * (.5 + .5 * signum(u[Id]));

    Gamma_x[Id][0] = coeff_right;
    Gamma_x[Id][1] = coeff_left;
  }

  for (unsigned int ii = 0; ii < idStaggeredInternalVectVertical.size(); ii++) {
    const auto &Id = idStaggeredInternalVectVertical[ii];

    const double coeff_right = c1 * alpha * std::pow(std::abs(S_y[Id]), beta) *
                               v[Id] * (.5 - .5 * signum(v[Id])),

                 coeff_left = c1 * alpha * std::pow(std::abs(S_y[Id]), beta) *
                              v[Id] * (.5 + .5 * signum(v[Id]));

    Gamma_y[Id] = std::array<double, 2>{{coeff_right, coeff_left}};
  }

  for (unsigned int ii = 0; ii < idStaggeredBoundaryVectNorth.size(); ii++) {
    const auto &Id = idStaggeredBoundaryVectNorth[ii];

    const double coeff_right = c1 * alpha * std::pow(std::abs(S_y[Id]), beta) *
                               v[Id] * (.5 - .5 * signum(v[Id])),

                 coeff_left = c1 * alpha * std::pow(std::abs(S_y[Id]), beta) *
                              v[Id] * (.5 + .5 * signum(v[Id]));

    Gamma_y[Id] = std::array<double, 2>{{coeff_right, coeff_left}};
  }

  for (unsigned int ii = 0; ii < idStaggeredBoundaryVectSouth.size(); ii++) {
    const auto &Id = idStaggeredBoundaryVectSouth[ii];

    const double coeff_right = c1 * alpha * std::pow(std::abs(S_y[Id]), beta) *
                               v[Id] * (.5 - .5 * signum(v[Id])),

                 coeff_left = c1 * alpha * std::pow(std::abs(S_y[Id]), beta) *
                              v[Id] * (.5 + .5 * signum(v[Id]));

    Gamma_y[Id] = std::array<double, 2>{{coeff_right, coeff_left}};
  }
}

//==============================================================================
