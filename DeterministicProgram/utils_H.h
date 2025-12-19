#pragma once

//! std library
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#ifndef ENABLE_CUDA

//! Eigen library
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#else

#include "cuda_utils_loop_H.cuh"

#endif

//! Parse library
#include "GetPot.hpp"

//! IML++ CG template
#include "cg.hpp"

namespace Eigen {

namespace internal {

template <typename Scalar>
inline void putVectorElt(Scalar value, std::ofstream &out,
                         const unsigned int &i) {
  out << i << " " << value << "\n";
}
template <typename Scalar>
inline void putVectorElt(std::complex<Scalar> value, std::ofstream &out,
                         const unsigned int &i) {
  out << i << " " << value.real << " " << value.imag() << "\n";
}

} // namespace internal

template <typename VectorType>
bool saveMarketVector_lis(const VectorType &vec, const std::string &filename) {
  typedef typename VectorType::Scalar Scalar;
  std::ofstream out(filename.c_str(), std::ios::out);
  if (!out)
    return false;

  out.flags(std::ios_base::scientific);
  out.precision(64);
  if (internal::is_same<Scalar, std::complex<float>>::value ||
      internal::is_same<Scalar, std::complex<double>>::value)
    out << "%%MatrixMarket vector coordinate complex general\n";
  else
    out << "%%MatrixMarket vector coordinate real general\n";
  out << vec.size() << "\n";
  for (int i = 0; i < vec.size(); i++) {
    internal::putVectorElt(vec(i), out, i + 1);
  }
  out.close();
  return true;
}

} // namespace Eigen

class Vector2D {

public:
  //! Empty constructor (all components are set to zero)
  Vector2D() {}

  Vector2D(
      std::array<double, 2> const &indices) // Note the use of double indices
      : M_coords(indices) {}

  //! Copy constructor
  Vector2D(Vector2D const &vector) { *this = vector; }

  ~Vector2D() = default;

  //! Operator +=
  Vector2D &operator+=(Vector2D const &vector) {
    for (unsigned int i = 0; i < 2; i++)
      M_coords[i] += vector.M_coords[i];
    return *this;
  }

  //! Assignment operator
  Vector2D &operator=(Vector2D const &vector) {
    for (unsigned int i = 0; i < 2; i++)
      M_coords[i] = vector.M_coords[i];
    return *this;
  }

  Vector2D operator+(Vector2D const &vector) const {
    Vector2D tmp(*this);
    return tmp += vector;
  }

  //! Operator -=
  Vector2D &operator-=(Vector2D const &vector) {
    for (unsigned int i = 0; i < 2; i++)
      M_coords[i] -= vector.M_coords[i];
    return *this;
  }

  //! Operator -
  Vector2D operator-(Vector2D const &vector) const {
    Vector2D tmp(*this);
    return tmp -= vector;
  }

  //! Operator *= (multiplication by scalar)
  Vector2D &operator*=(double const &factor) {
    for (unsigned int i = 0; i < 2; i++)
      M_coords[i] *= factor;
    return *this;
  }

  //! Operator /= (division by scalar)
  Vector2D &operator/=(double const &factor) {
    *this *= 1. / factor;
    return *this;
  }

  //! Operator / (division by scalar)
  Vector2D operator/(double const &factor) const {
    Vector2D tmp(*this);
    return tmp /= factor;
  }

  Vector2D operator*(double const &factor) {
    Vector2D tmp(*this);
    return tmp *= factor;
  }

  double dot(Vector2D const &vector) const {
    double scalarProduct = 0.;
    for (unsigned int i = 0; i < 2; i++)
      scalarProduct += M_coords[i] * vector.M_coords[i];
    return scalarProduct;
  }

  double norm() const { return std::sqrt(this->dot(*this)); }

  //! Operator ()b
  double const &operator()(unsigned int const &i) const { return M_coords[i]; }

  //! Operator ()
  double &operator()(unsigned int const &i) { return M_coords[i]; }

private:
  std::array<double, 2> M_coords;
};

//! Operator * (multiplication by scalar on the right)
Vector2D operator*(Vector2D const &vector, double const &factor);

//! Operator * (multiplication by scalar on the left)
Vector2D operator*(double const &factor, Vector2D const &vector);

std::map<int, std::array<double, 2>> createCN_map_Gav(const std::string &file);

std::map<std::array<int, 2>, int> createCN_map();

//==============================================================================

class Raster {

public:
  Raster(const std::string &file);

  ~Raster() = default;

  unsigned int ncols, nrows;

  double xllcorner, yllcorner, cellsize, NODATA_value;

  Eigen::SparseMatrix<double> Coords; // forse mettere una matrice densa
};

//==============================================================================

double signum(const double &x);

//==============================================================================

class Rain // Previous interpolation not Linear
{

public:
  Rain(const std::string &infiltrationModel, const unsigned int &N,
       const bool &isInitialLoss, const double &perc_initialLoss,
       const bool is_precipitation, const bool constant_precipitation,
       const std::string precipitation_file, const std::string file_dir,
       const double time_spacing_rain, const int number_stations,
       const double max_Days, const GetPot &dataFile, const double &xllcorner,
       const double &yllcorner, const double &pixel_size,
       const unsigned int &N_rows, const unsigned int &N_cols,
       const std::vector<unsigned int> &idBasinVect);

  Rain() = delete;
  ~Rain() = default;

  void constant_precipitation(const std::string &file,
                              const unsigned int &ndata,
                              const bool &is_precipitation,
                              const double &time_spacing);

  void IDW_precipitation(const std::vector<std::string> &file_vect,
                         const std::vector<unsigned int> &ndata_vec,
                         const std::vector<double> &time_spacing_vect,
                         const std::vector<double> &X,
                         const std::vector<double> &Y, const double &xllcorner,
                         const double &yllcorner, const double &pixel_size,
                         const unsigned int &N_rows, const unsigned int &N_cols,
                         const std::vector<unsigned int> &idBasinVect);

  void computePrecipitation(const double &time, const std::vector<double> &S,
                            const std::vector<double> &melt_mask,
                            const std::vector<double> &h_G,
                            const std::vector<double> &H,
                            const unsigned int &N_rows,
                            const unsigned int &N_cols,
                            const std::vector<unsigned int> &idBasinVect);

  std::vector<double> DP_total, DP_cumulative, DP_infiltrated;

  double dt_rain;

private:
  std::vector<std::vector<double>> Hyetograph, // # station times ndata
      IDW_weights;

  std::vector<double> M_time_spacing_vect;
  bool M_isInitialLoss;
  double rainfall_intensity = 0;
  double c;
};

//==============================================================================

class Temperature {

public:
  Temperature(const std::string &file, const unsigned int &N,
              const unsigned int &max_Days, const double &T_crit,
              const std::vector<double> &orography, const unsigned int &ndata,
              const unsigned int &steps_per_hour, const double &time_spacing,
              const double &height_thermometer, const std::string format_temp);

  Temperature() = delete;
  ~Temperature() = default;

  void computeTemperature(const unsigned int &i,
                          const std::vector<double> &orography,
                          const std::vector<unsigned int> &idBasinVect);

  std::vector<double> T_raster, melt_mask, T_dailyMean, T_dailyMin, T_dailyMax,
      J;

  const double T_crit;

private:
  std::vector<double> Temperature_Graph; // length: ndata
  static constexpr double Temp_diff = -6.5e-3;
  const double height_th;
};

//==============================================================================

class evapoTranspiration {

public:
  evapoTranspiration(const std::string &ET_model, const unsigned int &N,
                     const std::vector<double> &orography,
                     const std::vector<double> &J, const unsigned int &max_Days,
                     const double &phi_rad, const double &height_thermometer);

  evapoTranspiration() = delete;
  ~evapoTranspiration() = default;

  void ET(const std::vector<double>
              &T_mean, // length nstep: vector of temperature in deg Celsius
          const std::vector<double> &T_min, // length nstep
          const std::vector<double> &T_max, // length nstep
          const int &i, const std::vector<unsigned int> &idBasinVect,
          const std::vector<double> &orography);

  std::vector<double> ET_vec;

private:
  std::vector<double> Ra;
  unsigned int M_evapoTranspiration_model;
  static constexpr double M_Gsc = 0.082; // Solar constant
  const double height_th;
  static constexpr double Temp_diff = -6.5e-3;
};

//==============================================================================

template<class T_type, class U_type>
class frictionClass {

public:
  frictionClass(
      const T_type&H_interface_horizontal,
      const T_type&H_interface_vertical,
      const T_type&u, const T_type&v,
      const U_type&idStaggeredInternalVectHorizontal,
      const U_type&idStaggeredBoundaryVectWest,
      const U_type&idStaggeredBoundaryVectEast,
      const U_type&idStaggeredInternalVectVertical,
      const U_type&idStaggeredBoundaryVectNorth,
      const U_type&idStaggeredBoundaryVectSouth,
      const std::string &friction_model, const double &n_manning,
      const double &dt_DSV, const std::vector<double> &d_90,
      const std::vector<double> &rough, const double &H_min,
      const unsigned int &N_rows, const unsigned int &N_cols, 
#ifdef ENABLE_CUDA
      const thrust::host_vector<double> & S_x, const thrust::host_vector<double> &S_y
#else
      const std::vector<double> &S_x, const std::vector<double> &S_y
#endif      
     )  
     : H_interface_horizontal(H_interface_horizontal),
      H_interface_vertical(H_interface_vertical), u(u), v(v),
      idStaggeredInternalVectHorizontal(idStaggeredInternalVectHorizontal),
      idStaggeredBoundaryVectWest(idStaggeredBoundaryVectWest),
      idStaggeredBoundaryVectEast(idStaggeredBoundaryVectEast),
      idStaggeredInternalVectVertical(idStaggeredInternalVectVertical),
      idStaggeredBoundaryVectNorth(idStaggeredBoundaryVectNorth),
      idStaggeredBoundaryVectSouth(idStaggeredBoundaryVectSouth),
      M_n_manning(n_manning), M_dt_DSV(dt_DSV), N_rows(N_rows), N_cols(N_cols) {

    M_H_min = std::pow(H_min, M_expo);

    std::vector<double> M_fc0_lower_x(S_x.size()), M_fc0_greater_x(S_x.size()),
        M_fc0_lower_y(S_y.size()), M_fc0_greater_y(S_y.size());

    for (unsigned int i = 0; i < N_rows; i++) {
      for (unsigned int j = 0; j < N_cols; j++) {
        const unsigned int IDcell = j + i * N_cols, IDleft = IDcell + i,
                           IDright = IDleft + 1, IDup = IDcell,
                           IDdown = IDcell + N_cols;

        const auto d_90_cell = rough[IDcell] * d_90[IDcell];

        M_fc0_greater_x[IDleft] =
            std::pow(d_90_cell, .45) / (.56 * std::pow(M_g, .44));
        M_fc0_lower_x[IDleft] =
            std::pow(d_90_cell, .24) / (2.73 * std::pow(M_g, .49));
        M_fc0_greater_y[IDup] =
            std::pow(d_90_cell, .45) / (.56 * std::pow(M_g, .44));
        M_fc0_lower_y[IDup] =
            std::pow(d_90_cell, .24) / (2.73 * std::pow(M_g, .49));

        M_fc0_greater_x[IDright] =
            std::pow(d_90_cell, .45) / (.56 * std::pow(M_g, .44));
        M_fc0_lower_x[IDright] =
            std::pow(d_90_cell, .24) / (2.73 * std::pow(M_g, .49));
        M_fc0_greater_y[IDdown] =
            std::pow(d_90_cell, .45) / (.56 * std::pow(M_g, .44));
        M_fc0_lower_y[IDdown] =
            std::pow(d_90_cell, .24) / (2.73 * std::pow(M_g, .49));
      }
    }

    alfa_x.resize(S_x.size());
    alfa_y.resize(S_y.size());

    M_coeff = M_g * std::pow(M_n_manning, 2.);

    M_expo_r_x_vect.resize(S_x.size());
    M_gamma_dt_DSV_x_.resize(S_x.size());
    for (unsigned int k = 0; k < S_x.size(); k++) {
      M_expo_r_x_vect[k] = M_expo_r1 * (std::abs(S_x[k]) > .006) +
                           M_expo_r2 * (std::abs(S_x[k]) <= .006);
      const auto M_Rick_x =
          (M_fc0_greater_x[k] * std::pow(std::abs(S_x[k]), .33)) *
                (std::abs(S_x[k]) > .006) +
            (M_fc0_lower_x[k] * std::pow(std::abs(S_x[k]), .08)) *
                (std::abs(S_x[k]) <= .006);
      M_gamma_dt_DSV_x_[k] = M_g * std::pow(M_Rick_x, 2);
    }

    M_expo_r_y_vect.resize(S_y.size());
    M_gamma_dt_DSV_y_.resize(S_y.size());
    for (unsigned int k = 0; k < S_y.size(); k++) {
      M_expo_r_y_vect[k] = M_expo_r1 * (std::abs(S_y[k]) > .006) +
                           M_expo_r2 * (std::abs(S_y[k]) <= .006);
      const auto M_Rick_y =
          (M_fc0_greater_y[k] * std::pow(std::abs(S_y[k]), .33)) *
              (std::abs(S_y[k]) > .006) +
          (M_fc0_lower_y[k] * std::pow(std::abs(S_y[k]), .08)) *
              (std::abs(S_y[k]) <= .006);
      M_gamma_dt_DSV_y_[k] = M_g * std::pow(M_Rick_y, 2);
    }

    if (friction_model == "None") {
      M_frictionModel = 0;
    } else if (friction_model == "Manning") {
      M_frictionModel = 1;
    } else if (friction_model == "Rickenmann") {
      M_frictionModel = 2;
    } else {
      std::cout << "No friction model inserted!!" << std::endl;
      exit(-1.);
    }
  }

  frictionClass() = delete;
  ~frictionClass() = default;

  void f_x() {
#ifdef ENABLE_CUDA
    
    LAUNCH_KERNEL(
      computeHorizontalInternalKernel_friction,
      idStaggeredInternalVectHorizontal.size(),
      thrust::raw_pointer_cast(idStaggeredInternalVectHorizontal.data()),
      thrust::raw_pointer_cast(H_interface_horizontal.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(M_expo_r_x_vect.data()),
      N_cols
    );

    LAUNCH_KERNEL(
      computeHorizontalWestKernel_interface_friction,
      idStaggeredBoundaryVectWest.size(),
      thrust::raw_pointer_cast(idStaggeredBoundaryVectWest.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(horizontal.data()),
      N_cols
    );

    LAUNCH_KERNEL(
      computeHorizontalEastKernel_friction,
      idStaggeredBoundaryVectEast.size(),
      thrust::raw_pointer_cast(idStaggeredBoundaryVectEast.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(horizontal.data()),
      N_cols
    );

#else
    for (const auto &Id : idStaggeredInternalVectHorizontal) {
      double alfa = 1.;

      const auto &H_int = H_interface_horizontal[Id];
      const auto &exponent = M_expo_r_x_vect[Id];
      const auto den =
          std::pow(H_int, M_expo + exponent * (M_frictionModel == 2));

      if (den > M_H_min) {
        const auto u_abs = std::abs(u[Id]);
        double coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * u_abs / den *
                       (M_frictionModel > 0);
        coeff = std::max(
            coeff, M_dt_DSV * M_gamma_dt_DSV_x_[Id] *
                       std::pow(u_abs, 1. - exponent * (M_frictionModel == 2)) /
                       den);
        alfa = 1. / (1. + coeff);
      }
      alfa_x[Id] = alfa;
    }

    for (const auto &Id : idStaggeredBoundaryVectWest) {
      double alfa = 1.;

      const auto &H_int = H_interface_horizontal[Id];
      const auto &exponent = M_expo_r_x_vect[Id];
      const auto den =
          std::pow(H_int, M_expo + exponent * (M_frictionModel == 2));

      if (den > M_H_min) {
        const auto u_abs = std::abs(u[Id]);
        double coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * u_abs / den *
                       (M_frictionModel > 0);
        coeff = std::max(
            coeff, M_dt_DSV * M_gamma_dt_DSV_x_[Id] *
                       std::pow(u_abs, 1. - exponent * (M_frictionModel == 2)) /
                       den);
        alfa = 1. / (1. + coeff);
      }
      alfa_x[Id] = alfa;
    }

    for (const auto &Id : idStaggeredBoundaryVectEast) {
      double alfa = 1.;

      const auto &H_int = H_interface_horizontal[Id];
      const auto &exponent = M_expo_r_x_vect[Id];
      const auto den =
          std::pow(H_int, M_expo + exponent * (M_frictionModel == 2));

      if (den > M_H_min) {
        const auto u_abs = std::abs(u[Id]);
        double coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * u_abs / den *
                       (M_frictionModel > 0);
        coeff = std::max(
            coeff, M_dt_DSV * M_gamma_dt_DSV_x_[Id] *
                       std::pow(u_abs, 1. - exponent * (M_frictionModel == 2)) /
                       den);
        alfa = 1. / (1. + coeff);
      }
      alfa_x[Id] = alfa;
    }
#endif
  }


  void f_y() {
#ifdef ENABLE_CUDA
#error "Not yet implemented"
#else
    for (const auto &Id : idStaggeredInternalVectVertical) {
      double alfa = 1.;

      const auto &H_int = H_interface_vertical[Id];
      const auto &exponent = M_expo_r_y_vect[Id];
      const auto den =
        std::pow(H_int, M_expo + exponent * (M_frictionModel == 2));

    if (den > M_H_min) {
      const auto v_abs = std::abs(v[Id]);

      double coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * v_abs / den *
                     (M_frictionModel > 0);
      coeff = std::max(
          coeff, M_dt_DSV * M_gamma_dt_DSV_y_[Id] *
                     std::pow(v_abs, 1. - exponent * (M_frictionModel == 2)) /
                     den);
      alfa = 1. / (1. + coeff);
    }
    alfa_y[Id] = alfa;
  }

  for (const auto &Id : idStaggeredBoundaryVectNorth) {
    double alfa = 1.;

    const auto &H_int = H_interface_vertical[Id];
    const auto &exponent = M_expo_r_y_vect[Id];
    const auto den =
        std::pow(H_int, M_expo + exponent * (M_frictionModel == 2));

    if (den > M_H_min) {
      const auto v_abs = std::abs(v[Id]);
      double coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * v_abs / den *
                     (M_frictionModel > 0);
        coeff = std::max(
            coeff, M_dt_DSV * M_gamma_dt_DSV_y_[Id] *
                       std::pow(v_abs, 1. - exponent * (M_frictionModel == 2)) /
                       den);
        alfa = 1. / (1. + coeff);
      }
      alfa_y[Id] = alfa;
    }  

    for (const auto &Id : idStaggeredBoundaryVectSouth) {
      double alfa = 1.;

      const auto &H_int = H_interface_vertical[Id];
      const auto &exponent = M_expo_r_y_vect[Id];
      const auto den =
          std::pow(H_int, M_expo + exponent * (M_frictionModel == 2));

      if (den > M_H_min) {
        const auto v_abs = std::abs(v[Id]);
        double coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * v_abs / den *
                       (M_frictionModel > 0);
        coeff = std::max(
            coeff, M_dt_DSV * M_gamma_dt_DSV_y_[Id] *
                       std::pow(v_abs, 1. - exponent * (M_frictionModel == 2)) /
                       den);
        alfa = 1. / (1. + coeff);
      }
      alfa_y[Id] = alfa;
    }
#endif
  }

  T_type alfa_x, alfa_y;

private:
  unsigned int M_frictionModel;

  const double &M_n_manning;
  const double &M_dt_DSV;

  double M_coeff;
  std::function<double(double const &, double const &)> M_gamma_dt_DSV =
      [](double const &dt, double const &cc) { return dt * cc; };

  double M_H_min;

  const T_type&u;
  const T_type&v;
  const T_type&H_interface_horizontal;
  const T_type&H_interface_vertical;

  const U_type&idStaggeredInternalVectHorizontal;
  const U_type&idStaggeredBoundaryVectWest;
  const U_type&idStaggeredBoundaryVectEast;
  const U_type&idStaggeredInternalVectVertical;
  const U_type&idStaggeredBoundaryVectNorth;
  const U_type&idStaggeredBoundaryVectSouth;

  const unsigned int &N_rows;
  const unsigned int &N_cols;

  static constexpr double M_expo = 4. / 3.;
  static constexpr double M_expo_r1 = .11 * 2;
  static constexpr double M_expo_r2 = .03 * 2;
  static constexpr double M_g = 9.81;
  static constexpr double M_toll = 1.e-4;

  T_type M_expo_r_x_vect, M_expo_r_y_vect, M_gamma_dt_DSV_x_,
      M_gamma_dt_DSV_y_;
};

//==============================================================================

template<class T_type, class U_type>
class upwind {

public:
  upwind(const T_type&H, const T_type&u,
         const T_type&v,

         const U_type&idStaggeredInternalVectHorizontal,
         const U_type&idStaggeredBoundaryVectWest,
         const U_type&idStaggeredBoundaryVectEast,
         const U_type&idStaggeredInternalVectVertical,
         const U_type&idStaggeredBoundaryVectNorth,
         const U_type&idStaggeredBoundaryVectSouth,

         const unsigned int &N_rows, const unsigned int &N_cols)
      : H(H), u(u), v(v),
        idStaggeredInternalVectHorizontal(idStaggeredInternalVectHorizontal),
        idStaggeredBoundaryVectWest(idStaggeredBoundaryVectWest),
        idStaggeredBoundaryVectEast(idStaggeredBoundaryVectEast),
        idStaggeredInternalVectVertical(idStaggeredInternalVectVertical),
        idStaggeredBoundaryVectNorth(idStaggeredBoundaryVectNorth),
        idStaggeredBoundaryVectSouth(idStaggeredBoundaryVectSouth),
        N_cols(N_cols), N_rows(N_rows) {
    horizontal.resize(u.size());
    vertical.resize(v.size());
  }

  upwind() = delete;
  ~upwind() = default;

  void computeHorizontal() {
#ifdef ENABLE_CUDA

    LAUNCH_KERNEL(
      computeHorizontalInternalKernel_interface,
      idStaggeredInternalVectHorizontal.size(),
      thrust::raw_pointer_cast(idStaggeredInternalVectHorizontal.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(horizontal.data()),
      N_cols
    );

    LAUNCH_KERNEL(
      computeHorizontalWestKernel_interface,
      idStaggeredBoundaryVectWest.size(),
      thrust::raw_pointer_cast(idStaggeredBoundaryVectWest.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(horizontal.data()),
      N_cols
    );

    LAUNCH_KERNEL(
      computeHorizontalEastKernel_interface,
      idStaggeredBoundaryVectEast.size(),
      thrust::raw_pointer_cast(idStaggeredBoundaryVectEast.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(horizontal.data()),
      N_cols
    );

#else

    for (const auto &Id : idStaggeredInternalVectHorizontal) {
      const unsigned int i = Id / (N_cols + 1), 
          IDeast = Id - i,    
          IDwest = Id - i - 1; 
      const double &H_left = H[IDwest], &H_right = H[IDeast];
      horizontal[Id] =
          (H_left + H_right) * .5 + signum(u[Id]) * (H_left - H_right) * .5;
    }
    for (const auto &Id : idStaggeredBoundaryVectWest) {
      const unsigned int i = Id / (N_cols + 1);
      const double H_left = 0, &H_right = H[Id - i];
      horizontal[Id] = (H_left + H_right) * .5 +
                        signum(u[Id + 1]) * (H_left - H_right) *
                           .5; 
    }
    for (const auto &Id : idStaggeredBoundaryVectEast) {
      const unsigned int i = Id / (N_cols + 1);
      const double &H_left = H[Id - i - 1], H_right = 0;
      horizontal[Id] = (H_left + H_right) * .5 +
                       signum(u[Id - 1]) * (H_left - H_right) *
                           .5; 
    }

#endif
  }

  void computeVertical() {
#ifdef ENABLE_CUDA

    LAUNCH_KERNEL(
      computeVerticalInternalKernel_interface,
      idStaggeredInternalVectVertical.size(),
      thrust::raw_pointer_cast(idStaggeredInternalVectVertical.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(vertical.data()),
      N_cols
    );

    LAUNCH_KERNEL(
      computeVerticalNorthKernel_interface,
      idStaggeredBoundaryVectNorth.size(),
      thrust::raw_pointer_cast(idStaggeredBoundaryVectNorth.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(vertical.data()),
      N_cols
    );

    LAUNCH_KERNEL(
      computeVerticalSouthKernel_interface,
      idStaggeredBoundaryVectSouth.size(),
      thrust::raw_pointer_cast(idStaggeredBoundaryVectSouth.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(vertical.data()),
      N_cols
    );


#else

  for (const auto &Id : idStaggeredInternalVectVertical) {
      const auto IDsouth = Id,   
          IDnorth = Id - N_cols;
      const double &H_left = H[IDnorth], &H_right = H[IDsouth];
      vertical[Id] =
          (H_left + H_right) * .5 + signum(v[Id]) * (H_left - H_right) * .5;
    }
    for (const auto &Id : idStaggeredBoundaryVectNorth) {
      const double H_left = 0, &H_right = H[Id];
      vertical[Id] = (H_left + H_right) * .5 +
                     signum(v[Id + N_cols]) * (H_left - H_right) * .5;
    }
    for (const auto &Id : idStaggeredBoundaryVectSouth) {
      const double &H_left = H[Id - N_cols], H_right = 0;
      vertical[Id] = (H_left + H_right) * .5 +
                     signum(v[Id - N_cols]) * (H_left - H_right) * .5;
    }

#endif    
  }

  T_type horizontal, vertical;

private:
  const T_type&H;
  const T_type&u;
  const T_type&v;

  const U_type&idStaggeredInternalVectHorizontal;
  const U_type&idStaggeredBoundaryVectWest;
  const U_type&idStaggeredBoundaryVectEast;
  const U_type&idStaggeredInternalVectVertical;
  const U_type&idStaggeredBoundaryVectNorth;
  const U_type&idStaggeredBoundaryVectSouth;

  const unsigned int &N_cols;
  const unsigned int &N_rows;
};

//==============================================================================

bool is_file_exist(const char *fileName);

void bilinearInterpolation(const std::vector<double> &u,
                           const std::vector<double> &v,
                           std::vector<double> &u_star,
                           std::vector<double> &v_star,
                           const unsigned int &nrows, const unsigned int &ncols,
                           const double &dt, const double &pixel_size);

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
    const std::vector<unsigned int> &idStaggeredBoundaryVectSouth);

double bilinearInterpolation(const std::vector<double> &H,
                             const unsigned int &ncols,
                             const unsigned int &nrows,
                             const Vector2D &XX_gauges);

double bilinearInterpolation(const std::vector<double> &u,
                             const std::vector<double> &v,
                             const unsigned int &ncols,
                             const unsigned int &nrows,
                             const Vector2D &XX_gauges);

int computePourCell(const int &IDcell, const unsigned int &N_cols,
                    const std::vector<double> &oro,
                    const std::set<unsigned int> &idBasinVect,
                    const std::set<unsigned int> &idStaggeredBoundaryVectSouth,
                    const std::set<unsigned int> &idStaggeredBoundaryVectNorth,
                    const std::set<unsigned int> &idStaggeredBoundaryVectWest,
                    const std::set<unsigned int> &idStaggeredBoundaryVectEast);

void computeAdjacencies(
    const std::vector<double> &basin_mask_Vec_mpi,
    const std::vector<double> &basin_mask_Vec,

    std::vector<unsigned int> &idStaggeredBoundaryVectSouth_mpi,
    std::vector<unsigned int> &idStaggeredBoundaryVectNorth_mpi,
    std::vector<unsigned int> &idStaggeredBoundaryVectWest_mpi,
    std::vector<unsigned int> &idStaggeredBoundaryVectEast_mpi,

    std::vector<unsigned int> &idStaggeredInternalVectHorizontal_mpi,
    std::vector<unsigned int> &idStaggeredInternalVectVertical_mpi,

    std::vector<unsigned int> &idBasinVectReIndex_mpi,
    std::vector<unsigned int> &idBasinVectReIndex,

    const unsigned int &N_rows, const unsigned int &N_cols);

void computeAdjacencies(
    const std::vector<double> &basin_mask_Vec,

    std::vector<unsigned int> &idStaggeredBoundaryVectSouth,
    std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    std::vector<unsigned int> &idStaggeredBoundaryVectEast,

    std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    std::vector<unsigned int> &idStaggeredInternalVectVertical,

    std::vector<unsigned int> &idBasinVect,
    std::vector<unsigned int> &idBasinVectReIndex,

    const unsigned int &N_rows, const unsigned int &N_cols);

void computeAdjacencies(
    const std::vector<double> &basin_mask_Vec_input,
    const std::vector<std::tuple<bool, int>> &excluded_ids,

    std::vector<unsigned int> &idStaggeredBoundaryVectSouth,
    std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    std::vector<unsigned int> &idStaggeredBoundaryVectEast,

    std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    std::vector<unsigned int> &idStaggeredInternalVectVertical,

    std::vector<unsigned int> &idBasinVect,
    std::vector<unsigned int> &idBasinVectReIndex,

    const unsigned int &N_rows, const unsigned int &N_cols);

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

    const std::vector<std::tuple<bool, int>> &excluded_ids,
    std::vector<double> &additional_source_term,

    std::vector<Eigen::Triplet<double>> &coefficients, Eigen::VectorXd &rhs);

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
    const bool &isNonReflectingBC);

void putDry_excludedNodes(
    const std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    const std::vector<unsigned int> &idStaggeredInternalVectVertical,
    const std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    const std::vector<unsigned int> &idStaggeredBoundaryVectEast,
    const std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    const std::vector<unsigned int> &idStaggeredBoundaryVectSouth,
    const std::vector<unsigned int> &idBasinVect, const unsigned int &N_cols,
    const std::vector<std::tuple<bool, int>> &excluded_ids,

    Eigen::VectorXd &H, Eigen::VectorXd &eta,
    const std::vector<double> &orography, std::vector<double> &u,
    std::vector<double> &v);

void compute_dt_adaptive(const std::vector<double> &H,
                         const std::vector<double> &H_old,
                         const std::vector<double> &H_oldold,
                         const std::vector<unsigned int> &idBasinVect,
                         double &dt,
                         const double &local_estimator_time_tolerance,
                         const double &time, const double &timed,
                         const double &timedd);

double maxdt(const std::vector<double> &u, const std::vector<double> &v,
             const double &Hmax, const double &pixel_size);

double maxCourant(const std::vector<double> &u, const std::vector<double> &v,
                  const double &c1);

double maxCourant(const std::vector<double> &H, const double &c1);

double compute_dt_sediment(const double alpha, const double beta,
                           const double S_x, const double S_y,
                           const std::vector<double> &u,
                           const std::vector<double> &v,
                           const double pixel_size, const double dt_DSV,
                           unsigned int *numberOfSteps);

int current_start_chunk(const int &rank,
                        const std::vector<int> &chunk_length_vec);

void saveVector(const Eigen::VectorXd &b, const std::string &Name);

void saveMatrix(const Eigen::SparseMatrix<double> &A, const std::string &Name);

void saveSolution(const std::string &preName, const std::string &flag,
                  const unsigned int &N_rows, const unsigned int &N_cols,
                  const double &xllcorner, const double &yllcorner,
                  const double &cellsize, const double &NODATA_value,
                  const Eigen::VectorXd &H); // it is H or orography

void saveSolution(const std::string &preName, const std::string &flag,
                  const unsigned int &N_rows, const unsigned int &N_cols,
                  const double &xllcorner, const double &yllcorner,
                  const double &cellsize, const double &NODATA_value,
                  const std::vector<double> &H); // it is H or orography

void saveSolution(const std::string &preName, const std::string &flag,
                  const unsigned int &N_rows, const unsigned int &N_cols,
                  const double &xllcorner, const double &yllcorner,
                  const double &cellsize, const double &NODATA_value,
                  const std::vector<int> &H); // it is H or orography

void saveSolution(const std::string &preName, const std::string &flag,
                  const unsigned int &N_rows, const unsigned int &N_cols,
                  const double &xllcorner, const double &yllcorner,
                  const double &cellsize, const double &NODATA_value,
                  const unsigned int &n, const std::vector<double> &u,
                  const std::vector<double> &v,
                  const Eigen::VectorXd &H); // it is H or orography

void saveSolution(const std::string &preName, const unsigned int &N_rows,
                  const unsigned int &N_cols, const double &xllcorner,
                  const double &yllcorner, const double &cellsize,
                  const double &NODATA_value,
                  const std::vector<std::tuple<bool, int>>
                      excluded_ids); // excluded regions, high slopes I hope

void saveSolution(const std::string &preName, const std::string &flag,
                  const unsigned int &N_rows, const unsigned int &N_cols,
                  const double &xllcorner, const double &yllcorner,
                  const double &cellsize, const double &NODATA_value,
                  const unsigned int &n, const std::vector<double> &u,
                  const std::vector<double> &v,
                  const std::vector<double> &H); // it is H or orography

void saveTemporalSequence(const Vector2D &X_gauges, const double &time,
                          const std::string &preName, const double &H);

void saveTemporalSequence(const double &time, const std::string &preName,
                          const double &H);

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
    std::vector<double> &Res_x, std::vector<double> &Res_y);

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
    std::vector<std::array<double, 2>> &Gamma_y);

std::vector<double> compute_d_perc(const std::vector<double> &clay,
                                   const std::vector<double> &sand,
                                   const double &perc);
