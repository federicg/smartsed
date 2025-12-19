#pragma once

//! std library
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <complex>

//! pure CPU header
#include "code_init.h"

//! Eigen library
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/SparseCholesky>

//! Include GPU impl.
#ifdef ENABLE_CUDA
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

//==============================================================================

template <class T>
T signum(const T x) {
  return ((x > 0) ? 1.0 : (x < 0) ? -1.0 : 0.0);
}

//==============================================================================

template <class U_type>
class Rain // Previous interpolation not Linear
{
public:
#ifdef ENABLE_CUDA
  using T_type = thrust::device_vector<double>;
#else
  using T_type = std::vector<double>;
#endif

  Rain(const std::string &infiltrationModel, const unsigned int N,
       const bool isInitialLoss, const double perc_initialLoss,
       const bool is_precipitation, const bool constant_precipitation,
       const std::string precipitation_file, const std::string file_dir,
       const double time_spacing_rain, const int number_stations,
       const double max_Days, const GetPot &dataFile, const double xllcorner,
       const double yllcorner, const double pixel_size,
       const unsigned int N_rows, const unsigned int N_cols,
       const U_type &idBasinVect) {

    M_isInitialLoss = isInitialLoss;
    c = perc_initialLoss;

    if (infiltrationModel != "None" && infiltrationModel != "SCS-CN") {
      std::cout << "Insert a valid infiltration model, STOP!" << std::endl;
      exit(1.);
    }

    DP_total.resize(N);
    DP_cumulative.resize(N);
    DP_infiltrated.resize(N);
    IDW_weights.resize(N);

    dt_rain = 0;

    if (constant_precipitation || !is_precipitation) {
      dt_rain = time_spacing_rain * 3600;

      const auto ndata_rain = std::round(max_Days * 24 / time_spacing_rain);

      this->constant_precipitation(file_dir + precipitation_file, ndata_rain,
                                 is_precipitation, time_spacing_rain);
    } else // IDW
    {
      std::vector<std::string> precipitation_file;
      std::vector<double> time_spacing_rain, X, Y;
      std::vector<unsigned int> ndata_rain;

      for (int number = 1; number <= number_stations; number++) {
        std::string filename = "files/meteo_data/rain_file_";
        filename += std::to_string(number);
        const std::string precipitation_file_current =
            dataFile(filename.c_str(), " ");

        filename = "files/meteo_data/time_spacing_rain_";
        filename += std::to_string(number);
        const double time_spacing_rain_current = dataFile(filename.c_str(), 1.);

        filename = "files/meteo_data/X_";
        filename += std::to_string(number);
        const double X_current = dataFile(filename.c_str(), 1.);

        filename = "files/meteo_data/Y_";
        filename += std::to_string(number);
        const double Y_current = dataFile(filename.c_str(), 1.);

        const unsigned int ndata_rain_current =
            std::round(max_Days * 24 / time_spacing_rain_current);

        precipitation_file.push_back(file_dir + precipitation_file_current);
        time_spacing_rain.push_back(time_spacing_rain_current);
        X.push_back(X_current);
        Y.push_back(Y_current);
        ndata_rain.push_back(ndata_rain_current);
      }

      dt_rain =
          *std::min_element(time_spacing_rain.begin(), time_spacing_rain.end()) *
          3600;

      this->IDW_precipitation(precipitation_file, ndata_rain, time_spacing_rain,
                              X, Y, xllcorner, yllcorner, pixel_size, N_rows,
                              N_cols, idBasinVect);
    }

  }

  Rain() = delete;
  ~Rain() = default;

  void constant_precipitation(const std::string &file,
                              const unsigned int ndata,
                              const bool is_precipitation,
                              const double time_spacing) {
    M_time_spacing_vect.resize(1);
    M_time_spacing_vect[0] = time_spacing;

    Hyetograph.resize(1);

    if (is_precipitation) {
      std::ifstream ff(file);
      if (ff.is_open()) {
        std::string str;
        std::getline(ff, str);
        for (unsigned int i = 0; i < ndata; i++) {
          std::string str1;
          ff >> str1;

          std::string hour;
          ff >> hour;

          std::string hour_string(hour.begin(), hour.begin() + 2),
              minute_string(hour.begin() + 3, hour.begin() + 5),
              second_string(hour.begin() + 6, hour.begin() + 8);

          const int hour_number = std::stoi(hour_string),
                    minute_number = std::stoi(minute_string),
                    second_number = std::stoi(second_string);

          const double hour_full = double(hour_number) +
                                   double(minute_number) / 60 +
                                   double(second_number) / 3600;

          double value;
          ff >> value;

          const int rr = int(std::round(24. / time_spacing));
          if (std::abs((i % rr) * time_spacing - hour_full) >
              .1 * (i % rr) * time_spacing) {
            std::cout << str1 << " " << hour << " " << value << std::endl;
            std::cout << "Invalid rain file" << std::endl;
            exit(1.);
          }
          // mm/h  --> m/sec.
          Hyetograph[0].push_back(value * 1.e-3 / (time_spacing * 3600));
        }
        ff.close();
      } else {
        std::cout << "Unable to open the file, check rain_file in the "
                     "SMARTSED_input file"
                  << std::endl;
        exit(-1.);
      }
    } else {
      for (unsigned int i = 0; i < ndata; i++) {
        Hyetograph[0].push_back(0.);
      }
    }
    for (unsigned int i = 0; i < IDW_weights.size(); i++) {
      IDW_weights[i].push_back(1.);
    }
  }

  void IDW_precipitation(const std::vector<std::string> &file_vect,
                         const std::vector<unsigned int> &ndata_vec,
                         const std::vector<double> &time_spacing_vect,
                         const std::vector<double> &X,
                         const std::vector<double> &Y, const double xllcorner,
                         const double yllcorner, const double pixel_size,
                         const unsigned int N_rows, const unsigned int N_cols,
                         const U_type &idBasinVect) {
    M_time_spacing_vect = time_spacing_vect;

    Hyetograph.resize(ndata_vec.size());

    for (unsigned int k = 0; k < file_vect.size(); k++) {
      const auto file = file_vect[k];
      const auto time_spacing = time_spacing_vect[k];
      std::ifstream ff(file);

      if (ff.is_open()) {

        std::string str;
        std::getline(ff, str);

        for (unsigned int i = 0; i < ndata_vec[k]; i++) {
          std::string Id_sensor;
          ff >> Id_sensor;

          std::string day;
          ff >> day;

          std::string hour;
          ff >> hour;

          std::string hour_string(hour.begin(), hour.begin() + 2),
              minute_string(hour.begin() + 3, hour.begin() + 5);

          const int hour_number = std::stoi(hour_string),
                    minute_number = std::stoi(minute_string);

          const double hour_full =
              double(hour_number) + double(minute_number) / 60;

          double value;
          ff >> value;

          const int rr = int(std::round(24. / time_spacing));
          if (std::abs((i % rr) * time_spacing - hour_full) >
              .1 * (i % rr) * time_spacing) {
            std::cout << hour << std::endl;
            std::cout << "Invalid rain file" << std::endl;
            exit(1.);
          }

          // mm/h  --> m/sec.
          Hyetograph[k].push_back(value * 1.e-3 / (time_spacing * 3600));
        }

        ff.close();

      } else {
        std::cout << "Unable to open the " << k + 1
                  << "-th rain file, check the SMARTSED_input file" << std::endl;
        exit(-1.);
      }
    }

    // fill NO_DATA values
    unsigned int k = 0;
    for (auto &it : Hyetograph) {

      for (unsigned int i = 0; i < ndata_vec[k]; i++) {
        auto &value = it[i];

        if (value < 0) // -999.0 in ARPA files
        {
          std::vector<double> otherStationsRain;
          unsigned int kk = 0;
          for (const auto &itt : Hyetograph) {

            const int ii =
                std::floor(i * (time_spacing_vect[k] / time_spacing_vect[kk]));
            if (itt[ii] >= 0.) {
              otherStationsRain.push_back(itt[ii]);
            }
            kk++;
          }

          if (otherStationsRain.size() != 0) {
            double sum = 0;
            for (const auto &iter : otherStationsRain) {
              sum += iter;
            }
            value = sum / otherStationsRain.size();
          } else {
            value = 0.;
          }
        }
      }
      k++;
    }

    for (const auto &it : Hyetograph) {
      for (const auto &value : it) {
        if (value < 0.) {
          std::cout << value << " One negative value in class rain!" << std::endl;
          exit(1.);
        }
      }
    }

    // compute distances for IDW method
    for (const auto &IDcenter : idBasinVect) {
      const int i = IDcenter / N_cols, j = IDcenter % N_cols;
 
      const double X_cell = j * pixel_size + xllcorner,
                   Y_cell = -i * pixel_size + yllcorner + N_rows * pixel_size;

      for (unsigned int ii = 0; ii < X.size(); ii++) {
        const double p = 3.;
        // divide by 1000 to obtain better number, the result does not change
        const auto delta_x = (X_cell - X[ii]) / 1000,
                   delta_y = (Y_cell - Y[ii]) / 1000;
        const auto dist =
            std::pow(std::sqrt(std::pow(delta_x, 2) + std::pow(delta_y, 2)), p);
        if (dist == 0) {
          IDW_weights[IDcenter].push_back(1.);
        } else {
          IDW_weights[IDcenter].push_back(1. / dist);
        }
      }
      double sum = 0;
      for (const auto &iter : IDW_weights[IDcenter]) {
        sum += iter;
      }
      for (unsigned int ii = 0; ii < IDW_weights[IDcenter].size(); ii++) {
        IDW_weights[IDcenter][ii] /= sum;
      }
    }
  }

  template <class T_type, class Ud_type>
  void computePrecipitation(const double time, const T_type &S,
                            const T_type &melt_mask,
                            const T_type &h_G,
                            const T_type &H,
                            const unsigned int N_rows,
                            const unsigned int N_cols,
                            const Ud_type &idBasinVect) {
#ifdef ENABLE_CUDA
    computePrecipitationKernel(time, );
#else
    // SCS-CN method and Initial and Constant Loss Model
    for (unsigned int Id = 0; Id < Hyetograph.size(); Id++) {

      const unsigned int i_index =
          std::floor(time / (M_time_spacing_vect[Id] * 3600));

      for (const auto &IDcenter : idBasinVect) {

        if (Id == 0) {
          DP_total[IDcenter] = 0.;
          DP_cumulative[IDcenter] = 0.;
          DP_infiltrated[IDcenter] = 0.;
        }

        rainfall_intensity = Hyetograph[Id][i_index] * IDW_weights[IDcenter][Id];

        const double deltaSoilMoisture = h_G[IDcenter] - S[IDcenter];

        double weight = 0.;
        if (S[IDcenter] > 0 && deltaSoilMoisture < 0.) {
          weight = std::pow(deltaSoilMoisture / S[IDcenter], 2.);
        }

        if (weight > 1. || h_G[IDcenter] < 0) {
          std::cout << "Error in weight infiltration model\n"
                    << "h_G = " << h_G[IDcenter] << "\nweight = " << weight
                    << "\nmax soil moisture ret. = " << S[IDcenter] << std::endl;
          exit(-1);
        }

        double infiltrationRate =
                   weight * rainfall_intensity * melt_mask[IDcenter],
               potential_runoff = std::max(
                   rainfall_intensity * melt_mask[IDcenter] - infiltrationRate,
                   0.);

        if (((H[IDcenter] + h_G[IDcenter]) < (c * S[IDcenter])) &&
            M_isInitialLoss) // initial loss
        {
          potential_runoff = 0.;
          infiltrationRate = rainfall_intensity * melt_mask[IDcenter];
        }

        DP_total[IDcenter] += rainfall_intensity;
        DP_cumulative[IDcenter] += potential_runoff;
        DP_infiltrated[IDcenter] += infiltrationRate;
      }
    }
#endif
  }

  T_type DP_total, DP_cumulative, DP_infiltrated;

  double dt_rain;

private:
  std::vector<std::vector<double>> Hyetograph, // # station times ndata
      IDW_weights;

  T_type M_time_spacing_vect;
  bool M_isInitialLoss;
  double rainfall_intensity = 0;
  double c;
};

//==============================================================================

template <class T_type, class U_type>
class Temperature {
public:
  Temperature(const std::string &file, const unsigned int N,
              const unsigned int max_Days, const double T_crit,
              const unsigned int ndata,
              const unsigned int steps_per_hour, const double time_spacing,
              const double height_thermometer, const std::string& format_temp, 
	      const T_type &orography, const U_type &idBasinVect) 
    : T_crit(T_crit), height_th(height_thermometer), 
	orography(orography), idBasinVect(idBasinVect) {

    T_raster.resize(N);
    melt_mask.resize(N);

#ifdef ENABLE_CUDA
    thrust::host_vector<double> T_dailyMean_pot, T_dailyMin_pot, T_dailyMax_pot;
#else
    std::vector<double> T_dailyMean_pot, T_dailyMin_pot, T_dailyMax_pot;
#endif
    
    T_dailyMean_pot.resize(max_Days);
    T_dailyMin_pot.resize(max_Days);
    T_dailyMax_pot.resize(max_Days);
    J.resize(max_Days);

    Temperature_Graph.reserve(ndata);

    std::vector<double> J_ndata;
    J_ndata.reserve(ndata);

    std::ifstream ff(file);

    if (ff.is_open()) {

      if (format_temp == "comune") {
        std::string str;
        std::getline(ff, str);

        for (unsigned int i = 0; i < ndata; i++) {

          std::string str1;
          ff >> str1;

          std::vector<unsigned int> nn;

          if (str1.length() != 10) {
            std::cout << str1 << std::endl;
            std::cout << "Wrong Temperature file format" << std::endl;
            exit(1.);
          }

          std::string day_string(str1.begin(), str1.begin() + 2),
              month_string(str1.begin() + 3, str1.begin() + 5),
              year_string(str1.begin() + 6, str1.end());

          const int day = std::stoi(day_string), month = std::stoi(month_string),
                    year = std::stoi(year_string);

          for (unsigned int uu = 1; uu < month; uu++) {
            if (uu == 4 || uu == 6 || uu == 9 || uu == 11) {
              nn.push_back(30);
            } else if (uu != 2) {
              nn.push_back(31);
            } else if ((year % 4 == 0 && year % 100 != 0) ||
                       year % 400 == 0) // febbraio bisestile
            {
              nn.push_back(29);
            } else // febbraio non bisestile
            {
              nn.push_back(28);
            }
          }

          double scalar_result = 0;
          for (unsigned int uu = 0; uu < nn.size(); uu++) {
            scalar_result += nn[uu];
          }

          J_ndata.push_back(day + scalar_result);

          std::string hour;
          ff >> hour;

          std::string hour_string(hour.begin(), hour.begin() + 2),
              minute_string(hour.begin() + 3, hour.begin() + 5),
              second_string(hour.begin() + 6, hour.begin() + 8);

          const int hour_number = std::stoi(hour_string),
                    minute_number = std::stoi(minute_string),
                    second_number = std::stoi(second_string);

          const double hour_full = double(hour_number) +
                                   double(minute_number) / 60 +
                                   double(second_number) / 3600;

          double value; // temperature data
          ff >> value;

          const int rr = int(std::round(24. / time_spacing));
          if (std::abs((i % rr) * time_spacing - hour_full) >
              .1 * (i % rr) *
                  time_spacing) 
                      
          {
            std::cout << str1 << " " << hour << " " << value << " " << i + 1
                      << std::endl;
            std::cout << "Invalid temperature file, maybe check time spacing in "
                         "SMARTSED_input file"
                      << std::endl;
            exit(1.);
          }

          if (value == -999 && i == 0) {
            std::cout << "Correct temperature file to eliminate first element as "
                         "NODATA value"
                      << std::endl;
            exit(-1.);
          }

          if (value == -999) {
            value = Temperature_Graph[i - 1];
          }

          Temperature_Graph.push_back(value);
        }

        ff.close();
      } else if (format_temp == "arpa") {

        std::string str;
        std::getline(ff, str);

        for (unsigned int i = 0; i < ndata; i++) {
          unsigned int Id_sensor;
          ff >> Id_sensor;

          std::string str1;
          ff >> str1;

          std::vector<unsigned int> nn;

          if (str1.length() != 10) {
            std::cout << str1 << std::endl;
            std::cout << "Wrong Temperature file format" << std::endl;
            exit(1.);
          }

          std::string year_string(str1.begin(), str1.begin() + 4),
              month_string(str1.begin() + 5, str1.begin() + 7),
              day_string(str1.begin() + 8, str1.end());

          const int day = std::stoi(day_string), month = std::stoi(month_string),
                    year = std::stoi(year_string);

          for (unsigned int uu = 1; uu < month; uu++) {
            if (uu == 4 || uu == 6 || uu == 9 || uu == 11) {
              nn.push_back(30);
            } else if (uu != 2) {
              nn.push_back(31);
            } else if ((year % 4 == 0 && year % 100 != 0) ||
                       year % 400 == 0) // febbraio bisestile
            {
              nn.push_back(29);
            } else // febbraio non bisestile
            {
              nn.push_back(28);
            }
          }

          double scalar_result = 0;
          for (unsigned int uu = 0; uu < nn.size(); uu++) {
            scalar_result += nn[uu];
          }

          J_ndata.push_back(day + scalar_result);

          std::string hour;
          ff >> hour;

          std::string hour_string(hour.begin(), hour.begin() + 2),
              minute_string(hour.begin() + 3, hour.begin() + 5);

          const int hour_number = std::stoi(hour_string),
                    minute_number = std::stoi(minute_string);

          const double hour_full =
                double(hour_number) + double(minute_number) / 60;

          double value; // temperature data
          ff >> value;

          const int rr = int(std::round(24. / time_spacing));
          if (std::abs((i % rr) * time_spacing - hour_full) >
              .1 * (i % rr) * time_spacing) {
            std::cout << str1 << " " << hour << " " << value << std::endl;
            std::cout << "Invalid temperature file" << std::endl;
            exit(1.);
          }

          if (value == -999 && i == 0) {
            std::cout << "Correct temperature file to eliminate first element as "
                         "NODATA value"
                      << std::endl;
            exit(-1.);
          }

          if (value == -999) {
            value = Temperature_Graph[i - 1];
          }

          Temperature_Graph.push_back(value);
        }

        ff.close();
      } else {
        std::cout << "Temperature format non recognized" << std::endl;
        exit(-1.);
      }

    } else {
      std::cout << "Unable to open the file, check temperature_file in the "
                   "SMARTSED_input file"
                << std::endl;
      exit(-1.);
    }

    // --------------------------------------------- //
    for (unsigned int n = 1; n <= max_Days; n++) {

      const auto i = std::floor((n - 1) * (24. / time_spacing));

      unsigned int k = 0,
                   h = i; 

      while (k != static_cast<unsigned int>(std::round((24. / time_spacing)))) {

        T_dailyMean_pot[n - 1] += Temperature_Graph[h];

        if (k != 0) {

          if (Temperature_Graph[h] < T_dailyMin[n - 1]) {
            T_dailyMin_pot[n - 1] = Temperature_Graph[h];
          }

          if (Temperature_Graph[h] > T_dailyMax[n - 1]) {
            T_dailyMax_pot[n - 1] = Temperature_Graph[h];
          }
        } else {
          T_dailyMin_pot[n - 1] = Temperature_Graph[h];
          T_dailyMax_pot[n - 1] = Temperature_Graph[h];
        }

        h++;
        k++;
      }

      if (k == 0) {
        std::cout << "Something wrong in Temperature class constructor"
                  << std::endl;
        exit(-1.);
      }

      T_dailyMean_pot[n - 1] /= k;
    }
    // --------------------------------------------- //

    // Now fill J starting from J_ndata
    for (unsigned int n = 1; n <= max_Days; n++) {
      const unsigned int i = std::floor((n - 1) * (24. / time_spacing));
      J[n - 1] = J_ndata[i];
    }

#ifdef ENABLE_CUDA
    T_dailyMean = T_dailyMean_pot;
    T_dailyMin = T_dailyMin_pot;
    T_dailyMax = T_dailyMax_pot;
#else
    T_dailyMean = std::move(T_dailyMean_pot);
    T_dailyMin = std::move(T_dailyMin_pot);
    T_dailyMax = std::move(T_dailyMax_pot);
#endif
  }

  Temperature() = delete;
  ~Temperature() = default;

  void computeTemperature(const double time, const double dt_temp, const double dt_DSV) {
#ifdef ENABLE_CUDA
    computeTemperatureKernel();
#else
    // update only if necessary  --> governed by temperature dynamics, i.e.
    // time_spacing_temp
    if (std::floor(time / dt_temp) > std::floor((time - dt_DSV) / dt_temp)) {
      const unsigned int i_t = std::floor(time / dt_temp);
      const auto T = Temperature_Graph[i_t];
      if (std::isnan(T)) {
        std::cout << "NAN in computeTemperature" << std::endl;
        exit(-1);
      }
      for (const auto &j : idBasinVect) {
        T_raster[j] = T + Temp_diff * (orography[j] - height_th);
        melt_mask[j] = (T_raster[j] > T_crit); // melt_mask = 1 -\mu in the paper
      }
    }
#endif
  }

  T_type T_raster, melt_mask, T_dailyMean, T_dailyMin, T_dailyMax;
  std::vector<double> J;

private:
  const T_type&orography;
  const U_type&idBasinVect;
  T_type Temperature_Graph; // length: ndata
  static constexpr double Temp_diff = -6.5e-3;
  const double height_th;
  const double T_crit;
};

//==============================================================================

template <class T_type, class U_type>
class evapoTranspiration {
public:
  evapoTranspiration(const std::string &ET_model, const unsigned int &N,
                     const std::vector<double> &J, const unsigned int &max_Days,
                     const double &phi_rad, const double &height_thermometer,
		     const U_type &idBasinVect,
                     const T_type &orography) 
      : height_th(height_thermometer), idBasinVect(idBasinVect),
        orography(orography) {

#ifdef ENABLE_CUDA
    thrust::host_vector<double> Ra_pot;
#else
    std::vector<double> Ra_pot;
#endif

    ET_vec.resize(N);

    if (ET_model == "None") {
      M_evapoTranspiration_model = 0;
    } else if (ET_model == "Hargreaves") {
      M_evapoTranspiration_model = 1;
    } else {
      std::cout << "No evapo-transpiration model inserted!!" << std::endl;
      exit(-1.);
    }

    Ra_pot.resize(max_Days);
    for (unsigned int n = 1; n <= max_Days; n++) {
      const auto dr = 1 + .033 * std::cos(2 * M_PI * J[n - 1] / 365),
                 delta = .409 * std::sin(2 * M_PI * J[n - 1] / 365 - 1.39),
                 ws = std::acos(-std::tan(phi_rad) * std::tan(delta));
      Ra_pot[n - 1] = (24 * 60 / M_PI) * M_Gsc * dr *
                  (ws * std::sin(phi_rad) * std::sin(delta) +
                   std::cos(phi_rad) * std::cos(delta) * std::sin(ws));
    }

#ifdef ENABLE_CUDA
    Ra = Ra_pot;
#else
    Ra = std::move(Ra_pot);
#endif

  }

  evapoTranspiration() = delete;
  ~evapoTranspiration() = default;

  void ET(const T_type &T_mean, // length nstep: vector of temperature in deg Celsius
          const T_type &T_min, // length nstep
          const T_type &T_max, // length nstep
          const unsigned int i) {
#ifdef ENABLE_CUDA
    ETKernel();
#else
    switch (M_evapoTranspiration_model) {
    case 0:
      break;
    case 1:
      for (const auto &k : idBasinVect) {
        const auto t_mean = T_mean[i] + Temp_diff * (orography[k] - height_th),
                   t_max  = T_max [i] + Temp_diff * (orography[k] - height_th),
                   t_min  = T_min [i] + Temp_diff * (orography[k] - height_th);
        // unity: mm/day --> m/sec.
        ET_vec[k] = .0023 * Ra[i] * (t_mean + 17.8) *
                    std::pow((t_max - t_min), .5) * (1.e-3 / (24 * 3600));
      }
      break;
    }
#endif
  }

  T_type ET_vec;

private:
  const T_type&orography;
  const U_type&idBasinVect;
  T_type Ra;
  unsigned int M_evapoTranspiration_model;
  static constexpr double M_Gsc = 0.082; // Solar constant
  const double height_th;
  static constexpr double Temp_diff = -6.5e-3;
};

//==============================================================================

template<class T_type, class U_type, class V_type>
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
      const V_type& S_x, const V_type& S_y
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
    compute_horizontal_friction_wrapper(idStaggeredInternalVectHorizontal,
                idStaggeredBoundaryVectWest, idStaggeredBoundaryVectEast,
                H_interface_horizontal, u, M_expo_r_x_vect,
                alfa_x, M_gamma_dt_DSV_x_, M_dt_DSV, M_coeff, M_H_min,
                M_expo, M_frictionModel);
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
    compute_vertical_friction_wrapper(idStaggeredInternalVectVertical,
                idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectSouth,
                H_interface_vertical, v, M_expo_r_y_vect,
                alfa_y, M_gamma_dt_DSV_y_, M_dt_DSV, M_coeff, M_H_min,
                M_expo, M_frictionModel);
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

#ifndef ENABLE_CUDA
  std::function<double(double const &, double const &)> M_gamma_dt_DSV =
      [](double const &dt, double const &cc) { return dt * cc; };
#endif

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
    compute_horizontal_interface_wrapper(idStaggeredInternalVectHorizontal,
                idStaggeredBoundaryVectWest, idStaggeredBoundaryVectEast,
                H, u, horizontal, N_cols);
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
    compute_vertical_interface_wrapper(idStaggeredInternalVectVertical,
                idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectSouth,
                H, v, vertical, N_cols);
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

    const std::vector<std::tuple<bool, int>> &excluded_ids,
    std::vector<double> &additional_source_term,

    std::vector<Eigen::Triplet<double>> &coefficients, Eigen::VectorXd &rhs);

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
    const bool &isNonReflectingBC);

//==============================================================================

void compute_dt_adaptive(const std::vector<double> &H,
                         const std::vector<double> &H_old,
                         const std::vector<double> &H_oldold,
                         const std::vector<unsigned int> &idBasinVect,
                         double &dt,
                         const double &local_estimator_time_tolerance,
                         const double &time, const double &timed,
                         const double &timedd);

//==============================================================================

double maxdt(const std::vector<double> &u, const std::vector<double> &v,
             const double &Hmax, const double &pixel_size);

//==============================================================================

double maxCourant(const std::vector<double> &u, const std::vector<double> &v,
                  const double &c1);

//==============================================================================

double maxCourant(const std::vector<double> &H, const double &c1);

//==============================================================================

double compute_dt_sediment(const double alpha, const double beta,
                           const double S_x, const double S_y,
                           const std::vector<double> &u,
                           const std::vector<double> &v,
                           const double pixel_size, const double dt_DSV,
                           unsigned int *numberOfSteps);

//==============================================================================

int current_start_chunk(const int &rank,
                        const std::vector<int> &chunk_length_vec);

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
    std::vector<double> &Res_x, std::vector<double> &Res_y);

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
    std::vector<std::array<double, 2>> &Gamma_y);

//==============================================================================
