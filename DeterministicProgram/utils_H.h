#ifndef UTILS_H_H
#define UTILS_H_H

//! std library
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <complex>
#include <numeric>

//! pure CPU header
#include "code_init.h"

//! Eigen library
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/SparseCholesky>

//! Include GPU impl.
#ifdef ENABLE_CUDA
#include "cuda_utils_loop_H.cuh"
#include <thrust/extrema.h>
#endif

//! Parse library
#include "GetPot.hpp"

//! IML++ CG template
#include "cg.hpp"

//==============================================================================

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
  using T_type = thrust::host_vector<double>;
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

    DP_total_pot.resize(N);
    DP_cumulative_pot.resize(N);
    DP_infiltrated_pot.resize(N);
    IDW_weights.resize(N);
   
#ifdef ENABLE_CUDA
    DP_total = DP_total_pot;
    DP_cumulative = DP_cumulative_pot;
    DP_infiltrated = DP_infiltrated_pot;
#else
    DP_total = std::move(DP_total_pot);
    DP_cumulative = std::move(DP_cumulative_pot);
    DP_infiltrated = std::move(DP_infiltrated_pot);
#endif

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

#ifdef ENABLE_CUDA

    thrust::host_vector<unsigned int> size_hy(Hyetograph.size()), 
	    size_IDW(IDW_weights.size());
    unsigned int size_hy_s = 0, size_IDW_s = 0;
    for (int i = 0; i < Hyetograph.size(); i++) {
	    size_hy[i] = Hyetograph[i].size();
	    size_hy_s += Hyetograph[i].size();
    }

    for (int i = 0; i < IDW_weights.size(); i++) {
	    size_IDW[i] = IDW_weights[i].size();
	    size_IDW_s += IDW_weights[i].size();
    }

    thrust::host_vector<double> flat_hy, flat_IDW;

    flat_hy.resize(size_hy_s);
    flat_IDW.resize(size_IDW_s);

    for (int i = 0; i < Hyetograph.size(); i++)
	    for (int j = 0; j < Hyetograph[i].size(); j++)
		    flat_hy[std::accumulate(size_hy.begin(), size_hy.begin()+i, 0) + j] 
			    = Hyetograph[i][j];

     for (int i = 0; i < IDW_weights.size(); i++)
	    for (int j = 0; j < IDW_weights[i].size(); j++)
		    flat_IDW[std::accumulate(size_IDW.begin(), size_IDW.begin()+i, 0) + j] 
			    = IDW_weights[i][j];

     Hyetograph_gpu = flat_hy;
     IDW_weights_gpu = flat_IDW;

     M_time_spacing_vect_gpu = M_time_spacing_vect;

     thrust::host_vector<unsigned int> offset_hy(Hyetograph.size() + 1, 0);
     thrust::inclusive_scan(size_hy.begin(), size_hy.end(), offset_hy.begin() + 1);

     d_offset_hy = offset_hy;

#endif

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

#ifdef ENABLE_CUDA

    thrust::host_vector<unsigned int> size_hy(Hyetograph.size()), 
	    size_IDW(IDW_weights.size());
    unsigned int size_hy_s = 0, size_IDW_s = 0;
    for (int i = 0; i < Hyetograph.size(); i++) {
	    size_hy[i] = Hyetograph[i].size();
	    size_hy_s += Hyetograph[i].size();
    }

    for (int i = 0; i < IDW_weights.size(); i++) {
	    size_IDW[i] = IDW_weights[i].size();
	    size_IDW_s += IDW_weights[i].size();
    }

    thrust::host_vector<double> flat_hy, flat_IDW;

    flat_hy.resize(size_hy_s);
    flat_IDW.resize(size_IDW_s);

    for (int i = 0; i < Hyetograph.size(); i++)
	    for (int j = 0; j < Hyetograph[i].size(); j++)
		    flat_hy[std::accumulate(size_hy.begin(), size_hy.begin()+i, 0) + j] 
			    = Hyetograph[i][j];

     for (int i = 0; i < IDW_weights.size(); i++)
	    for (int j = 0; j < IDW_weights[i].size(); j++)
		    flat_IDW[std::accumulate(size_IDW.begin(), size_IDW.begin()+i, 0) + j] 
			    = IDW_weights[i][j];

     Hyetograph_gpu = flat_hy;
     IDW_weights_gpu = flat_IDW;

     M_time_spacing_vect_gpu = M_time_spacing_vect;

     thrust::host_vector<unsigned int> offset_hy(Hyetograph.size() + 1, 0);
     thrust::inclusive_scan(size_hy.begin(), size_hy.end(), offset_hy.begin() + 1);

     d_offset_hy = offset_hy;

#endif

  }

  template <class T_type, class Ud_type>
  void computePrecipitation(const double time, const T_type &S,
                            const T_type &melt_mask,
                            const T_type &h_G,
                            const T_type &H,
                            const unsigned int N_cols,
                            const Ud_type &idBasinVect
#ifdef ENABLE_CUDA
			    , cudaStream_t stream = 0
#endif

		  ) {
#ifdef ENABLE_CUDA

      computePrecipitation_wrapper(
        idBasinVect,
        time, c, M_isInitialLoss,
	M_time_spacing_vect_gpu,
	S,
        melt_mask,
	Hyetograph_gpu,
	IDW_weights_gpu,
	DP_total, DP_cumulative, DP_infiltrated,
	h_G, H, N_cols, d_offset_hy, Hyetograph.size(), stream);

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

        const double rainfall_intensity = Hyetograph[Id][i_index] * IDW_weights[IDcenter][Id];

        const double deltaSoilMoisture = h_G[IDcenter] - S[IDcenter];

        double weight = 0.;
        if (S[IDcenter] > 0 && deltaSoilMoisture < 0.) {
          weight = (deltaSoilMoisture / S[IDcenter])*(deltaSoilMoisture / S[IDcenter]);
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

  T_type DP_total_pot, DP_cumulative_pot, DP_infiltrated_pot;

#ifdef ENABLE_CUDA
  thrust::device_vector<unsigned int> d_offset_hy;
  thrust::device_vector<double> DP_total, DP_cumulative, DP_infiltrated,
	  M_time_spacing_vect_gpu;
#else
  T_type DP_total, DP_cumulative, DP_infiltrated;
#endif

  double dt_rain;

private:
#ifdef ENABLE_CUDA
  thrust::device_vector<double> Hyetograph_gpu, IDW_weights_gpu;
  thrust::host_vector<thrust::host_vector<double>> Hyetograph, // # station times ndata
	  IDW_weights;
#else
  std::vector<std::vector<double>> Hyetograph, // # station times ndata
	  IDW_weights;
#endif

  T_type M_time_spacing_vect;
  bool M_isInitialLoss;
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
    
    T_dailyMean.resize(max_Days);
    T_dailyMin.resize(max_Days);
    T_dailyMax.resize(max_Days);
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
                       year % 400 == 0) // leap year February
            {
              nn.push_back(29);
            } else // non-leap year February
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
                       year % 400 == 0) // leap year February
            {
              nn.push_back(29);
            } else // non-leap year February
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

        T_dailyMean[n - 1] += Temperature_Graph[h];

        if (k != 0) {

          if (Temperature_Graph[h] < T_dailyMin[n - 1]) {
            T_dailyMin[n - 1] = Temperature_Graph[h];
          }

          if (Temperature_Graph[h] > T_dailyMax[n - 1]) {
            T_dailyMax[n - 1] = Temperature_Graph[h];
          }
        } else {
          T_dailyMin[n - 1] = Temperature_Graph[h];
          T_dailyMax[n - 1] = Temperature_Graph[h];
        }

        h++;
        k++;
      }

      if (k == 0) {
        std::cout << "Something wrong in Temperature class constructor"
                  << std::endl;
        exit(-1.);
      }

      T_dailyMean[n - 1] /= k;
    }
    // --------------------------------------------- //

    // Now fill J starting from J_ndata
    for (unsigned int n = 1; n <= max_Days; n++) {
      const unsigned int i = std::floor((n - 1) * (24. / time_spacing));
      J[n - 1] = J_ndata[i];
    }
  }

  Temperature() = delete;
  ~Temperature() = default;

  void computeTemperature(const double time, const double dt_temp, const double dt_DSV
#ifdef ENABLE_CUDA
         , cudaStream_t stream = 0
#endif

		  ) {
#ifdef ENABLE_CUDA
    if (std::floor(time / dt_temp) > std::floor((time - dt_DSV) / dt_temp)) {
      const unsigned int i_t = std::floor(time / dt_temp);
      const double T = Temperature_Graph[i_t];
      if (std::isnan(T)) {
        std::cout << "NAN in computeTemperature" << std::endl;
        exit(-1);
      }
      computeTemperature_wrapper(
        idBasinVect,
        T_raster,
        melt_mask,
        orography,
        T, Temp_diff, height_th, T_crit, stream);
    }
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

  T_type T_raster, melt_mask; 
  std::vector<double> J, T_dailyMean, T_dailyMin, T_dailyMax;

private:
  const T_type&orography;
  const U_type&idBasinVect;
  std::vector<double> Temperature_Graph; // length: ndata
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

    ET_vec.resize(N);

    if (ET_model == "None") {
      M_evapoTranspiration_model = 0;
    } else if (ET_model == "Hargreaves") {
      M_evapoTranspiration_model = 1;
    } else {
      std::cout << "No evapo-transpiration model inserted!!" << std::endl;
      exit(-1.);
    }

    Ra.resize(max_Days);
    for (unsigned int n = 1; n <= max_Days; n++) {
      const auto dr = 1 + .033 * std::cos(2 * M_PI * J[n - 1] / 365),
                 delta = .409 * std::sin(2 * M_PI * J[n - 1] / 365 - 1.39),
                 ws = std::acos(-std::tan(phi_rad) * std::tan(delta));
      Ra[n - 1] = (24 * 60 / M_PI) * M_Gsc * dr *
                  (ws * std::sin(phi_rad) * std::sin(delta) +
                   std::cos(phi_rad) * std::cos(delta) * std::sin(ws));
    }
  }

  evapoTranspiration() = delete;
  ~evapoTranspiration() = default;

  void ET(const T_type &T_mean, // length nstep: vector of temperature in deg Celsius
          const T_type &T_min, // length nstep
          const T_type &T_max, // length nstep
          const unsigned int i
#ifdef ENABLE_CUDA
	  , cudaStream_t stream = 0
#endif

	  ) {
#ifdef ENABLE_CUDA
    if (M_evapoTranspiration_model == 0) return;  // early exit for case 0

    computeET_wrapper(
      idBasinVect, ET_vec, orography, 
      Ra[i],
      T_mean[i],   // scalar base temps — host access
      T_max[i],
      T_min[i],
      Temp_diff, height_th,
      M_evapoTranspiration_model,
      stream);
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
      const double &dt_DSV, 
#ifndef ENABLE_CUDA
      const std::vector<double> &d_90,
      const std::vector<double> &rough, 
#else
      const thrust::host_vector<double> &d_90,
      const thrust::host_vector<double> &rough,
#endif
      const double &H_min,
      const unsigned int &N_rows, const unsigned int &N_cols, 
      const T_type& S_x, const T_type& S_y
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

#ifdef ENABLE_CUDA // copy to GPU the vectors
    M_expo_r_x_vect_gpu = M_expo_r_x_vect;
    M_expo_r_y_vect_gpu = M_expo_r_y_vect; 
    M_gamma_dt_DSV_x_gpu_ = M_gamma_dt_DSV_x_;
    M_gamma_dt_DSV_y_gpu_ = M_gamma_dt_DSV_y_;
#endif
  }

  frictionClass() = delete;
  ~frictionClass() = default;

  void f_x(
#ifdef ENABLE_CUDA
    cudaStream_t stream = 0
#endif
		  ) {
#ifdef ENABLE_CUDA
    compute_horizontal_friction_wrapper(idStaggeredInternalVectHorizontal,
                idStaggeredBoundaryVectWest, idStaggeredBoundaryVectEast,
                H_interface_horizontal, u, M_expo_r_x_vect_gpu,
                alfa_x, M_gamma_dt_DSV_x_gpu_, M_dt_DSV, M_coeff, M_H_min,
                M_expo, M_frictionModel, stream);
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


  void f_y(
#ifdef ENABLE_CUDA
    cudaStream_t stream = 0
#endif
		  ) {
#ifdef ENABLE_CUDA
    compute_vertical_friction_wrapper(idStaggeredInternalVectVertical,
                idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectSouth,
                H_interface_vertical, v, M_expo_r_y_vect_gpu,
                alfa_y, M_gamma_dt_DSV_y_gpu_, M_dt_DSV, M_coeff, M_H_min,
                M_expo, M_frictionModel, stream);
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

#ifndef ENABLE_CUDA
  std::vector<double> M_expo_r_x_vect, M_expo_r_y_vect, M_gamma_dt_DSV_x_,
      M_gamma_dt_DSV_y_;
#else
  thrust::host_vector<double> M_expo_r_x_vect, M_expo_r_y_vect, M_gamma_dt_DSV_x_,
      M_gamma_dt_DSV_y_;
  thrust::device_vector<double> M_expo_r_x_vect_gpu, M_expo_r_y_vect_gpu, M_gamma_dt_DSV_x_gpu_,
      M_gamma_dt_DSV_y_gpu_;
#endif
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

  void computeHorizontal(
#ifdef ENABLE_CUDA
    cudaStream_t stream = 0  // default stream if not specified
#endif
		  ) {
#ifdef ENABLE_CUDA
    compute_horizontal_interface_wrapper(idStaggeredInternalVectHorizontal,
                idStaggeredBoundaryVectWest, idStaggeredBoundaryVectEast,
                H, u, horizontal, N_cols, stream);
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

  void computeVertical(
#ifdef ENABLE_CUDA
    cudaStream_t stream = 0  // default stream if not specified
#endif
		  ) {
#ifdef ENABLE_CUDA
    compute_vertical_interface_wrapper(idStaggeredInternalVectVertical,
                idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectSouth,
                H, v, vertical, N_cols, stream);
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

template<class T_type, class U_type>
void bilinearInterpolation(
    const T_type &u, const T_type &v,
    T_type &u_star, T_type &v_star,
    const unsigned int nrows, const unsigned int ncols, const double dt_DSV,
    const double pixel_size,
    const U_type &idStaggeredInternalVectHorizontal,
    const U_type &idStaggeredInternalVectVertical,
    const U_type &idStaggeredBoundaryVectWest,
    const U_type &idStaggeredBoundaryVectEast,
    const U_type &idStaggeredBoundaryVectNorth,
    const U_type &idStaggeredBoundaryVectSouth
#ifdef ENABLE_CUDA
    , cudaStream_t stream = 0  // default stream if not specified
#endif
    ) {

#ifdef ENABLE_CUDA

  bilinearInterpolationHorizontal_wrapper(idStaggeredInternalVectHorizontal,
                idStaggeredBoundaryVectWest, idStaggeredBoundaryVectEast,
                u, v, u_star, dt_DSV/pixel_size, nrows, ncols, stream);

  bilinearInterpolationVertical_wrapper(idStaggeredInternalVectVertical,
                idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectSouth,
                u, v, v_star, dt_DSV/pixel_size, nrows, ncols, stream);

#else

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
      y_1 += 1;
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
#endif
}

//==============================================================================

template <class T_type, class U_type>
void buildMatrix(
    const T_type &H_int_x, const T_type &H_int_y,
    const T_type &orography, const T_type &u_star,
    const T_type &v_star, const T_type &u,
    const T_type &v, const T_type &H,
    const unsigned int N_cols,
    const double c1, const double c3,
    const double H_min, const T_type &precipitation,
    const double dt_DSV, const T_type &alfa_x,
    const T_type &alfa_y,
    const U_type &idStaggeredInternalVectHorizontal,
    const U_type &idStaggeredInternalVectVertical,
    const U_type &idStaggeredBoundaryVectWest,
    const U_type &idStaggeredBoundaryVectEast,
    const U_type &idStaggeredBoundaryVectNorth,
    const U_type &idStaggeredBoundaryVectSouth,
    const U_type &idBasinVect,
#ifndef ENABLE_CUDA
    const U_type &idBasinVect_not_excluded,
    const U_type
        &idStaggeredInternalVectHorizontal_not_excluded,
    const U_type
        &idStaggeredInternalVectVertical_not_excluded,
#endif
    const U_type &idBasinVectReIndex,
    const bool isNonReflectingBC, const bool isH,

#ifdef ENABLE_CUDA
    const int *d_A_rows, const int *d_A_columns, double *d_A_values, Vec &rhs, const int nnz,
    cudaStream_t stream = 0  // default stream if not specified
#else
    T_type &additional_source_term, const std::vector<std::tuple<bool, int>> &excluded_ids,
    std::vector<Eigen::Triplet<double>> &coefficients, Eigen::VectorXd &rhs
#endif
    ) {

#ifdef ENABLE_CUDA

  buildMatrix_wrapper(H_int_x, H_int_y, orography, 
		  u_star, v_star, u, v, H, N_cols, c1, c3, H_min, precipitation,
		  dt_DSV, alfa_x, alfa_y, idStaggeredInternalVectHorizontal, 
		  idStaggeredInternalVectVertical, idStaggeredBoundaryVectWest, 
		  idStaggeredBoundaryVectEast, idStaggeredBoundaryVectNorth, 
		  idStaggeredBoundaryVectSouth, idBasinVect, idBasinVectReIndex,
		  isNonReflectingBC, nnz, d_A_rows, d_A_columns, d_A_values, rhs, stream);

#else
	    
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
#endif
}

//==============================================================================

template <class T_type, class U_type>
void updateVel(
    T_type &u, T_type &v,
    const T_type &u_star, const T_type &v_star,
    const T_type &alfa_x, const T_type &alfa_y,
    const double N_rows, const double N_cols, const double c2,
    const double H_min, const T_type &eta,
    const T_type &H, const T_type &orography,
    const U_type &idStaggeredInternalVectHorizontal,
    const U_type &idStaggeredInternalVectVertical,
    const U_type &idStaggeredBoundaryVectWest,
    const U_type &idStaggeredBoundaryVectEast,
    const U_type &idStaggeredBoundaryVectNorth,
    const U_type &idStaggeredBoundaryVectSouth,
    const bool isNonReflectingBC
#ifdef ENABLE_CUDA
    , cudaStream_t stream = 0  // default stream if not specified
#endif
    ) {

#ifdef ENABLE_CUDA

#else
  // +-----------------------------------------------+
  // |           Update Vertical Velocity            |
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
  // |           Update Horizontal Velocity          |
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
#endif
}

//==============================================================================

template <class T_type, class U_type>
void compute_dt_adaptive(const T_type &H,
                         const T_type &H_old,
                         const T_type &H_oldold,
                         const U_type &idBasinVect,
                         double& dt,
                         const double local_estimator_time_tolerance,
                         const double time, const double timed,
                         const double timedd) {

  // capture all constants for the functor
  const double time_ = time, timed_ = timed, timedd_ = timedd;	

  // get raw pointers before the lambda
  const double* H_ptr      = thrust::raw_pointer_cast(H.data());
  const double* H_old_ptr  = thrust::raw_pointer_cast(H_old.data());
  const double* H_oldd_ptr = thrust::raw_pointer_cast(H_oldold.data());

  // functor: given an index, compute the Nu_hmean_cell contribution
  auto compute_nu = [=] __host__ __device__ (unsigned int Id) -> double {

	  const double hcell      = H_ptr[Id];
	  const double hcell_old  = H_old_ptr[Id];
	  const double hcell_oldd = H_oldd_ptr[Id];

	  const double dh_t = (hcell - hcell_old) / (time_ - timed_);

	  const double h1 = hcell_oldd / ((timedd_ - timed_) * (timedd_ - time_));
	  const double h2 = hcell_old  / ((timed_  - timedd_) * (timed_  - time_));
	  const double h3 = hcell      / ((time_   - timedd_) * (time_   - timed_));

	  const double a_coeff = h1 + h2 + h3;
	  const double b_coeff = -(h1 * (time_ + timed_) +
			  h2 * (time_ + timedd_) +
			  h3 * (timed_ + timedd_));

	  const double Nu_hmean_cell =
		  (1. / 3. * a_coeff * a_coeff *
		   (time_ * time_ + time_ * timed_ + timed_ * timed_) +
		   a_coeff * (b_coeff - dh_t) * (time_ + timed_) +
		   (b_coeff - dh_t) * (b_coeff - dh_t));

	  return Nu_hmean_cell * (time_ - timed_) * (time_ - timed_);
  };

  // reduce over idBasinVect
  double nu_htot = thrust::transform_reduce(
		  idBasinVect.begin(),
		  idBasinVect.end(),
		  compute_nu,
		  0.0,
		  thrust::plus<double>());


  double dt_candidate =
      local_estimator_time_tolerance / std::sqrt(nu_htot) * (time - timed);

  // compute new dt
  dt = (nu_htot > 0 && dt_candidate < dt) ? dt_candidate : dt;

}

//==============================================================================

template <class T>
static inline T maxdt_compute(const T vel_max_x, const T vel_max_y,
                              const T Hmax,      const T pixel_size,
			      const T gravity) {

  static_assert(std::is_same<T, double>::value, 
		  "maxdt_compute only supports double");

  const T Co     = 2.;
  const T Co_cel = 1e4;
  const T cel    = std::sqrt(Hmax * gravity);

  T dt = Co * pixel_size /
              (std::max(vel_max_x, vel_max_y) + std::numeric_limits<T>::epsilon());
  dt = std::min(dt, Co_cel * pixel_size /
                    (cel + std::numeric_limits<T>::epsilon()));

  return dt;
}

//==============================================================================

template <class T_type>
double maxdt(const T_type &u, const T_type &v,
             const double Hmax, const double pixel_size,
             const double gravity) {

#ifdef ENABLE_CUDA

  const auto v_mm = deviceMinMax(v);
  const auto u_mm = deviceMinMax(u);

  const double vel_max_y = std::max(v_mm.max_val, std::abs(v_mm.min_val));
  const double vel_max_x = std::max(u_mm.max_val, std::abs(u_mm.min_val));

#else
 
  const double vel_max_y = std::max(*std::max_element(v.begin(), v.end()),
                                    std::abs(*std::min_element(v.begin(), v.end())));
  const double vel_max_x = std::max(*std::max_element(u.begin(), u.end()),
                                    std::abs(*std::min_element(u.begin(), u.end())));

#endif

  return maxdt_compute(vel_max_x, vel_max_y, Hmax, pixel_size, gravity);

}

//==============================================================================

template <class T_type>
double maxCourant(const T_type &u,
                  const T_type &v, const double c1) {
#ifdef ENABLE_CUDA
  const auto v_mm = deviceMinMax(v);
  const auto u_mm = deviceMinMax(u);

  const double Courant_y = std::max(v_mm.max_val, std::abs(v_mm.min_val));
  const double Courant_x = std::max(u_mm.max_val, std::abs(u_mm.min_val));

  return std::max(Courant_y, Courant_x) * c1;
#else
  const double Courant_y = std::max(*std::max_element(v.begin(), v.end()),
                                    std::abs(*std::min_element(v.begin(), v.end())));
  const double Courant_x = std::max(*std::max_element(u.begin(), u.end()),
                                    std::abs(*std::min_element(u.begin(), u.end())));

  return std::max(Courant_y, Courant_x) * c1;
#endif
}

template <class T_type>
double maxCourant(const T_type &H, const double c1, 
		const double gravity) {
#ifdef ENABLE_CUDA
  return std::sqrt(deviceMax(H) * gravity) * c1;  // deviceMax from before 
#else
  const double Courant_cel =
      std::sqrt(*std::max_element(H.begin(), H.end()) * gravity);
  return (Courant_cel * c1);
#endif
}

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
template <class T_type, class U_type>
void computeResiduals(
    const T_type &n_x, const T_type &n_y,
    const unsigned int N_cols,
    const T_type &coeff, // hydraulic conductivity
    const U_type &idStaggeredInternalVectHorizontal,
    const U_type &idStaggeredInternalVectVertical,
    const U_type &idStaggeredBoundaryVectWest,
    const U_type &idStaggeredBoundaryVectEast,
    const U_type &idStaggeredBoundaryVectNorth,
    const U_type &idStaggeredBoundaryVectSouth,
    const U_type &idBasinVect,
    const T_type &T_raster, const T_type &melt_mask,
    T_type &h_sn,
    T_type &h, 
    const T_type &ET_vec, 
    const T_type &DP_infiltrated,
    const T_type &DP_total,
    const double c1_min,
    const double dt_min, 
    const double T_thr, 
    T_type &h_interface_x, T_type &h_interface_y
#ifdef ENABLE_CUDA
    , cudaStream_t stream = 0  // default stream if not specified
#endif    
    ) {

#ifdef ENABLE_CUDA

      computeResidualsHorizontal_wrapper(
        idStaggeredInternalVectHorizontal,
	idStaggeredBoundaryVectWest, idStaggeredBoundaryVectEast,
	coeff, n_x, h,
        h_interface_x, N_cols, stream);
 
      computeResidualsVertical_wrapper(
        idStaggeredInternalVectVertical,
	idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectSouth,
	coeff, n_y, h,
        h_interface_y, N_cols, stream);

      updateSnowGravLayers_wrapper(
        idBasinVect, T_raster, melt_mask, h_sn, 
	h, ET_vec, DP_infiltrated, DP_total, c1_min, dt_min, T_thr,
	N_cols, h_interface_x, h_interface_y, stream);

#else

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

  // +-----------------------------------------------+
  // |                   Internal                    |
  // +-----------------------------------------------+

  for (const auto &k : idBasinVect) {
    const unsigned int i = k / N_cols;	
    const double S_coeff = 4.62e-10 * h_sn[k] * (T_raster[k] - T_thr) * melt_mask[k];
    const double Res_x_cell = h_interface_x[k + 1 + i] - h_interface_x[k + i];
    const double Res_y_cell = h_interface_y[k + N_cols] - h_interface_y[k];	  
    h[k] = std::max(h[k] + (S_coeff - ET_vec[k]) * dt_min + DP_infiltrated[k] * dt_min - c1_min * (Res_x_cell + Res_y_cell), 0.);
    h_sn[k] += DP_total[k] * (1. - melt_mask[k]) * dt_min - S_coeff * dt_min;
  }

#endif
}

//==============================================================================
// For sediment transport
template <class T_type, class U_type>
void computeResidualsTruncated(
    const T_type &u, const T_type &v,
    const unsigned int N_cols, const unsigned int N_rows,
    const unsigned int N, const double c1, const T_type &S_x,
    const T_type &S_y, const double alpha, const double beta,
    const double gamma,
    const U_type &idStaggeredInternalVectHorizontal,
    const U_type &idStaggeredInternalVectVertical,
    const U_type &idStaggeredBoundaryVectWest,
    const U_type &idStaggeredBoundaryVectEast,
    const U_type &idStaggeredBoundaryVectNorth,
    const U_type &idStaggeredBoundaryVectSouth,
#ifdef ENABLE_CUDA
    T_type &Gamma_x_1, T_type &Gamma_x_2,
    T_type &Gamma_y_1, T_type &Gamma_y_2, 
    const T_type& S_x_mod, const double& S_y_mod, cudaStream_t stream = 0
#else
    std::vector<std::array<double, 2>> &Gamma_x,
    std::vector<std::array<double, 2>> &Gamma_y
#endif
 
    ) {

#ifdef ENABLE_CUDA

      computeResidualsTruncatedHorizontal_wrapper(
        idStaggeredInternalVectHorizontal,
	idStaggeredBoundaryVectWest, idStaggeredBoundaryVectEast,
        Gamma_x_1, Gamma_x_2,
	S_x_mod,
        u, c1, stream);
 
      computeResidualsTruncatedVertical_wrapper(
        idStaggeredInternalVectVertical,
	idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectSouth,
        Gamma_y_1, Gamma_y_2,
	S_y_mod,
        v, c1, stream);

#else
  for (unsigned int ii = 0; ii < idStaggeredInternalVectHorizontal.size(); ii++) {
    const auto &Id = idStaggeredInternalVectHorizontal[ii];

    const double coeff_right = c1 * alpha * std::pow(std::abs(S_x[Id]), beta) * // here it is better to compute this vector outside the main loop!
                               u[Id] * (.5 - .5 * signum(u[Id])),

                 coeff_left = c1 * alpha * std::pow(std::abs(S_x[Id]), beta) *
                              u[Id] * (.5 + .5 * signum(u[Id]));

    Gamma_x[Id] = std::array<double, 2>{{coeff_right, coeff_left}};
  }

  for (unsigned int ii = 0; ii < idStaggeredBoundaryVectWest.size(); ii++) {
    const auto &Id = idStaggeredBoundaryVectWest[ii];

    const double coeff_right = c1 * alpha * std::pow(std::abs(S_x[Id]), beta) *
                               u[Id] * (.5 - .5 * signum(u[Id])),

                 coeff_left = c1 * alpha * std::pow(std::abs(S_x[Id]), beta) *
                              u[Id] * (.5 + .5 * signum(u[Id]));

    Gamma_x[Id] = std::array<double, 2>{{coeff_right, coeff_left}};
  }
 
  for (unsigned int ii = 0; ii < idStaggeredBoundaryVectEast.size(); ii++) {
    const auto &Id = idStaggeredBoundaryVectEast[ii];

    const double coeff_right = c1 * alpha * std::pow(std::abs(S_x[Id]), beta) *
                               u[Id] * (.5 - .5 * signum(u[Id])),

                 coeff_left = c1 * alpha * std::pow(std::abs(S_x[Id]), beta) *
                              u[Id] * (.5 + .5 * signum(u[Id]));

    Gamma_x[Id] = std::array<double, 2>{{coeff_right, coeff_left}};
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
#endif
}

//==============================================================================

void make_sparsity_pattern(thrust::host_vector<unsigned int>& idBasinVect,
		std::vector<unsigned int>& basin_mask, 
		thrust::host_vector<unsigned int>& idBasinVectReIndex,
		int *row_offsets, int* columns,
		unsigned int const N_rows, unsigned int const N_cols);

//==============================================================================
#endif
