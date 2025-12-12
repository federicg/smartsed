/*
*******************************************************************************

    Copyright (C) 2019 Politecnico di Milano

    This file is part of SMART-SED.

    SMART-SED is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SMART-SED is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with SMART-SED.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/

/*!
    @file main.cpp
    @brief

    @author      Federico     Gatti        MOX Politecnico di Milano
   <federico.gatti@polimi.it>
    @mantainer

    @supervisors Luca         Bonaventura  MOX Politecnico di Milano
   <luca.bonaventura@polimi.it> Alessandra   Menafoglio   MOX Politecnico di
   Milano               <alessandra.menafoglio@polimi.it> Laura        Longoni
   Applied Geology Politecnico di Milano   <laura.longoni@polimi.it>

    @date 20-06-2020


 */

// i: row    index
// j: column index


#include "utils_H.h"
#include "code_init.h"

//! Parse library
#include "GetPot.hpp"

//! for simple profiling
#include "timing.h"

#include <mpi.h>

#ifndef NO_GPU

// Thrust
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>

// cuSPARSE

#endif

int main(int argc, char **argv) {

  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Reading parameters through GetPot
  GetPot command_line(argc, argv);
  const std::string dataFileName =
      command_line.follow("SMARTSED_input", 2, "-f", "--file");
  GetPot dataFile(dataFileName);

  const std::string file_dir = "../Inputs/";
  std::string temperature_file, ET_model;
  const int totSimNumber = command_line.follow(2, "-sim");
  const int scaling_factor = command_line.follow(2, "-scale");
  const std::string friction_model = dataFile("physics/friction_model", "None");
  const double n_manning = dataFile("physics/n_manning", 0.01);

  const double height_thermometer =
      dataFile("files/meteo_data/height_thermometer", 200.);

  const unsigned int steps_per_hour = dataFile("discretization/steps_per_hour", 10);
  const double max_Days = dataFile("discretization/max_Days", 20.);
  const double starting_day = dataFile("discretization/starting_day", 0);
  const double H_min = dataFile("discretization/H_min", 0.001);
  const double T_thr = dataFile("discretization/T_thr", 0);

  const double time_spacing_temp =
        dataFile("files/meteo_data/time_spacing_temp", 1.);
  const std::string format_temp =
        dataFile("files/meteo_data/format_temp", "arpa");
  const double dt_temp = time_spacing_temp * 3600;

  const bool direct_method = dataFile("linear_solver/direct_method", true);

  const bool is_precipitation =
        dataFile("files/meteo_data/precipitation", true);
  const bool constant_precipitation =
        dataFile("files/meteo_data/constant_precipitation", true);

  const bool save_temporal_sequence =
      dataFile("discretization/save_temporal_sequence", false);
  const bool isNonReflectingBC =
      dataFile("discretization/isNonReflectingBC", false);

  const bool is_sediment_transport =
      dataFile("physics/is_sediment_transport", true);

  const bool spit_out_matrix = dataFile("debug/spit_out_matrix", false);
  const std::string matrix_name = dataFile("debug/matrix_name", "/tmp/matrix_");
  const std::string vector_name = dataFile("debug/vector_name", "/tmp/vector_");
  std::string tmpname = "";

  const bool spit_out_solutions_each_time_step =
      dataFile("debug/spit_out_solutions_each_time_step", false);

  const double frequency_save = dataFile("debug/frequency_save", 24.);

  const std::string precipitation_file =
          dataFile("files/meteo_data/rain_file", " ");

  const double time_spacing_rain =
          dataFile("files/meteo_data/time_spacing_rain", 1.);

  const int number_stations =
          dataFile("files/meteo_data/number_stations", 1);

  const int number_gauges = dataFile("discretization/number_gauges", 1);

  const double nstep = steps_per_hour * max_Days * 24;
  const double t_final = max_Days * 24 * 3600;
  const double dt_DSV_given = t_final / double(nstep);
  double dt_DSV = dt_DSV_given;

  if ((size > totSimNumber && totSimNumber > 0) ||
      (size != 1 && totSimNumber <= 1)) {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  const int chunk_length = totSimNumber / size;
  const int residual =
      std::round((double(totSimNumber) / size - chunk_length) * size);

  std::vector<int> chunk_sim_vec(size);
  chunk_sim_vec.assign(size, chunk_length);

  for (unsigned int i = 0; i < residual; i++) {
    chunk_sim_vec[i] += 1;
  }

  // execute R script
  if (totSimNumber >= 0 && rank == 0) {
    const std::string bashCommand =
        std::string("Rscript "
                    "../Geostatistics/Downscaling_Simulation_SoilGrids/"
                    "Downscaling/DownscalingAitchisonSmartSed_2020.R ") +
        std::to_string(totSimNumber) + " " +
        std::to_string(scaling_factor);
    std::system(bashCommand.c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0)
    std::cout << "mean # of simulations per rank, " << chunk_length
              << std::endl;

  if (rank == 0) {
    std::cout << "------------------------ " << std::endl;
    std::cout << "friction_model         = " << friction_model << std::endl;
    std::cout << "n_manning              = " << n_manning << std::endl;
    std::cout << "steps_per_hour         = " << steps_per_hour << std::endl;
    std::cout << "max_Days               = " << max_Days << std::endl;
    std::cout << "t_final                = " << t_final << " sec."
              << std::endl;
    std::cout << "dt_DSV_given           = " << dt_DSV << " sec."
              << std::endl;
    std::cout << "H_min                  = " << H_min << std::endl;
    std::cout << "------------------------ " << std::endl;
  }


  // +-----------------------------------------------+
  // |                Start MC sim                   |
  // +-----------------------------------------------+


  // Starts the Monte Carlo simulation
  for (int currentSimNumber =
           std::min(totSimNumber, current_start_chunk(rank, chunk_sim_vec) + 1);
       currentSimNumber <=
       std::min(totSimNumber, current_start_chunk(rank + 1, chunk_sim_vec));
       currentSimNumber++) {

    /* Variables living all over the code */

    std::string output_dir, infiltrationModel;

    unsigned int N_rows, N_cols, N, numberOfSteps = 1;

    std::vector<unsigned int> idStaggeredBoundaryVectSouth,
        idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectWest,
        idStaggeredBoundaryVectEast, idStaggeredInternalVectHorizontal,
        idStaggeredInternalVectVertical, idBasinVect, idBasinVectReIndex,

        idStaggeredBoundaryVectSouth_excluded,
        idStaggeredBoundaryVectNorth_excluded,
        idStaggeredBoundaryVectWest_excluded,
        idStaggeredBoundaryVectEast_excluded,
        idStaggeredInternalVectHorizontal_excluded,
        idStaggeredInternalVectVertical_excluded, idBasinVect_excluded,
        idBasinVectReIndex_excluded;

    std::vector<std::array<double, 2>> Gamma_vect_x, Gamma_vect_y;

    std::vector<std::tuple<bool, int>> excluded_ids;
    std::vector<double> additional_source_term;
    std::vector<std::vector<unsigned int>> kk_gauges;

    std::vector<double> basin_mask_Vec, orography, h_G, h_sd, h_sn, S_coeff,
        W_Gav, W_Gav_cum, hydraulic_conductivity, Z_Gav, d_90, Res_x, Res_y, u,
        v, n_x, n_y, u_star, v_star, h_interface_x, h_interface_y, slope_x,
        slope_y, slope_cell, soilMoistureRetention, roughness_vect, eta, H;

    Eigen::VectorXd H_basin, rhs;

    double pixel_size, // meter/pixel
        xllcorner, yllcorner, xllcorner_staggered_u, yllcorner_staggered_u,
        xllcorner_staggered_v, yllcorner_staggered_v, NODATA_value, phi_rad, dt_min, c1_sed, dt_sed, slope_x_max, slope_y_max;

    int iter = 0;
    double c1_DSV_, c2_DSV_, c3_DSV_, minH, maxH;

    /*                   */


    resize_rasters( dataFile, N_rows, N_cols, N,

                    idStaggeredBoundaryVectSouth,
    idStaggeredBoundaryVectNorth, 
    idStaggeredBoundaryVectWest,
    idStaggeredBoundaryVectEast, 
    idStaggeredInternalVectHorizontal,
    idStaggeredInternalVectVertical, 
    idBasinVect, 
    idBasinVectReIndex,

    idStaggeredBoundaryVectSouth_excluded,
    idStaggeredBoundaryVectNorth_excluded,
    idStaggeredBoundaryVectWest_excluded,
    idStaggeredBoundaryVectEast_excluded,
    idStaggeredInternalVectHorizontal_excluded,
    idStaggeredInternalVectVertical_excluded, 
    idBasinVect_excluded,
    idBasinVectReIndex_excluded,

    Gamma_vect_x, 
    Gamma_vect_y,

    excluded_ids,
    additional_source_term,

    basin_mask_Vec, 
    orography, 
    h_G, 
    h_sd, 
    h_sn, 
    S_coeff,
    W_Gav, 
    W_Gav_cum, 
    hydraulic_conductivity, 
    Z_Gav, 
    d_90, 
    Res_x, 
    Res_y, 
    u,
    v, 
    n_x, 
    n_y, 
    u_star, 
    v_star, 
    h_interface_x, 
    h_interface_y, 
    slope_x,
    slope_y, 
    slope_cell, 
    soilMoistureRetention, 
    roughness_vect, 
    eta, 
    H,

    H_basin, rhs,

    pixel_size, // meter/pixel
    xllcorner, 
    yllcorner, 
    xllcorner_staggered_u, 
    yllcorner_staggered_u,
    xllcorner_staggered_v, 
    yllcorner_staggered_v, 
    NODATA_value,
   
    currentSimNumber,
    rank,
    scaling_factor,
    friction_model,
    output_dir, file_dir,
    save_temporal_sequence,
    max_Days,
    temperature_file,
    ET_model,
    infiltrationModel,
    phi_rad,
    dt_min,
    dt_temp, slope_x_max, slope_y_max, number_gauges, kk_gauges );

    // +-----------------------------------------------+
    // |                    Rain                       |
    // +-----------------------------------------------+

    Rain precipitation(infiltrationModel, N,
                       dataFile("files/infiltration/isInitialLoss", false),
                       dataFile("files/infiltration/perc_initialLoss", 0.05),
		       is_precipitation, constant_precipitation,
		       precipitation_file, file_dir, time_spacing_rain,
                       number_stations, max_Days, dataFile,
		       xllcorner, yllcorner, pixel_size,
                       N_rows, N_cols, idBasinVect);

    dt_min = std::min(precipitation.dt_rain, dt_temp);
    for (const auto &kk : hydraulic_conductivity) {
      if (kk * (dt_min / pixel_size) > 1.) {
        if (rank == 0)
          std::cout << "Error! Courant number for gravitational layer is "
                       "greater than 1!"
                    << std::endl;
        exit(-1.);
      }
    }

    // +-----------------------------------------------+
    // |                 Temperature                   |
    // +-----------------------------------------------+

    Temperature temp(file_dir + temperature_file, N, max_Days, T_thr, orography,
                     std::round(max_Days * 24 / time_spacing_temp),
                     steps_per_hour, time_spacing_temp, height_thermometer,
                     format_temp);

    // +-----------------------------------------------+
    // |              Evapotranspiration               |
    // +-----------------------------------------------+

    evapoTranspiration ET(ET_model, N, orography, temp.J, max_Days, phi_rad,
                          height_thermometer);

    // +-----------------------------------------------+
    // |                   Runoff                      |
    // +-----------------------------------------------+
    
    upwind H_interface(H, u, v, idStaggeredInternalVectHorizontal_excluded,
                       idStaggeredBoundaryVectWest_excluded,
                       idStaggeredBoundaryVectEast_excluded,
                       idStaggeredInternalVectVertical_excluded,
                       idStaggeredBoundaryVectNorth_excluded,
                       idStaggeredBoundaryVectSouth_excluded, N_rows, N_cols);

    frictionClass alfa(H_interface.horizontal, H_interface.vertical, u, v,
                       idStaggeredInternalVectHorizontal_excluded,
                       idStaggeredBoundaryVectWest_excluded,
                       idStaggeredBoundaryVectEast_excluded,
                       idStaggeredInternalVectVertical_excluded,
                       idStaggeredBoundaryVectNorth_excluded,
                       idStaggeredBoundaryVectSouth_excluded, friction_model,
                       n_manning, dt_DSV, d_90, roughness_vect, 0., N_rows,
                       N_cols, slope_x, slope_y);

    H_basin.resize(idBasinVect_excluded.size());
    rhs.resize(idBasinVect_excluded.size()); 

    const double c1_min = dt_min / pixel_size;
    static constexpr double gravity = 9.81;

    auto c1_DSV = [](double dt_DSV, double pixel_size) {
      return dt_DSV / pixel_size;
    };

    auto c2_DSV = [](double c1) {
      return gravity * c1;
    };

    auto c3_DSV = [](double c1) {
      return gravity * c1 * c1;
    };

    const double area = std::pow(pixel_size, 2) * 1.e-6; // km^2

    for (unsigned int i = 0; i < H_basin.size(); i++)
      H_basin(i) = 0.;
    for (unsigned int i = 0; i < rhs.size(); i++)
      rhs(i) = 0.;

    Eigen::SparseMatrix<double> A(idBasinVect_excluded.size(), idBasinVect_excluded.size());

    // row, column and value in the Triplet
    std::vector<Eigen::Triplet<double>> coefficients;
    maxH = *std::max_element(H.begin(), H.end());

    dt_DSV = maxdt(u, v, maxH, pixel_size);
    dt_DSV = dt_DSV < dt_DSV_given ? dt_DSV : dt_DSV_given;
    dt_DSV = dt_DSV < t_final ? dt_DSV : t_final;

    c1_DSV_ = c1_DSV(dt_DSV, pixel_size);
    c2_DSV_ = c2_DSV(c1_DSV_);
    c3_DSV_ = c3_DSV(c1_DSV_);

    // +-----------------------------------------------+
    // |             Start the time loop               |
    // +-----------------------------------------------+

    // set initial quantities for the time step adaptation
    auto H_old = H;
    auto H_oldold = H;

    double time = 0., timed = -dt_DSV, timedd = -2. * dt_DSV;
    bool is_last_step = false, check_last = false;

    // +-----------------------------------------------+
    // |   Copy to Device all the needed objects       |
    // +-----------------------------------------------+
    
    // From this part on we just perform all the computations on the device, 
    // the CPU is used only to save the current solution on the disk
    /*
    auto size = 100*sizeof(double); // size in bytes of 100 doubles
    double *v_d;
    cudaMalloc(&v_d, size); // allocate on the device
    double* v_h = (double*)malloc(size);
    cudaMemcpy(v_d, v_h, size, cudaMemcpyHostToDevice);
    cudaMemcpy(v_h, v_d, size, cudaMemcpyDeviceToHost);

    // Host -> Device, if h_std is a std::vector<float>
    thrust::device_vector<float> d_vec(h_std.begin(), h_std.end());
    */

    // Host -> Device, if h_std is a std::vector<float>
    //thrust::device_vector<double> H_device(H.begin(), H.end());


    while (!is_last_step) {

      if (rank == 0) {
        std::cout << "Simulation progress: " << time / t_final * 100 << " %"
                  << " max surface run-off vel. based Courant: "
                  << maxCourant(u, v, c1_DSV_)
                  << " max surface run-off cel. based Courant: "
                  << maxCourant(H, c1_DSV_) << std::endl;
        std::cout << "Current dt, " << dt_DSV << ", given dt, " << dt_DSV_given
                  << std::endl;
      }

      // Compute interface fluxes via upwind method
      H_interface.computeHorizontal();
      H_interface.computeVertical();

      // Compute alfa coefficients
      alfa.f_x();
      alfa.f_y();

      // update only if necessary  --> governed by temperature dynamics, i.e.
      // time_spacing_temp
      if (std::floor(time / dt_temp) > std::floor((time - dt_DSV) / dt_temp)) {
        temp.computeTemperature(std::floor(time / dt_temp), orography,
                                idBasinVect);
      }

      // ET varies daily
      if (std::floor(time / (24. * 3600)) >
          std::floor((time - dt_DSV) / (24. * 3600))) {
        ET.ET(temp.T_dailyMean, temp.T_dailyMin, temp.T_dailyMax,
              std::floor(time / (24 * 3600)), idBasinVect, orography);
      }

      // update only if necessary
      if (std::floor(time / dt_min) > std::floor((time - dt_DSV) / dt_min)) {

        // +-----------------------------------------------+
        // |            Gravitational Layer                |
        // +-----------------------------------------------+

        // h_G
        // vertical and horizontal residuals for Gravitational Layer
        computeResiduals(
            n_x, n_y, N_cols, N_rows, h_G, hydraulic_conductivity,
            idStaggeredInternalVectHorizontal, idStaggeredInternalVectVertical,
            idStaggeredBoundaryVectWest, idStaggeredBoundaryVectEast,
            idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectSouth,
            idBasinVect, h_interface_x, h_interface_y, Res_x, Res_y);
      }

      precipitation.computePrecipitation(time, soilMoistureRetention,
                                         temp.melt_mask, h_G, H, N_rows, N_cols,
                                         idBasinVect);

      // update only if necessary
      if (std::floor(time / dt_min) > std::floor((time - dt_DSV) / dt_min)) {
        for (const unsigned int &k : idBasinVect) {
          S_coeff[k] = 4.62e-10 * h_sn[k] * (temp.T_raster[k] - T_thr) *
                       temp.melt_mask[k];
        }

        for (const auto &k : idBasinVect) {
          h_G[k] += (S_coeff[k] - ET.ET_vec[k]) * dt_min +
                    precipitation.DP_infiltrated[k] * dt_min -
                    c1_min * (Res_x[k] + Res_y[k]);
          h_G[k] *= (h_G[k] >= 0); // to account for evapotranspiration
        }

        // +-----------------------------------------------+
        // |              Snow Accumulation                |
        // +-----------------------------------------------+

        for (const auto &k : idBasinVect) {
          const auto snow_acc =
              precipitation.DP_total[k] * (1. - temp.melt_mask[k]) * dt_min -
              S_coeff[k] * dt_min;
          h_sn[k] += snow_acc;
        }
      }

      // +-----------------------------------------------+
      // |               De Saint Venant                 |
      // +-----------------------------------------------+

      // tic();
      //  fill u_star and v_star with a Bilinear Interpolation
      bilinearInterpolation(u, v, u_star, v_star, N_rows, N_cols, dt_DSV,
                            pixel_size,
                            idStaggeredInternalVectHorizontal_excluded,
                            idStaggeredInternalVectVertical_excluded,
                            idStaggeredBoundaryVectWest_excluded,
                            idStaggeredBoundaryVectEast_excluded,
                            idStaggeredBoundaryVectNorth_excluded,
                            idStaggeredBoundaryVectSouth_excluded);

      // toc ("bilinearInterpolation");

      // tic();

      coefficients.reserve(
          idBasinVect_excluded.size() +
          4 * idStaggeredInternalVectHorizontal_excluded.size() +
          4 * idStaggeredInternalVectVertical_excluded.size());

      additional_source_term.assign(N, 0.);

      buildMatrix(H_interface.horizontal, H_interface.vertical, orography,
                  u_star, v_star, u, v, H, N_cols, N_rows, N, c1_DSV_, c3_DSV_,
                  0, 
                  precipitation.DP_cumulative, dt_DSV, alfa.alfa_x, alfa.alfa_y,
                  idStaggeredInternalVectHorizontal_excluded,
                  idStaggeredInternalVectVertical_excluded,
                  idStaggeredBoundaryVectWest_excluded,
                  idStaggeredBoundaryVectEast_excluded,
                  idStaggeredBoundaryVectNorth_excluded,
                  idStaggeredBoundaryVectSouth_excluded, idBasinVect_excluded,
                  idBasinVect, idStaggeredInternalVectHorizontal,
                  idStaggeredInternalVectVertical, idBasinVectReIndex_excluded,
                  isNonReflectingBC, true, excluded_ids, additional_source_term,
                  coefficients, rhs);

      A.setFromTriplets(coefficients.begin(), coefficients.end());
      A.makeCompressed();
      coefficients.clear();

      // toc ("assemble matrix");

      if (spit_out_matrix) {
        tmpname = matrix_name + std::to_string(iter);
        std::cout << "saving " << tmpname << std::endl;
        Eigen::saveMarket(A, tmpname, false);
        tmpname = vector_name + std::to_string(iter);
        std::cout << "saving " << tmpname << std::endl;
        Eigen::saveMarketVector_lis(rhs, tmpname);
      }

      if (direct_method) // Direct Sparse method: Cholesky being A spd
      {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
          std::cout << "Decomposition Failed" << std::endl;
          exit(-1.);
        }
        H_basin = solver.solve(rhs);
      } else {

        double tol = 1.e-6;
        int result, maxit = 15000;

        // IML++
        Eigen::IncompleteCholesky<double> IC(
            A); // Create I cholesky preconditioner
        result =
            LinearAlgebra::CG(A, H_basin, rhs, IC, maxit, tol); // Solve system

        if (rank == 0) {
          std::cout << " IML++ CG " << std::endl;
          std::cout << "CG flag = " << result << std::endl;
          std::cout << "iterations performed: " << maxit << std::endl;
          std::cout << "tolerance achieved  : " << tol << std::endl;
        }
      }

      minH = H_basin.minCoeff();
      maxH = H_basin.maxCoeff();
      if (rank == 0)
        std::cout << "min H: " << minH << " max H: " << maxH << std::endl;

      if (minH < -H_min) {
        dt_DSV = dt_DSV / 10.;
        if (dt_DSV == 0) {
          if (rank == 0)
            std::cout << "dt has gone to zero, sorry, STOP!" << std::endl;
          exit(-1.);
        }
        c1_DSV_ = c1_DSV(dt_DSV, pixel_size);
        c2_DSV_ = c2_DSV(c1_DSV_);
        c3_DSV_ = c3_DSV(c1_DSV_);
        continue;
      }

      for (const unsigned int &Id : idBasinVect_excluded) {
        const unsigned int IDreIndex = idBasinVectReIndex_excluded[Id];
        H_oldold[Id] = H_old[Id];
        H_old[Id] = H[Id];
        H[Id] = std::abs(H_basin(IDreIndex));
        eta[Id] = H[Id] + orography[Id];
      }

      // +-----------------------------------------------+
      // |                  Update u, v                  |
      // +-----------------------------------------------+

      updateVel(u, v, u_star, v_star, alfa.alfa_x, alfa.alfa_y, N_rows, N_cols,
                c2_DSV_, 0, eta, H, orography,
                idStaggeredInternalVectHorizontal_excluded,
                idStaggeredInternalVectVertical_excluded,
                idStaggeredBoundaryVectWest_excluded,
                idStaggeredBoundaryVectEast_excluded,
                idStaggeredBoundaryVectNorth_excluded,
                idStaggeredBoundaryVectSouth_excluded, isNonReflectingBC);

      // +-----------------------------------------------+
      // |             Sediment Transport                |
      // +-----------------------------------------------+

      if (is_sediment_transport) {
        // tic();
        constexpr double alfa_coeff = 2.5, beta_coeff = 1.6, gamma_coeff = 1.;

        dt_sed = compute_dt_sediment(alfa_coeff, beta_coeff, slope_x_max,
                                     slope_y_max, u, v, pixel_size, dt_DSV,
                                     &numberOfSteps);
        c1_sed = dt_sed / pixel_size;

        // vertical and horizontal residuals truncated for Sediment Transport
        computeResidualsTruncated(
            u, v, N_cols, N_rows, N, c1_sed, slope_x, slope_y,
            alfa_coeff,  // alfa
            beta_coeff,  // beta
            gamma_coeff, // gamma
            idStaggeredInternalVectHorizontal_excluded,
            idStaggeredInternalVectVertical_excluded,
            idStaggeredBoundaryVectWest_excluded,
            idStaggeredBoundaryVectEast_excluded,
            idStaggeredBoundaryVectNorth_excluded,
            idStaggeredBoundaryVectSouth_excluded, Gamma_vect_x, Gamma_vect_y);

        additional_source_term.assign(N, 0.);

        for (const auto &k : idBasinVect) {
          const auto &current_tuple = excluded_ids[k];
          if (std::get<0>(current_tuple)) {
            const auto &k_pour = std::get<1>(current_tuple);
            if (k_pour >= 0) {
              additional_source_term[k_pour] +=
                  1.e-3 * M_PI * Z_Gav[k] *
                  std::sqrt(std::abs((.1 + .1 * temp.T_raster[k]) *
                                     temp.melt_mask[k])) *
                  precipitation.DP_total[k] * dt_sed;
            }
          }
        }

        for (const unsigned int &k : idBasinVect_excluded) {
          // 1.e-3 is the conversion factor, look at EPM theory
          W_Gav[k] = 1.e-3 * M_PI * Z_Gav[k] *
                         std::sqrt(std::abs((.1 + .1 * temp.T_raster[k]) *
                                            temp.melt_mask[k])) *
                         precipitation.DP_total[k] * dt_sed +
                     additional_source_term[k];
        }

        if (rank == 0)
          std::cout << "# steps for solid transport, " << numberOfSteps
                    << std::endl;

        for (unsigned int kk = 0; kk < numberOfSteps; kk++) {

          additional_source_term.assign(N, 0.);

          // assemble horizontal fluxes
          for (const unsigned int &Id : idStaggeredInternalVectHorizontal_excluded) {
            const unsigned int i = Id / (N_cols + 1), IDeast = Id - i,
                       IDwest = Id - i - 1;

            const double &h_left = h_sd[IDwest], &h_right = h_sd[IDeast];

            h_interface_x[Id] =
                Gamma_vect_x[Id][0] * h_right + Gamma_vect_x[Id][1] * h_left;
          }

          for (const unsigned int &Id : idStaggeredBoundaryVectWest_excluded) {
            const unsigned int i = Id / (N_cols + 1);

            const double h_left = 0, // 0,
                h_right = h_sd[Id - i];

            h_interface_x[Id] =
                Gamma_vect_x[Id][0] * h_right + Gamma_vect_x[Id][1] * h_left;
          }

          for (const unsigned int &Id : idStaggeredBoundaryVectEast_excluded) {

            const unsigned int i = Id / (N_cols + 1);

            const double h_left = h_sd[Id - i - 1], h_right = 0;

            h_interface_x[Id] =
                Gamma_vect_x[Id][0] * h_right + Gamma_vect_x[Id][1] * h_left;
          }

          for (const unsigned int &Id : idBasinVect_excluded) {
            const unsigned int i = Id / N_cols;

            Res_x[Id] = h_interface_x[Id + 1 + i] - h_interface_x[Id + i];
          }

          // assemble vertical fluxes
          for (const unsigned int &Id : idStaggeredInternalVectVertical_excluded) {
            const unsigned int IDsouth = Id, IDnorth = Id - N_cols;

            const double &h_left = h_sd[IDnorth], &h_right = h_sd[IDsouth];

            h_interface_y[Id] =
                Gamma_vect_y[Id][0] * h_right + Gamma_vect_y[Id][1] * h_left;
          }

          for (const unsigned int &Id : idStaggeredBoundaryVectNorth_excluded) {

            const double h_left = 0, // 0
                h_right = h_sd[Id];

            h_interface_y[Id] =
                Gamma_vect_y[Id][0] * h_right + Gamma_vect_y[Id][1] * h_left;
          }

          for (const unsigned int &Id : idStaggeredBoundaryVectSouth_excluded) {

            const double h_left = h_sd[Id - N_cols], h_right = 0;

            h_interface_y[Id] =
                Gamma_vect_y[Id][0] * h_right + Gamma_vect_y[Id][1] * h_left;
          }

          for (const unsigned int &Id : idBasinVect_excluded) {
            Res_y[Id] = h_interface_y[Id + N_cols] - h_interface_y[Id];
          }

          // assemble the source term
          for (const auto &Id : idStaggeredInternalVectHorizontal) {

            const unsigned int i = Id / (N_cols + 1), IDleft = Id - i - 1,
                       IDright = Id - i;

            if (std::get<0>(excluded_ids[IDleft])) {
              const auto &k_pour = std::get<1>(excluded_ids[IDleft]);
              if (k_pour >= 0) {
                additional_source_term[k_pour] +=
                    h_interface_x[Id] * std::abs(u[Id]) * c1_sed;
              }
            }

            if (std::get<0>(excluded_ids[IDright])) {
              const auto &k_pour = std::get<1>(excluded_ids[IDright]);
              if (k_pour >= 0) {
                additional_source_term[k_pour] +=
                    h_interface_x[Id] * std::abs(u[Id]) * c1_sed;
              }
            }
          }

          for (const auto &Id : idStaggeredInternalVectVertical) {

            const unsigned int IDleft = Id - N_cols, IDright = Id;

            if (std::get<0>(excluded_ids[IDleft])) {
              const auto &k_pour = std::get<1>(excluded_ids[IDleft]);
              if (k_pour >= 0) {
                additional_source_term[k_pour] +=
                    h_interface_y[Id] * std::abs(v[Id]) * c1_sed;
              }
            }

            if (std::get<0>(excluded_ids[IDright])) {
              const auto &k_pour = std::get<1>(excluded_ids[IDright]);
              if (k_pour >= 0) {
                additional_source_term[k_pour] +=
                    h_interface_y[Id] * std::abs(v[Id]) * c1_sed;
              }
            }
          }

          // Update the solution
          for (const unsigned int &Id : idBasinVect_excluded) {
            h_sd[Id] += -(Res_x[Id] + Res_y[Id]) + W_Gav[Id] +
                        additional_source_term[Id];
          }
        }
        // toc ("advance");

      } // End of if(sediment_transport)

      // +-----------------------------------------------+
      // |               Update time                     |
      // +-----------------------------------------------+

      timedd = timed;
      timed = time;
      time += dt_DSV;

      if (check_last) {
        time = t_final;
        is_last_step = true;
      }

      // +-----------------------------------------------+
      // |             Save The Raster Solution          |
      // +-----------------------------------------------+

      if (save_temporal_sequence) {

        for (int number = 1; number <= number_gauges; number++) {
          std::string filename_x = "discretization/X_gauges_",
                      filename_y = "discretization/Y_gauges_";

          filename_x += std::to_string(number);
          filename_y += std::to_string(number);

          const double X_gauges = dataFile(filename_x.c_str(), 0.);
          const double Y_gauges = dataFile(filename_y.c_str(), 0.);

          const Vector2D XX_gauges(std::array<double, 2>{{X_gauges, Y_gauges}});

          double H_current = 0., H_candidate = 0., mass_flux_candidate = 0.,
                 solid_flux_candidate = 0.;
          unsigned int kk_gauges_max = 0;
          for (const auto & candidate : kk_gauges[number - 1]) {
            const unsigned int i = candidate / N_cols;
            const auto & cc = H[candidate];

            const auto velo = std::sqrt(
                std::pow(((v[candidate] + v[candidate + N_cols]) * .5), 2.) +
                std::pow(((u[candidate + i] + u[candidate + i + 1]) * .5), 2.));

            H_candidate += cc;
            mass_flux_candidate += cc * velo;
            solid_flux_candidate += h_sd[candidate] * velo;

            if (cc > H_current) {
              kk_gauges_max = candidate;
              H_current = cc;
            }
          }
          const unsigned int i = kk_gauges_max / N_cols;

          H_candidate /= kk_gauges[number - 1].size();
          mass_flux_candidate /= kk_gauges[number - 1].size();
          solid_flux_candidate /= kk_gauges[number - 1].size();

          saveTemporalSequence(
              XX_gauges, time,
              output_dir + "timesteps_" + std::to_string(number), dt_DSV);

          saveTemporalSequence(XX_gauges, time,
                               output_dir + "waterSurfaceHeight_" +
                                   std::to_string(number),
                               H_candidate);
          saveTemporalSequence(XX_gauges, time,
                               output_dir + "waterSurfaceMassFlux_" +
                                   std::to_string(number),
                               mass_flux_candidate);
          saveTemporalSequence(XX_gauges, time,
                               output_dir + "SolidFlux_" +
                                   std::to_string(number),
                               solid_flux_candidate);

          saveTemporalSequence(XX_gauges, time,
                               output_dir + "waterSurfaceHeightmax_" +
                                   std::to_string(number),
                               H[kk_gauges_max]);

          saveTemporalSequence(
              XX_gauges, time,
              output_dir + "waterSurfaceMassFluxmax_" + std::to_string(number),
              H[kk_gauges_max] *
                  std::sqrt(
                      std::pow(
                          ((v[kk_gauges_max] + v[kk_gauges_max + N_cols]) / 2.),
                          2.) +
                      std::pow(
                          ((u[kk_gauges_max + i] + u[kk_gauges_max + i + 1]) /
                           2.),
                          2.)));

          saveTemporalSequence(
              XX_gauges, time,
              output_dir + "SolidFluxmax_" + std::to_string(number),
              h_sd[kk_gauges_max] *
                  std::sqrt(
                      std::pow(
                          ((v[kk_gauges_max] + v[kk_gauges_max + N_cols]) / 2.),
                          2.) +
                      std::pow(
                          ((u[kk_gauges_max + i] + u[kk_gauges_max + i + 1]) /
                           2.),
                          2.)));
        }
      }

      // sediment production zones,
      for (const unsigned int &k : idBasinVect) {
        W_Gav_cum[k] += W_Gav[k] * (dt_DSV / dt_sed);
      }

      if (spit_out_solutions_each_time_step) {
        iter++;

        saveSolution(output_dir + "u_", "u", N_rows, N_cols,
                     xllcorner_staggered_u, yllcorner_staggered_u, pixel_size,
                     NODATA_value, iter, u, v, H);
        saveSolution(output_dir + "v_", "v", N_rows, N_cols,
                     xllcorner_staggered_v, yllcorner_staggered_v, pixel_size,
                     NODATA_value, iter, u, v, H);
        saveSolution(output_dir + "H_", " ", N_rows, N_cols, xllcorner,
                     yllcorner, pixel_size, NODATA_value, iter, u, v, H);
        saveSolution(output_dir + "hsd_", " ", N_rows, N_cols, xllcorner,
                     yllcorner, pixel_size, NODATA_value, iter, u, v, h_sd);
        saveSolution(output_dir + "w_cum_", " ", N_rows, N_cols, xllcorner,
                     yllcorner, pixel_size, NODATA_value, iter, u, v,
                     W_Gav_cum);
        saveSolution(output_dir + "hG_", " ", N_rows, N_cols, xllcorner,
                     yllcorner, pixel_size, NODATA_value, iter, u, v, h_G);
        saveSolution(output_dir + "hsn_", " ", N_rows, N_cols, xllcorner,
                     yllcorner, pixel_size, NODATA_value, iter, u, v, h_sn);

      } else {

        if (std::floor(time / (frequency_save * 3600)) >
            std::floor((time - dt_DSV) / (frequency_save * 3600))) {
          iter++;

          const auto &currentDay =
              iter; // std::floor ( time / (frequency_save * 3600) )

          if (rank == 0)
            std::cout << "Saving solution ..., current saving " << iter
                      << std::endl;

          saveSolution(output_dir + "u_", "u", N_rows, N_cols,
                       xllcorner_staggered_u, yllcorner_staggered_u, pixel_size,
                       NODATA_value, currentDay, u, v, H);
          saveSolution(output_dir + "v_", "v", N_rows, N_cols,
                       xllcorner_staggered_v, yllcorner_staggered_v, pixel_size,
                       NODATA_value, currentDay, u, v, H);
          saveSolution(output_dir + "H_", " ", N_rows, N_cols, xllcorner,
                       yllcorner, pixel_size, NODATA_value, currentDay, u, v,
                       H);
          saveSolution(output_dir + "hsd_", " ", N_rows, N_cols, xllcorner,
                       yllcorner, pixel_size, NODATA_value, currentDay, u, v,
                       h_sd);
          saveSolution(output_dir + "w_cum_", " ", N_rows, N_cols, xllcorner,
                       yllcorner, pixel_size, NODATA_value, currentDay, u, v,
                       W_Gav_cum);
          saveSolution(output_dir + "hG_", " ", N_rows, N_cols, xllcorner,
                       yllcorner, pixel_size, NODATA_value, currentDay, u, v,
                       h_G);
          saveSolution(output_dir + "hsn_", " ", N_rows, N_cols, xllcorner,
                       yllcorner, pixel_size, NODATA_value, currentDay, u, v,
                       h_sn);
        }
      }

      // +-----------------------------------------------+
      // |               Update time step                |
      // +-----------------------------------------------+

      dt_DSV = maxdt(u, v, maxH, pixel_size);
      double dt_DSV_min = dt_DSV * .5;

      compute_dt_adaptive(H, H_old, H_oldold, idBasinVect_excluded, dt_DSV,
                          1.e-5, time, timed, timedd);

      dt_DSV = std::min(std::max(dt_DSV, dt_DSV_min), dt_DSV_given);

      c1_DSV_ = c1_DSV(dt_DSV, pixel_size);
      c2_DSV_ = c2_DSV(c1_DSV_);
      c3_DSV_ = c3_DSV(c1_DSV_);

      if ((time + dt_DSV) > t_final) {
        dt_DSV = t_final - time;
        check_last = true;
      }

    } // End Time Loop

  } // End Monte Carlo loop

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    print_timing_report();
  }
  MPI_Finalize();

  return (EXIT_SUCCESS);
}
