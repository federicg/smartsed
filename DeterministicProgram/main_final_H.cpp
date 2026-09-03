/*
    @author Federico Gatti   <federico.gatti@math.ethz.ch>
    The GPU implementation works just for CUDA, so NVIDIA GPUs,
    in the future there could be an implemnetation to work with HIP for AMD GPUs
*/

//! Utility CPU header
#include "utils_H.h"

//! for simple profiling
#include "timing.h"

//! MPI
#include <mpi.h>

//! wall-clock benchmarking of the time loop
#include <chrono>

//! std::sort (CSR row sorting for the red-black reindex)
#include <algorithm>

#ifdef ENABLE_CUDA
//! NetCDF
#include <netcdf.h>
#endif

//
#ifdef ENABLE_CUDA
#define CUDA_STREAM(s)        , s   // for functions that have other args before stream
#define CUDA_STREAM_ONLY(s)   s     // for functions where stream is the only arg
#else
#define CUDA_STREAM(s)
#define CUDA_STREAM_ONLY(s)
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

  const unsigned int steps_per_hour =
      dataFile("discretization/steps_per_hour", 10);
  const unsigned int max_Days = dataFile("discretization/max_Days", 20);
  const double starting_day = dataFile("discretization/starting_day", 0);
  const double H_min = dataFile("discretization/H_min", 0.001);
  const double T_thr = dataFile("discretization/T_thr", 0);

  const double time_spacing_temp =
      dataFile("files/meteo_data/time_spacing_temp", 1.);
  const std::string format_temp =
      dataFile("files/meteo_data/format_temp", "arpa");
  const double dt_temp = time_spacing_temp * 3600;

  const bool direct_method = dataFile("linear_solver/direct_method", true);

  // GPU-only: whether gpu_CG applies the IC(0) preconditioner (cusparseDcsric02
  // + the two SpSV triangular solves) or runs plain CG. The CPU reference
  // (LinearAlgebra::CG, cg.hpp) always uses the analogous Eigen
  // IncompleteCholesky preconditioner, so this should normally stay true --
  // it exists mainly to reproduce/compare against an unpreconditioned run.
  const bool use_preconditioner =
      dataFile("linear_solver/use_preconditioner", true);

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

  const bool restart_soilMoisture =
      dataFile("files/initial_conditions/restart_soilMoisture", false);

  const bool spit_out_matrix = dataFile("debug/spit_out_matrix", false);
  const std::string matrix_name = dataFile("debug/matrix_name", "/tmp/matrix_");
  const std::string vector_name = dataFile("debug/vector_name", "/tmp/vector_");
  std::string tmpname = "";

  const bool spit_out_solutions_each_time_step =
      dataFile("debug/spit_out_solutions_each_time_step", false);

  // DEBUG: if true, exit right after completing the first time step
  const bool stop_after_first_step =
      dataFile("debug/stop_after_first_step", false);

  // BENCH: if > 0, stop the time loop after this many accepted steps and
  //        print a wall-clock / field-fingerprint report (0 = run normally).
  //        Command-line "-steps N" overrides the datafile value.
  const unsigned int max_bench_steps = command_line.follow(
      static_cast<int>(dataFile("debug/max_steps", 0)), "-steps");

  const double frequency_save = dataFile("debug/frequency_save", 24.);

  const std::string precipitation_file =
      dataFile("files/meteo_data/rain_file", " ");

  const double time_spacing_rain =
      dataFile("files/meteo_data/time_spacing_rain", 1.);

  const int number_stations = dataFile("files/meteo_data/number_stations", 1);

  const int number_gauges = dataFile("discretization/number_gauges", 1);

  // BENCH counters (see max_bench_steps / stop_after_first_step)
  unsigned int bench_completed_steps = 0;
  unsigned int bench_loop_iters = 0;

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
  if (totSimNumber >= 0 && rank == 0 && !restart_soilMoisture) {
    const std::string bashCommand =
        std::string("Rscript "
                    "../Geostatistics/"
                    "DownscalingAitchisonSmartSed_2020.R ") +
        std::to_string(totSimNumber) + " " + std::to_string(scaling_factor);
    std::system(bashCommand.c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0)
    std::cout << "mean # of simulations per MPI rank, " << chunk_length
              << std::endl;

  if (rank == 0) {
    std::cout << "------------------------ " << std::endl;
    std::cout << "friction_model         = " << friction_model << std::endl;
    std::cout << "n_manning              = " << n_manning << std::endl;
    std::cout << "steps_per_hour         = " << steps_per_hour << std::endl;
    std::cout << "max_Days               = " << max_Days << std::endl;
    std::cout << "t_final                = " << t_final << " sec." << std::endl;
    std::cout << "dt_DSV_given           = " << dt_DSV << " sec." << std::endl;
    std::cout << "H_min                  = " << H_min << std::endl;
    std::cout << "------------------------ " << std::endl;
  }

  // +-----------------------------------------------+
  // |                Start MC sim                   |
  // +-----------------------------------------------+

#ifdef ENABLE_CUDA
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    static constexpr int N_STREAMS = 6;  // one per independent kernel launch
    const int max_streams = prop.concurrentKernels ? N_STREAMS : 1;  // fallback to 1 if no concurrency

    std::vector<cudaStream_t> streams(max_streams);
    for (auto &s : streams)
      cudaStreamCreate(&s);

    // use by index
    auto S = [&](int i) { return streams[i % max_streams]; };
#endif


  // Starts the Monte Carlo simulation
  for (int currentSimNumber =
           std::min(totSimNumber, current_start_chunk(rank, chunk_sim_vec) + 1);
       currentSimNumber <=
       std::min(totSimNumber, current_start_chunk(rank + 1, chunk_sim_vec));
       currentSimNumber++) {

    /* Variables living all over the code */

    std::string output_dir, infiltrationModel;

    unsigned int N_rows, N_cols, N, numberOfSteps = 1; 
    size_t iter = 0;

    std::vector<unsigned int> basin_mask_Vec;

    std::vector<std::tuple<bool, int>> excluded_ids;
    std::vector<std::vector<unsigned int>> kk_gauges;

#ifdef ENABLE_CUDA
/* define here the linear system stuffs
 * use the cuSPARSE lib
*/
    thrust::host_vector<double> Gamma_vect_x_1_pot, Gamma_vect_x_2_pot, Gamma_vect_y_1_pot, Gamma_vect_y_2_pot;

    thrust::host_vector<double> orography_pot, h_G_pot, h_sd_pot, h_sn_pot,
        W_Gav_pot, W_Gav_cum_pot, hydraulic_conductivity_pot, Z_Gav_pot, d_90, u_pot,
        v_pot, n_x_pot, n_y_pot, u_star_pot, v_star_pot, h_interface_x_pot, h_interface_y_pot, slope_x_pot,
        slope_y_pot, slope_x_mod_pot, slope_y_mod_pot, soilMoistureRetention_pot, roughness_vect, eta_pot, 
	H_pot, additional_source_term_pot;

    thrust::host_vector<unsigned int> idStaggeredBoundaryVectSouth_pot,
        idStaggeredBoundaryVectNorth_pot, idStaggeredBoundaryVectWest_pot,
        idStaggeredBoundaryVectEast_pot, idStaggeredInternalVectHorizontal_pot,
        idStaggeredInternalVectVertical_pot, idBasinVect_pot, 

        idStaggeredBoundaryVectSouth_excluded_pot,
        idStaggeredBoundaryVectNorth_excluded_pot,
        idStaggeredBoundaryVectWest_excluded_pot,
        idStaggeredBoundaryVectEast_excluded_pot,
        idStaggeredInternalVectHorizontal_excluded_pot,
        idStaggeredInternalVectVertical_excluded_pot, idBasinVect_excluded_pot,
        idBasinVectReIndex_excluded_pot;

#else
    std::vector<std::array<double, 2>> Gamma_vect_x, Gamma_vect_y;

    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd H_basin, rhs;
    std::vector<Eigen::Triplet<double>> coefficients;

    std::vector<double> orography_pot, h_G_pot, h_sd_pot, h_sn_pot,
        W_Gav_pot, W_Gav_cum_pot, hydraulic_conductivity_pot, Z_Gav_pot, d_90, u_pot,
        v_pot, n_x_pot, n_y_pot, u_star_pot, v_star_pot, h_interface_x_pot, h_interface_y_pot, slope_x_pot,
        slope_y_pot, soilMoistureRetention_pot, roughness_vect, eta_pot, 
	H_pot, additional_source_term_pot;

    std::vector<unsigned int> idStaggeredBoundaryVectSouth_pot,
        idStaggeredBoundaryVectNorth_pot, idStaggeredBoundaryVectWest_pot,
        idStaggeredBoundaryVectEast_pot, idStaggeredInternalVectHorizontal_pot,
        idStaggeredInternalVectVertical_pot, idBasinVect_pot,

        idStaggeredBoundaryVectSouth_excluded_pot,
        idStaggeredBoundaryVectNorth_excluded_pot,
        idStaggeredBoundaryVectWest_excluded_pot,
        idStaggeredBoundaryVectEast_excluded_pot,
        idStaggeredInternalVectHorizontal_excluded_pot,
        idStaggeredInternalVectVertical_excluded_pot, idBasinVect_excluded_pot,
        idBasinVectReIndex_excluded_pot;

#endif

    double pixel_size, // meter/pixel
        xllcorner, yllcorner, xllcorner_staggered_u, yllcorner_staggered_u,
        xllcorner_staggered_v, yllcorner_staggered_v, NODATA_value, phi_rad,
        dt_min, c1_sed, dt_sed, slope_x_max, slope_y_max,
        c1_DSV_, c2_DSV_, c3_DSV_, minH, maxH, time, timed, timedd, area, c1_min;

    bool is_last_step = false, check_last = false;

    constexpr double alfa_coeff = 2.5, beta_coeff = 1.6, gamma_coeff = 1., gravity = 9.81;
    auto c1_DSV = [](double dt_DSV, double pixel_size) { return dt_DSV / pixel_size; };
    auto c2_DSV = [](double c1) { return gravity * c1; };
    auto c3_DSV = [](double c1) { return gravity * c1 * c1; };

/*a naive implementation of transposition
   for (int i=0; i<N_rows; i++) {
	   for (int j=0; j<N_cols; j++) {
		   H_transpose[j*N_rows + i] = H[i*N_cols + j]; // H is in row-major
	   }
   }

*/

    /* End of (Variables living all over the code) */

    resize_rasters(
        dataFile, N_rows, N_cols, N,

        idStaggeredBoundaryVectSouth_pot, idStaggeredBoundaryVectNorth_pot,
        idStaggeredBoundaryVectWest_pot, idStaggeredBoundaryVectEast_pot,
        idStaggeredInternalVectHorizontal_pot, idStaggeredInternalVectVertical_pot,
        idBasinVect_pot, 

        idStaggeredBoundaryVectSouth_excluded_pot,
        idStaggeredBoundaryVectNorth_excluded_pot,
        idStaggeredBoundaryVectWest_excluded_pot,
        idStaggeredBoundaryVectEast_excluded_pot,
        idStaggeredInternalVectHorizontal_excluded_pot,
        idStaggeredInternalVectVertical_excluded_pot, idBasinVect_excluded_pot,
        idBasinVectReIndex_excluded_pot,

#ifdef ENABLE_CUDA
	Gamma_vect_x_1_pot, Gamma_vect_x_2_pot, Gamma_vect_y_1_pot, Gamma_vect_y_2_pot,
	slope_x_mod_pot, slope_y_mod_pot, alfa_coeff, beta_coeff,
#else	
        Gamma_vect_x, Gamma_vect_y,
#endif

        excluded_ids, additional_source_term_pot,

        basin_mask_Vec, orography_pot, h_G_pot, h_sd_pot, h_sn_pot, W_Gav_pot, W_Gav_cum_pot,
        hydraulic_conductivity_pot, Z_Gav_pot, d_90, u_pot, v_pot, n_x_pot, n_y_pot,
        u_star_pot, v_star_pot, h_interface_x_pot, h_interface_y_pot, slope_x_pot, slope_y_pot,
        soilMoistureRetention_pot, roughness_vect, eta_pot, H_pot,

        pixel_size, // meter/pixel
        xllcorner, yllcorner, xllcorner_staggered_u, yllcorner_staggered_u,
        xllcorner_staggered_v, yllcorner_staggered_v, NODATA_value,

        currentSimNumber, rank, scaling_factor, friction_model, output_dir,
        file_dir, save_temporal_sequence, max_Days, temperature_file, ET_model,
        infiltrationModel, phi_rad, dt_min, dt_temp, slope_x_max, slope_y_max,
        number_gauges, kk_gauges);

    // +-----------------------------------------------+
    // |   Copy to Device all the needed objects       |
    // +-----------------------------------------------+

    // From this part on we just perform all the computations on the device,
    // the CPU is used only to save the current solution on the disk
#ifdef ENABLE_CUDA

    // Device copy of dynamic variables
    thrust::device_vector<double> H = H_pot;
    thrust::device_vector<double> u = u_pot;
    thrust::device_vector<double> v = v_pot;
    thrust::device_vector<double> h_G = h_G_pot;
    thrust::device_vector<double> h_sd = h_sd_pot;
    thrust::device_vector<double> h_sn = h_sn_pot;
    thrust::device_vector<double> W_Gav = W_Gav_pot;
    thrust::device_vector<double> u_star = u_star_pot;
    thrust::device_vector<double> v_star = v_star_pot;
    thrust::device_vector<double> h_interface_x = h_interface_x_pot;
    thrust::device_vector<double> h_interface_y = h_interface_y_pot;
    thrust::device_vector<double> eta = eta_pot;

    // Move to const double vectors
    const thrust::device_vector<double> orography = orography_pot;
    const thrust::device_vector<double> Z_Gav = Z_Gav_pot; 
    const thrust::device_vector<double> n_x = n_x_pot; 
    const thrust::device_vector<double> n_y = n_y_pot;
    const thrust::device_vector<double> slope_x = slope_x_pot;
    const thrust::device_vector<double> slope_y = slope_y_pot;  
    const thrust::device_vector<double> hydraulic_conductivity = hydraulic_conductivity_pot;
    const thrust::device_vector<double> soilMoistureRetention = soilMoistureRetention_pot;

    // Move to const ids
    const thrust::device_vector<unsigned int> idStaggeredInternalVectHorizontal_excluded = idStaggeredInternalVectHorizontal_excluded_pot;
    const thrust::device_vector<unsigned int> idStaggeredBoundaryVectWest_excluded = idStaggeredBoundaryVectWest_excluded_pot;
    const thrust::device_vector<unsigned int> idStaggeredBoundaryVectEast_excluded = idStaggeredBoundaryVectEast_excluded_pot;
    const thrust::device_vector<unsigned int> idStaggeredInternalVectVertical_excluded = idStaggeredInternalVectVertical_excluded_pot;
    const thrust::device_vector<unsigned int> idStaggeredBoundaryVectNorth_excluded = idStaggeredBoundaryVectNorth_excluded_pot;
    const thrust::device_vector<unsigned int> idStaggeredBoundaryVectSouth_excluded = idStaggeredBoundaryVectSouth_excluded_pot;
    const thrust::device_vector<unsigned int> idBasinVect_excluded = idBasinVect_excluded_pot;
    const thrust::device_vector<unsigned int> idBasinVectReIndex_excluded = idBasinVectReIndex_excluded_pot;

    // Concatenated id list + a per-entry membership tag (0=internal,
    // 1=west/north, 2=east/south) for each staggered direction, built once
    // here from the already-classified _excluded lists. The DSV time-loop
    // kernels that used to take three separate id lists (upwind, gravitational
    // residuals, sediment residuals, updateVel, buildMatrix boundary) instead
    // take idAll+tag and run a single tag-dispatched launch. Membership is
    // pure concatenation -- nothing is reclassified -- which matters because
    // idStaggeredBoundaryVectNorth in particular mixes two index conventions
    // (interior mask boundary vs. literal top grid edge) that are not safe to
    // re-infer per face from basin_mask.
    thrust::device_vector<unsigned int> idHorizontalAll_excluded, horizontalTag_excluded,
                                        idVerticalAll_excluded, verticalTag_excluded;
    build_merged_ids_and_tags(idStaggeredInternalVectHorizontal_excluded,
                              idStaggeredBoundaryVectWest_excluded,
                              idStaggeredBoundaryVectEast_excluded,
                              idHorizontalAll_excluded, horizontalTag_excluded);
    build_merged_ids_and_tags(idStaggeredInternalVectVertical_excluded,
                              idStaggeredBoundaryVectNorth_excluded,
                              idStaggeredBoundaryVectSouth_excluded,
                              idVerticalAll_excluded, verticalTag_excluded);

    // set initial quantities for the time step adaptation
    thrust::device_vector<double> H_old = H_pot;
    thrust::device_vector<double> H_oldold = H_pot;

    thrust::device_vector<double> Gamma_vect_x_1 = Gamma_vect_x_1_pot;
    thrust::device_vector<double> Gamma_vect_x_2 = Gamma_vect_x_2_pot;
    thrust::device_vector<double> Gamma_vect_y_1 = Gamma_vect_y_1_pot;
    thrust::device_vector<double> Gamma_vect_y_2 = Gamma_vect_y_2_pot;

    // build now the slope_x_mod_pot and slope_y_mod_pot to be used in the computeResidualsTruncated function
    const thrust::device_vector<double> slope_x_mod = slope_x_mod_pot;
    const thrust::device_vector<double> slope_y_mod = slope_y_mod_pot;

#else

    // Move to dynamic quantities
    auto H = std::move(H_pot);
    auto u = std::move(u_pot);
    auto v = std::move(v_pot);
    auto additional_source_term = std::move(additional_source_term_pot);
    auto h_G = std::move(h_G_pot); 
    auto h_sd = std::move(h_sd_pot); 
    auto h_sn = std::move(h_sn_pot); 
    auto W_Gav = std::move(W_Gav_pot); 
    auto u_star = std::move(u_star_pot); 
    auto v_star = std::move(v_star_pot); 
    auto h_interface_x = std::move(h_interface_x_pot);
    auto h_interface_y = std::move(h_interface_y_pot); 
    auto eta = std::move(eta_pot); 
  
    // set initial quantities for the time step adaptation
    auto H_old = H;
    auto H_oldold = H;
   
    // Move to const double vectors
    const auto orography = std::move(orography_pot); 
    const auto hydraulic_conductivity = std::move(hydraulic_conductivity_pot);
    const auto Z_Gav = std::move(Z_Gav_pot); 
    const auto n_x = std::move(n_x_pot); 
    const auto n_y = std::move(n_y_pot);
    const auto slope_x = std::move(slope_x_pot); 
    const auto slope_y = std::move(slope_y_pot);
    const auto soilMoistureRetention = std::move(soilMoistureRetention_pot);

    // Move to const ids
    const auto idStaggeredInternalVectHorizontal_excluded = std::move(idStaggeredInternalVectHorizontal_excluded_pot);
    const auto idStaggeredBoundaryVectWest_excluded = std::move(idStaggeredBoundaryVectWest_excluded_pot);
    const auto idStaggeredBoundaryVectEast_excluded = std::move(idStaggeredBoundaryVectEast_excluded_pot);
    const auto idStaggeredInternalVectVertical_excluded = std::move(idStaggeredInternalVectVertical_excluded_pot);
    const auto idStaggeredBoundaryVectNorth_excluded = std::move(idStaggeredBoundaryVectNorth_excluded_pot);
    const auto idStaggeredBoundaryVectSouth_excluded = std::move(idStaggeredBoundaryVectSouth_excluded_pot);
    const auto idBasinVect_excluded = std::move(idBasinVect_excluded_pot);
    const auto idBasinVectReIndex_excluded = std::move(idBasinVectReIndex_excluded_pot);

    const auto idStaggeredInternalVectHorizontal = std::move(idStaggeredInternalVectHorizontal_pot);
    const auto idStaggeredBoundaryVectWest = std::move(idStaggeredBoundaryVectWest_pot);
    const auto idStaggeredBoundaryVectEast = std::move(idStaggeredBoundaryVectEast_pot);
    const auto idStaggeredInternalVectVertical = std::move(idStaggeredInternalVectVertical_pot);
    const auto idStaggeredBoundaryVectNorth = std::move(idStaggeredBoundaryVectNorth_pot);
    const auto idStaggeredBoundaryVectSouth = std::move(idStaggeredBoundaryVectSouth_pot);
    const auto idBasinVect = std::move(idBasinVect_pot);

#endif

    // +-----------------------------------------------+
    // |                    Rain                       |
    // +-----------------------------------------------+

    Rain precipitation(infiltrationModel, N,
                       dataFile("files/infiltration/isInitialLoss", false),
                       dataFile("files/infiltration/perc_initialLoss", 0.05),
                       is_precipitation, constant_precipitation,
                       precipitation_file, file_dir, time_spacing_rain,
                       number_stations, max_Days, dataFile, xllcorner,
                       yllcorner, pixel_size, N_rows, N_cols, 
#ifdef ENABLE_CUDA
		       idBasinVect_pot
#else		       
		       idBasinVect
#endif
		       );

    dt_min = std::min(precipitation.dt_rain, dt_temp);
#ifndef ENABLE_CUDA
    for (const auto &kk : hydraulic_conductivity) {
      if (kk * (dt_min / pixel_size) > 1.) {
        if (rank == 0)
          std::cout << "Error! Courant number for gravitational layer is "
                       "greater than 1!"
                    << std::endl;
        exit(-1.);
      }
    }
#endif

    // +-----------------------------------------------+
    // |                 Temperature                   |
    // +-----------------------------------------------+

    Temperature temp(file_dir + temperature_file, N, max_Days, T_thr, 
                     std::round(max_Days * 24 / time_spacing_temp),
                     steps_per_hour, time_spacing_temp, height_thermometer,
                     format_temp, 
		     orography, 
#ifdef ENABLE_CUDA
		     idBasinVect_excluded
#else
		     idBasinVect
#endif
		     ); // here we pass the device vectors, because are not used during initialization 

    // +-----------------------------------------------+
    // |              Evapotranspiration               |
    // +-----------------------------------------------+

    evapoTranspiration ET(ET_model, N, temp.J, max_Days, phi_rad,
                          height_thermometer, 
#ifdef ENABLE_CUDA
			  idBasinVect_excluded,
#else
			  idBasinVect, 
#endif
			  orography); // here we pass the device vectors, because are not used during initialization

    // +-----------------------------------------------+
    // |                   Runoff                      |
    // +-----------------------------------------------+

    // here we pass the device vectors, because are not used during initialization
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
    
    // +-----------------------------------------------+
    // |             Start the time loop               |
    // +-----------------------------------------------+
   
#ifdef ENABLE_CUDA
    // probably add here some resize for the cuSPARSE objects
    int *d_A_rows, *d_A_columns;
    double *d_A_values, *d_L_values;
    Vec     d_B, d_X, d_R, d_R_aux, d_P, d_T, d_tmp;

    const int m = idBasinVect_excluded.size();
    const int num_offsets = m + 1;

    //--------------------------------------------------------------------------
    // Build the CSR sparsity pattern first and take nnz from the pattern itself.
    // The old closed-form  nnz = 5*m - <domain-boundary faces>  is wrong as soon
    // as the basin contains excluded cells (pour points / steep cells): the
    // star-stencil then writes a different number of entries, the CSR is
    // under-allocated and corrupts, and cuSPARSE aborts in the IC factorization.
    // This bit only at the finest grid (10 m) where excluded cells are common.
    int*    h_A_rows    = (int*)malloc(num_offsets * sizeof(int));
    int*    h_A_columns = (int*)malloc((size_t)5 * m * sizeof(int)); // safe upper bound

    // Red-black reordering of the basin-cell row/column numbering, used only
    // for the sparse system (buildMatrix / IC(0) / CG / updateH). IC(0)'s
    // triangular solve is inherently sequential along the matrix's
    // dependency graph; cuSPARSE's SpSV analysis builds that graph via
    // level-scheduling, and for the natural raster-scan ordering of a
    // 5-point-stencil 2D grid the number of levels scales with the grid
    // dimension. Red-black ordering (all "red" cells i+j even first, then
    // all "black" cells) makes every red cell's neighbors black and vice
    // versa, collapsing the dependency graph to 2 levels regardless of grid
    // size. This is purely a relabeling (same nnz, same physics, same
    // solution) so it is applied GPU-only: CPU keeps the natural ordering
    // fed to Eigen, which reorders internally (AMD) anyway.
    thrust::host_vector<unsigned int> idBasinVect_rb(m);
    thrust::host_vector<unsigned int> idBasinVectReIndex_rb(
        idBasinVectReIndex_excluded_pot.size());
    {
      unsigned int n_red = 0;
      for (const auto Id : idBasinVect_excluded_pot) {
        const unsigned int i = Id / N_cols, j = Id % N_cols;
        if (((i + j) & 1u) == 0) n_red++;
      }
      unsigned int next_red = 0, next_black = n_red;
      for (const auto Id : idBasinVect_excluded_pot) {
        const unsigned int i = Id / N_cols, j = Id % N_cols;
        const unsigned int pos =
            (((i + j) & 1u) == 0) ? next_red++ : next_black++;
        idBasinVect_rb[pos] = Id;
        idBasinVectReIndex_rb[Id] = pos;
      }
    }

    make_sparsity_pattern(idBasinVect_rb, basin_mask_Vec,
                          idBasinVectReIndex_rb,
                          h_A_rows, h_A_columns,
                          N_rows, N_cols);

    // make_sparsity_pattern inserts a row's neighbors in fixed geometric
    // order (N, W, self, E, S). Under the natural reindex that happens to
    // come out ascending, but red-black scrambles it (a cell's neighbors are
    // all the opposite color, with no relation between geometric direction
    // and their reindex value) -- and cuSPARSE's CSR routines (csric02,
    // SpMV, SpSV) require sorted column indices per row. Sort each row here;
    // findPosition's linear search doesn't care about order either way.
    for (int row = 0; row < m; row++) {
      std::sort(h_A_columns + h_A_rows[row], h_A_columns + h_A_rows[row + 1]);
    }

    const int nnz = h_A_rows[m];
    if (rank == 0)
      std::cout << "[sparsity] m = " << m << ", nnz = " << nnz << std::endl;

    // Uploaded once: buildMatrix_cell_gather uses this to tell an internal
    // neighbor (basin-basin) from a boundary one; the mask never changes
    // during the run.
    const thrust::device_vector<unsigned int> basin_mask_gpu(
        basin_mask_Vec.begin(), basin_mask_Vec.end());

    // Device copy of the red-black reindex, used in place of
    // idBasinVectReIndex_excluded for buildMatrix/updateH on the GPU path.
    const thrust::device_vector<unsigned int> idBasinVectReIndex_rb_gpu(
        idBasinVectReIndex_rb.begin(), idBasinVectReIndex_rb.end());

    //--------------------------------------------------------------------------
    // allocate device memory for CSR matrices
    CHECK_CUDA( cudaMalloc((void**) &d_A_rows,    num_offsets * sizeof(int)) )
    CHECK_CUDA( cudaMalloc((void**) &d_A_columns, nnz * sizeof(int)) )
    CHECK_CUDA( cudaMalloc((void**) &d_A_values,  nnz * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_L_values,  nnz * sizeof(double)) )

    CHECK_CUDA( cudaMalloc((void**) &d_B.ptr,     m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_X.ptr,     m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_R.ptr,     m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_R_aux.ptr, m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_P.ptr,     m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_T.ptr,     m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_tmp.ptr,   m * sizeof(double)) )

    // ### cuSPARSE Handle and descriptors initialization ###
    cublasHandle_t   cublasHandle   = NULL;
    cusparseHandle_t cusparseHandle = NULL;
    CHECK_CUBLAS( cublasCreate(&cublasHandle) )
    CHECK_CUSPARSE( cusparseCreate(&cusparseHandle) )
    // Create dense vectors
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_B.vec,     m, d_B.ptr, CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_X.vec,     m, d_X.ptr, CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_R.vec,     m, d_R.ptr, CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_R_aux.vec, m, d_R_aux.ptr, CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_P.vec,     m, d_P.ptr,   CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_T.vec,     m, d_T.ptr,   CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_tmp.vec,   m, d_tmp.ptr, CUDA_R_64F) )

    // copy the CSR matrices and vectors into device memory
    CHECK_CUDA( cudaMemcpy(d_A_rows, h_A_rows, num_offsets * sizeof(int),
                           cudaMemcpyHostToDevice) )
    CHECK_CUDA( cudaMemcpy(d_A_columns, h_A_columns, nnz *  sizeof(int),
                           cudaMemcpyHostToDevice) )

    CHECK_CUDA( cudaMemset(d_A_values, 0, nnz * sizeof(double)) )
    CHECK_CUDA( cudaMemset(d_L_values, 1, nnz * sizeof(double)) )
    CHECK_CUDA( cudaMemset(d_X.ptr,    0, m   * sizeof(double)) )

    // free host memory once uploaded — no longer needed
    free(h_A_rows);
    free(h_A_columns);

    cusparseIndexBase_t  baseIdx = CUSPARSE_INDEX_BASE_ZERO;
    cusparseSpMatDescr_t matA, matL;
    int*                 d_L_rows      = d_A_rows;
    int*                 d_L_columns   = d_A_columns;
    cusparseFillMode_t   fill_lower    = CUSPARSE_FILL_MODE_LOWER;
    cusparseDiagType_t   diag_non_unit = CUSPARSE_DIAG_TYPE_NON_UNIT;
    // A
    CHECK_CUSPARSE( cusparseCreateCsr(&matA, m, m, nnz, d_A_rows,
                                      d_A_columns,
				      d_A_values,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      baseIdx, CUDA_R_64F) )
    // L
    CHECK_CUSPARSE( cusparseCreateCsr(&matL, m, m, nnz, d_L_rows,
                                      d_L_columns, d_L_values,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      baseIdx, CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseSpMatSetAttribute(matL,
                                              CUSPARSE_SPMAT_FILL_MODE,
                                              &fill_lower, sizeof(fill_lower)) )
    CHECK_CUSPARSE( cusparseSpMatSetAttribute(matL,
                                              CUSPARSE_SPMAT_DIAG_TYPE,
                                              &diag_non_unit,
                                              sizeof(diag_non_unit)) )

    //--------------------------------------------------------------------------
    // ### Preparation ### b = A * X
    const double alpha = 1.0;
    size_t       bufferSizeMV;
    void*        d_bufferMV;
    double       beta = 0.0;
    CHECK_CUSPARSE( cusparseSpMV_bufferSize(
                        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha, matA, d_X.vec, &beta, d_B.vec, CUDA_R_64F,
                        CUSPARSE_SPMV_ALG_DEFAULT, &bufferSizeMV) )
    CHECK_CUDA( cudaMalloc(&d_bufferMV, bufferSizeMV) )

    // X0 = 0
    cusparseMatDescr_t descrM;
    csric02Info_t      infoM        = NULL;
    int                bufferSizeIC = 0;
    void*              d_bufferIC;
    CHECK_CUSPARSE( cusparseCreateMatDescr(&descrM) )
    CHECK_CUSPARSE( cusparseSetMatIndexBase(descrM, baseIdx) )
    CHECK_CUSPARSE( cusparseSetMatType(descrM, CUSPARSE_MATRIX_TYPE_GENERAL) )
    CHECK_CUSPARSE( cusparseSetMatFillMode(descrM, CUSPARSE_FILL_MODE_LOWER) )
    CHECK_CUSPARSE( cusparseSetMatDiagType(descrM,
                                           CUSPARSE_DIAG_TYPE_NON_UNIT) )
    CHECK_CUSPARSE( cusparseCreateCsric02Info(&infoM) )

    CHECK_CUSPARSE( cusparseDcsric02_bufferSize(
                        cusparseHandle, m, nnz, descrM, d_L_values,
                        d_A_rows, d_A_columns, infoM, &bufferSizeIC) )
    CHECK_CUDA( cudaMalloc(&d_bufferIC, bufferSizeIC) )


    CHECK_CUSPARSE( cusparseDcsric02_analysis(
			    cusparseHandle, m, nnz, descrM, d_L_values,
			    d_A_rows, d_A_columns, infoM,
			    CUSPARSE_SOLVE_POLICY_NO_LEVEL, d_bufferIC) )
    int structural_zero;
    CHECK_CUSPARSE( cusparseXcsric02_zeroPivot(cusparseHandle, infoM,
			    &structural_zero) )

    cusparseSpSVDescr_t spsvDescrL, spsvDescrLT;
    void *d_bufferL, *d_bufferLT;
    size_t bufferSizeL, bufferSizeLT;
    const double one = 1.0;

    CHECK_CUSPARSE( cusparseSpSV_createDescr(&spsvDescrL) )
    CHECK_CUSPARSE( cusparseSpSV_createDescr(&spsvDescrLT) )

    CHECK_CUSPARSE( cusparseSpSV_bufferSize(
			    cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
			    &one, matL, d_R.vec, d_tmp.vec, CUDA_R_64F,
			    CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL, &bufferSizeL) )
    CHECK_CUSPARSE( cusparseSpSV_bufferSize(
			    cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE,
			    &one, matL, d_tmp.vec, d_R_aux.vec, CUDA_R_64F,
			    CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrLT, &bufferSizeLT) )

    CHECK_CUDA( cudaMalloc(&d_bufferL,  bufferSizeL) )
    CHECK_CUDA( cudaMalloc(&d_bufferLT, bufferSizeLT) )

    // NOTE: cusparseSpSV_analysis is NOT purely structural for
    // CUSPARSE_SPSV_ALG_DEFAULT -- it depends on the actual matL VALUES
    // (which row-dependency/level-scheduling it builds can differ depending
    // on which entries are numerically zero). d_L_values here is still the
    // cudaMemset(...,1,...) placeholder from above, not a real factorization,
    // so analysis is deliberately NOT run yet -- it is redone every time step
    // right after the real cusparseDcsric02 factorization, using the actual L.

    //--------------------------------------------------------------------------
    // Set the NetCDF library
    int ncid;
    int t_dimid, y_dimid, x_dimid, yu_dimid, xu_dimid, yv_dimid, xv_dimid;
    int // coordinate variable ids
	    xc_varid, yc_varid,
	    xu_varid, yu_varid,   // u face coordinates
	    xv_varid, yv_varid,   // v face coordinates
	    t_varid;
    int // field variable ids
	    H_varid, h_sn_varid, h_G_varid, h_sd_varid,
	    u_varid, v_varid;

    nc_create((output_dir + "output.nc").c_str(), NC_NETCDF4, &ncid);

    // dimensions
    nc_def_dim(ncid, "time", NC_UNLIMITED, &t_dimid);
    nc_def_dim(ncid, "y",    N_rows,       &y_dimid);    // cell center rows
    nc_def_dim(ncid, "x",    N_cols,       &x_dimid);    // cell center cols
    nc_def_dim(ncid, "yu",   N_rows,       &yu_dimid);   // u rows  (same as center)
    nc_def_dim(ncid, "xu",   N_cols + 1,   &xu_dimid);   // u cols  (one extra)
    nc_def_dim(ncid, "yv",   N_rows + 1,   &yv_dimid);   // v rows  (one extra)
    nc_def_dim(ncid, "xv",   N_cols,       &xv_dimid);   // v cols  (same as center)

    // time coordinate
    nc_def_var(ncid, "time", NC_DOUBLE, 1, &t_dimid, &t_varid);
    nc_put_att_text(ncid, t_varid, "units",         21, "seconds since 2000-01-01");
    nc_put_att_text(ncid, t_varid, "standard_name",  4, "time");
    nc_put_att_text(ncid, t_varid, "calendar",       9, "gregorian");

    // cell center coordinates
    nc_def_var(ncid, "xc", NC_DOUBLE, 1, &x_dimid,  &xc_varid);
    nc_def_var(ncid, "yc", NC_DOUBLE, 1, &y_dimid,  &yc_varid);
    nc_put_att_text(ncid, xc_varid, "units",         1, "m");
    nc_put_att_text(ncid, xc_varid, "standard_name", 21, "projection_x_coordinate");
    nc_put_att_text(ncid, yc_varid, "units",         1, "m");
    nc_put_att_text(ncid, yc_varid, "standard_name", 21, "projection_y_coordinate");

    // u face coordinates
    nc_def_var(ncid, "xu", NC_DOUBLE, 1, &xu_dimid, &xu_varid);
    nc_def_var(ncid, "yu", NC_DOUBLE, 1, &yu_dimid, &yu_varid);
    nc_put_att_text(ncid, xu_varid, "units",         1, "m");
    nc_put_att_text(ncid, xu_varid, "standard_name", 21, "projection_x_coordinate");
    nc_put_att_text(ncid, yu_varid, "units",         1, "m");
    nc_put_att_text(ncid, yu_varid, "standard_name", 21, "projection_y_coordinate");

    // v face coordinates
    nc_def_var(ncid, "xv", NC_DOUBLE, 1, &xv_dimid, &xv_varid);
    nc_def_var(ncid, "yv", NC_DOUBLE, 1, &yv_dimid, &yv_varid);
    nc_put_att_text(ncid, xv_varid, "units",         1, "m");
    nc_put_att_text(ncid, xv_varid, "standard_name", 21, "projection_x_coordinate");
    nc_put_att_text(ncid, yv_varid, "units",         1, "m");
    nc_put_att_text(ncid, yv_varid, "standard_name", 21, "projection_y_coordinate");

    double fill = -9999.0;

    // ── cell centered variables (t, y, x) ──────────────────
    int dims_center[3] = {t_dimid, y_dimid, x_dimid};

    nc_def_var(ncid, "H", NC_DOUBLE, 3, dims_center, &H_varid);
    nc_put_att_text(ncid, H_varid, "long_name",   11, "water depth");
    nc_put_att_text(ncid, H_varid, "units",        1, "m");
    nc_put_att_text(ncid, H_varid, "coordinates",  5, "yc xc");
    nc_put_att_double(ncid, H_varid, "_FillValue", NC_DOUBLE, 1, &fill);

    nc_def_var(ncid, "h_sn", NC_DOUBLE, 3, dims_center, &h_sn_varid);
    nc_put_att_text(ncid, h_sn_varid, "long_name",   10, "snow depth");
    nc_put_att_text(ncid, h_sn_varid, "units",        1, "m");
    nc_put_att_text(ncid, h_sn_varid, "coordinates",  5, "yc xc");
    nc_put_att_double(ncid, h_sn_varid, "_FillValue", NC_DOUBLE, 1, &fill);

    nc_def_var(ncid, "h_G", NC_DOUBLE, 3, dims_center, &h_G_varid);
    nc_put_att_text(ncid, h_G_varid, "long_name",   12, "ground water");
    nc_put_att_text(ncid, h_G_varid, "units",        1, "m");
    nc_put_att_text(ncid, h_G_varid, "coordinates",  5, "yc xc");
    nc_put_att_double(ncid, h_G_varid, "_FillValue", NC_DOUBLE, 1, &fill);

    nc_def_var(ncid, "h_sd", NC_DOUBLE, 3, dims_center, &h_sd_varid);
    nc_put_att_text(ncid, h_sd_varid, "long_name",   13, "soil moisture");
    nc_put_att_text(ncid, h_sd_varid, "units",        1, "m");
    nc_put_att_text(ncid, h_sd_varid, "coordinates",  5, "yc xc");
    nc_put_att_double(ncid, h_sd_varid, "_FillValue", NC_DOUBLE, 1, &fill);

    // ── staggered variables ─────────────────────────────────
    int dims_u[3] = {t_dimid, yu_dimid, xu_dimid};  // staggered in x
    int dims_v[3] = {t_dimid, yv_dimid, xv_dimid};  // staggered in y

    nc_def_var(ncid, "u", NC_DOUBLE, 3, dims_u, &u_varid);
    nc_put_att_text(ncid, u_varid, "long_name",   10, "x-velocity");
    nc_put_att_text(ncid, u_varid, "units",        3, "m/s");
    nc_put_att_text(ncid, u_varid, "coordinates",  5, "yu xu");
    nc_put_att_double(ncid, u_varid, "_FillValue", NC_DOUBLE, 1, &fill);

    nc_def_var(ncid, "v", NC_DOUBLE, 3, dims_v, &v_varid);
    nc_put_att_text(ncid, v_varid, "long_name",   10, "y-velocity");
    nc_put_att_text(ncid, v_varid, "units",        3, "m/s");
    nc_put_att_text(ncid, v_varid, "coordinates",  5, "yv xv");
    nc_put_att_double(ncid, v_varid, "_FillValue", NC_DOUBLE, 1, &fill);

    // global attributes — add CRS info
    nc_put_att_text(ncid, NC_GLOBAL, "Conventions",  6, "CF-1.8");
    nc_put_att_text(ncid, NC_GLOBAL, "title",        20, "Hydrological simulation");
    nc_put_att_text(ncid, NC_GLOBAL, "crs",          22, "EPSG:32632");  // UTM zone 32N

    nc_enddef(ncid);

    // separate coordinate arrays for each grid
    std::vector<double> xc(N_cols),     yc(N_rows);
    std::vector<double> xu(N_cols + 1), yu(N_rows);
    std::vector<double> xv(N_cols),     yv(N_rows + 1);

    // top-left origins
    const double originY_c = yllcorner             + N_rows       * pixel_size;
    const double originY_u = yllcorner_staggered_u + N_rows       * pixel_size;
    const double originY_v = yllcorner_staggered_v + (N_rows + 1) * pixel_size;

    // cell center coordinates
    for (int i = 0; i < N_cols;     i++) xc[i] = xllcorner             + (i + 0.5) * pixel_size;
    for (int j = 0; j < N_rows;     j++) yc[j] = originY_c             - (j + 0.5) * pixel_size;

    // u face coordinates
    for (int i = 0; i < N_cols + 1; i++) xu[i] = xllcorner_staggered_u + (i + 0.5) * pixel_size;
    for (int j = 0; j < N_rows;     j++) yu[j] = originY_u             - (j + 0.5) * pixel_size;

    // v face coordinates
    for (int i = 0; i < N_cols;     i++) xv[i] = xllcorner_staggered_v + (i + 0.5) * pixel_size;
    for (int j = 0; j < N_rows + 1; j++) yv[j] = originY_v             - (j + 0.5) * pixel_size;

    // write coordinate arrays once
    nc_put_var_double(ncid, xc_varid, xc.data());
    nc_put_var_double(ncid, yc_varid, yc.data());
    nc_put_var_double(ncid, xu_varid, xu.data());
    nc_put_var_double(ncid, yu_varid, yu.data());
    nc_put_var_double(ncid, xv_varid, xv.data());
    nc_put_var_double(ncid, yv_varid, yv.data());

    // helper lambda to avoid duplication ──
    auto saveToNetCDF = [&](size_t t_idx, double t_val) {
	    // copy all fields to host
	    H_pot    = H;
	    h_sn_pot = h_sn;
	    h_G_pot  = h_G;
	    h_sd_pot = h_sd;
	    u_pot    = u;
	    v_pot    = v;

	    // write time value
	    nc_put_var1_double(ncid, t_varid, &t_idx, &t_val);

	    // write fields
	    size_t start[3]   = {t_idx, 0, 0};
	    size_t count_c[3] = {1, (size_t)N_rows,     (size_t)N_cols    };
	    size_t count_u[3] = {1, (size_t)N_rows,     (size_t)N_cols + 1};
	    size_t count_v[3] = {1, (size_t)N_rows + 1, (size_t)N_cols    };

	    nc_put_vara_double(ncid, H_varid,    start, count_c, H_pot.data());
	    nc_put_vara_double(ncid, h_sn_varid, start, count_c, h_sn_pot.data());
	    nc_put_vara_double(ncid, h_G_varid,  start, count_c, h_G_pot.data());
	    nc_put_vara_double(ncid, h_sd_varid, start, count_c, h_sd_pot.data());
	    nc_put_vara_double(ncid, u_varid,    start, count_u, u_pot.data());
	    nc_put_vara_double(ncid, v_varid,    start, count_v, v_pot.data());
    };

#else
    A.resize(idBasinVect_excluded.size(), idBasinVect_excluded.size());
    A.setZero();
    H_basin.setZero(idBasinVect_excluded.size()); 
    rhs.setZero(idBasinVect_excluded.size());
#endif

#ifdef ENABLE_CUDA
    maxH = deviceMax(H);
#else
    maxH = *std::max_element(H.begin(), H.end());
#endif

    c1_min = dt_min / pixel_size;
    area = pixel_size*pixel_size * 1.e-6; // km^2

    dt_DSV = maxdt(u, v, maxH, pixel_size, gravity);
    dt_DSV = dt_DSV < dt_DSV_given ? dt_DSV : dt_DSV_given;
    dt_DSV = dt_DSV < t_final ? dt_DSV : t_final;

    c1_DSV_ = c1_DSV(dt_DSV, pixel_size);
    c2_DSV_ = c2_DSV(c1_DSV_);
    c3_DSV_ = c3_DSV(c1_DSV_);

    time = 0.; 
    timed = -dt_DSV;
    timedd = -2. * dt_DSV;

    // Up to this point, no computations are carried out on the device.
    // Just standard copies to Device are done (which are blocking)

#ifdef ENABLE_CUDA
    // create events once ──────────
    cudaEvent_t e_ET = NULL, e_temperature = NULL, 
		e_precip = NULL, e_H_interface_x = NULL, e_alfa_x = NULL,
		e_H_interface_y = NULL, e_alfa_y = NULL;
    cudaEventCreate(&e_ET);
    cudaEventCreate(&e_temperature);
    cudaEventCreate(&e_precip);
    cudaEventCreate(&e_H_interface_x);
    cudaEventCreate(&e_alfa_x);
    cudaEventCreate(&e_H_interface_y);
    cudaEventCreate(&e_alfa_y);
#endif

    const auto bench_t0 = std::chrono::steady_clock::now();

    while (!is_last_step) {

      ++bench_loop_iters;

      if (rank == 0) {
        std::cout << "Simulation progress: " << time / t_final * 100 << " %"
                  << " max surface run-off vel. based Courant: "
                  << maxCourant(u, v, c1_DSV_)
                  << " max surface run-off cel. based Courant: "
                  << maxCourant(H, c1_DSV_, gravity) << std::endl;
        std::cout << "Current dt, " << dt_DSV << ", given dt, " << dt_DSV_given
                  << std::endl;
      }

      //TODO: check if the host/device calls below are all independent!
      //TODO: check how many SM do we need to allocate!
     
      // Compute interface fluxes via upwind method
      H_interface.computeHorizontal(CUDA_STREAM_ONLY(S(0)));
      H_interface.computeVertical(CUDA_STREAM_ONLY(S(1)));
#ifdef ENABLE_CUDA
      // record when stream finishes
      cudaEventRecord(e_H_interface_x, S(0));
      cudaEventRecord(e_H_interface_y, S(1));
#endif

      // Compute alfa coefficients
      alfa.f_x(CUDA_STREAM_ONLY(S(2)));
      alfa.f_y(CUDA_STREAM_ONLY(S(3)));
#ifdef ENABLE_CUDA
      // record when stream finishes
      cudaEventRecord(e_alfa_x, S(2));
      cudaEventRecord(e_alfa_y, S(3));
#endif

      // Update the temperature 
      temp.computeTemperature(time, dt_temp, dt_DSV CUDA_STREAM(S(4)));
#ifdef ENABLE_CUDA
      // record when stream finishes
      cudaEventRecord(e_temperature, S(4));
#endif

      // ET varies daily
      if (std::floor(time / (24. * 3600)) >
          std::floor((time - dt_DSV) / (24. * 3600))) {

        ET.ET(temp.T_dailyMean, temp.T_dailyMin, temp.T_dailyMax,
              std::floor(time / (24 * 3600)) CUDA_STREAM(S(5)));

      }
#ifdef ENABLE_CUDA
	// record when stream finishes
	cudaEventRecord(e_ET, S(5));
#endif

      // +-----------------------------------------------+
      // |    Gravitational Layer, Snow Accumulation     |
      // +-----------------------------------------------+

      precipitation.computePrecipitation(time, soilMoistureRetention,
                                         temp.melt_mask, h_G, H, N_cols,
                                         idBasinVect_excluded CUDA_STREAM(S(6)));

      // update only if necessary
      if (std::floor(time / dt_min) > std::floor((time - dt_DSV) / dt_min)) {

#ifdef ENABLE_CUDA
	cudaStreamWaitEvent(S(6), e_temperature, 0);
	cudaStreamWaitEvent(S(6), e_ET, 0);
#endif

        // vertical and horizontal fluxes for Gravitational Layer
        computeResiduals(
            n_x, n_y, N_cols, hydraulic_conductivity,
#ifdef ENABLE_CUDA
	    idHorizontalAll_excluded, horizontalTag_excluded,
	    idVerticalAll_excluded, verticalTag_excluded,
	    idBasinVect_excluded,
#else
            idStaggeredInternalVectHorizontal, idStaggeredInternalVectVertical,
            idStaggeredBoundaryVectWest, idStaggeredBoundaryVectEast,
            idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectSouth,
	    idBasinVect,
#endif
	    temp.T_raster, temp.melt_mask,
	    h_sn, h_G, ET.ET_vec,
            precipitation.DP_infiltrated,
	    precipitation.DP_total,
	    c1_min,
	    dt_min, 
	    T_thr, 
            h_interface_x, h_interface_y CUDA_STREAM(S(6)));

      }
#ifdef ENABLE_CUDA
	// record when stream finishes
	cudaEventRecord(e_precip, S(6));
#endif

      // +-----------------------------------------------+
      // |               de Saint Venant                 |
      // +-----------------------------------------------+

      //  fill u_star and v_star with a Bilinear Interpolation
      bilinearInterpolation(u, v, u_star, v_star, N_rows, N_cols, dt_DSV,
                            pixel_size,
#ifdef ENABLE_CUDA
                            idHorizontalAll_excluded, horizontalTag_excluded,
                            idVerticalAll_excluded, verticalTag_excluded
#else
                            idStaggeredInternalVectHorizontal_excluded,
                            idStaggeredInternalVectVertical_excluded,
                            idStaggeredBoundaryVectWest_excluded,
                            idStaggeredBoundaryVectEast_excluded,
                            idStaggeredBoundaryVectNorth_excluded,
                            idStaggeredBoundaryVectSouth_excluded
#endif
                            CUDA_STREAM(S(7)));

#ifndef ENABLE_CUDA
      coefficients.reserve(
          idBasinVect_excluded.size() +
          4 * idStaggeredInternalVectHorizontal_excluded.size() +
          4 * idStaggeredInternalVectVertical_excluded.size());

      additional_source_term.assign(N, 0.);
#endif

      //TODO: insert a check that we do not deal with excluded domains in CUDA
#ifdef ENABLE_CUDA
	cudaStreamWaitEvent(S(7), e_alfa_x, 0);
	cudaStreamWaitEvent(S(7), e_alfa_y, 0);
	cudaStreamWaitEvent(S(7), e_precip, 0);
	cudaStreamWaitEvent(S(7), e_H_interface_x, 0);
	cudaStreamWaitEvent(S(7), e_H_interface_y, 0);
#endif

      buildMatrix(H_interface.horizontal, H_interface.vertical, orography,
                  u_star, v_star, u, v, H, N_cols, c1_DSV_, c3_DSV_,
                  0, precipitation.DP_cumulative, dt_DSV, alfa.alfa_x,
                  alfa.alfa_y,
#ifdef ENABLE_CUDA
                  idHorizontalAll_excluded, horizontalTag_excluded,
                  idVerticalAll_excluded, verticalTag_excluded,
                  idBasinVect_excluded,
#else
                  idStaggeredInternalVectHorizontal_excluded,
                  idStaggeredInternalVectVertical_excluded,
                  idStaggeredBoundaryVectWest_excluded,
                  idStaggeredBoundaryVectEast_excluded,
                  idStaggeredBoundaryVectNorth_excluded,
                  idStaggeredBoundaryVectSouth_excluded, idBasinVect_excluded,
                  idBasinVect, idStaggeredInternalVectHorizontal,
                  idStaggeredInternalVectVertical,
#endif
#ifdef ENABLE_CUDA
                  idBasinVectReIndex_rb_gpu,
#else
                  idBasinVectReIndex_excluded,
#endif
                  isNonReflectingBC, true,
#ifdef ENABLE_CUDA
		  basin_mask_gpu, d_A_rows, d_A_columns, d_A_values, d_B, nnz
#endif
		  CUDA_STREAM(S(7))
#ifndef ENABLE_CUDA
		  additional_source_term, excluded_ids, coefficients, rhs
#endif
		  );
/*
      // +-----------------------------------------------+
      // |   TEMP: buildMatrix validation (CPU vs GPU)    |
      // |   fingerprint the assembled A / rhs before the |
      // |   (still-disabled) linear solve below          |
      // +-----------------------------------------------+
      {
        auto fp_report = [&](const char *nm, long n, long double s,
                              long double s2, long double mx) {
          if (rank == 0)
            printf("[MATRIX-FP] %-3s n=%ld  sum=%+.10Le  l2=%.10Le  max=%.10Le\n",
                   nm, n, s, std::sqrt(s2), mx);
        };

#ifndef ENABLE_CUDA
        A.setFromTriplets(coefficients.begin(), coefficients.end());
        A.makeCompressed();
        coefficients.clear();

        {
          const double *vals = A.valuePtr();
          const long nnz_check = A.nonZeros();
          long double s = 0.0L, s2 = 0.0L, mx = 0.0L;
          for (long kk = 0; kk < nnz_check; kk++) {
            const long double x = static_cast<long double>(vals[kk]);
            const long double a = x < 0.0L ? -x : x;
            s += x; s2 += x * x; if (a > mx) mx = a;
          }
          fp_report("A", nnz_check, s, s2, mx);
        }
        {
          const long n_rhs = static_cast<long>(rhs.size());
          long double s = 0.0L, s2 = 0.0L, mx = 0.0L;
          for (long kk = 0; kk < n_rhs; kk++) {
            const long double x = static_cast<long double>(rhs(kk));
            const long double a = x < 0.0L ? -x : x;
            s += x; s2 += x * x; if (a > mx) mx = a;
          }
          fp_report("rhs", n_rhs, s, s2, mx);
        }
#else
        CHECK_CUDA( cudaStreamSynchronize(S(7)) )
        std::vector<double> h_A_values_check(nnz), h_rhs_check(m);
        CHECK_CUDA( cudaMemcpy(h_A_values_check.data(), d_A_values,
                                nnz * sizeof(double), cudaMemcpyDeviceToHost) )
        CHECK_CUDA( cudaMemcpy(h_rhs_check.data(), d_B.ptr,
                                m * sizeof(double), cudaMemcpyDeviceToHost) )

        {
          long double s = 0.0L, s2 = 0.0L, mx = 0.0L;
          for (int kk = 0; kk < nnz; kk++) {
            const long double x = static_cast<long double>(h_A_values_check[kk]);
            const long double a = x < 0.0L ? -x : x;
            s += x; s2 += x * x; if (a > mx) mx = a;
          }
          fp_report("A", nnz, s, s2, mx);
        }
        {
          long double s = 0.0L, s2 = 0.0L, mx = 0.0L;
          for (int kk = 0; kk < m; kk++) {
            const long double x = static_cast<long double>(h_rhs_check[kk]);
            const long double a = x < 0.0L ? -x : x;
            s += x; s2 += x * x; if (a > mx) mx = a;
          }
          fp_report("rhs", m, s, s2, mx);
        }
#endif
        if (rank == 0) fflush(stdout);
      }
*/

#ifndef ENABLE_CUDA

      A.setFromTriplets(coefficients.begin(), coefficients.end());
      A.makeCompressed();
      coefficients.clear();

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

      // use swap here!
      H_oldold.swap(H_old);
      H_old.swap(H);
      for (const auto &Id : idBasinVect_excluded) {
        const auto IDreIndex = idBasinVectReIndex_excluded[Id];
        H[Id] = std::abs(H_basin(IDreIndex));
        eta[Id] = H[Id] + orography[Id];
      }

#else

      // reset L = A for IC factorization
      CHECK_CUDA( cudaMemcpy(d_L_values, d_A_values,
			      nnz * sizeof(double),
			      cudaMemcpyDeviceToDevice) )

      //--------------------------------------------------------------------------
      // Warm start: d_X is intentionally NOT reset to zero here. It carries
      // over the previous solve's result (zeroed once at setup) as the
      // initial guess for CG, exactly like the CPU reference's H_basin,
      // which persists across steps/retries -- H changes smoothly enough
      // step to step that this sharply cuts the iteration count needed
      // (matches the CPU behavior we already see: 420/61/8/2/1 CG iterations
      // across successive adaptive-dt retries within one accepted step).
      // M = L * L^T
      CHECK_CUSPARSE( cusparseDcsric02(
			      cusparseHandle, m, nnz, descrM, d_L_values,
			      d_A_rows, d_A_columns, infoM,
			      CUSPARSE_SOLVE_POLICY_NO_LEVEL, d_bufferIC) )
      // Find numerical zero. NOTE: cusparseXcsric02_zeroPivot's own return
      // status is about whether the *query* succeeded (always SUCCESS in
      // practice) -- CHECK_CUSPARSE does NOT catch a zero/negative pivot,
      // that information is only in the output value, which must be checked
      // explicitly.
      int numerical_zero = -1;
      CHECK_CUSPARSE( cusparseXcsric02_zeroPivot(cusparseHandle, infoM,
			      &numerical_zero) )
      if (numerical_zero >= 0 && rank == 0)
        std::cout << "[CG] WARNING: IC(0) factorization hit a non-positive "
                     "pivot at reindexed row " << numerical_zero
                  << " -- preconditioner is invalid this step." << std::endl;

      // Re-analyze the triangular solves against the REAL, just-factorized L
      // (cusparseSpSV_analysis is not purely structural for ALG_DEFAULT; A/L
      // change every step, so a one-time analysis against a placeholder is
      // not valid here).
      CHECK_CUSPARSE( cusparseSpSV_analysis(
			      cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
			      &one, matL, d_R.vec, d_tmp.vec, CUDA_R_64F,
			      CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL, d_bufferL) )
      CHECK_CUSPARSE( cusparseSpSV_analysis(
			      cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE,
			      &one, matL, d_tmp.vec, d_R_aux.vec, CUDA_R_64F,
			      CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrLT, d_bufferLT) )

      //--------------------------------------------------------------------------
      // ### Run CG computation ###
      printf("CG loop:\n"); //TODO: check the stream used here after, because it should depend on the stream used in buildMatrix
      gpu_CG(cublasHandle, cusparseHandle, m,
		      matA, matL, d_B, d_X, d_R, d_R_aux, d_P, d_T,
		      d_tmp, d_bufferMV, spsvDescrL, spsvDescrLT,
		      1000,  // max iterations
		      1e-6,  // tolerance
		      use_preconditioner);
      //--------------------------------------------------------------------------

      minH = deviceMin(d_X.ptr, m);
      maxH = deviceMax(d_X.ptr, m);

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

      // Preserve previous H for adaptive time-step estimator
      H_oldold.swap(H_old);
      H_old.swap(H);

      // Scatter solver result into H and update eta = H + orography
      updateH_wrapper(idBasinVect_excluded, idBasinVectReIndex_rb_gpu,
                      H, eta, orography, d_X.ptr, S(7));

#endif

      // +-----------------------------------------------+
      // |                  Update u, v                  |
      // +-----------------------------------------------+

      updateVel(u, v, u_star, v_star, alfa.alfa_x, alfa.alfa_y, N_rows, N_cols,
                c2_DSV_, 0, eta, H, orography,
#ifndef ENABLE_CUDA
                idStaggeredInternalVectHorizontal_excluded,
                idStaggeredInternalVectVertical_excluded,
                idStaggeredBoundaryVectWest_excluded,
                idStaggeredBoundaryVectEast_excluded,
                idStaggeredBoundaryVectNorth_excluded,
                idStaggeredBoundaryVectSouth_excluded,
#else
                idHorizontalAll_excluded, horizontalTag_excluded,
                idVerticalAll_excluded, verticalTag_excluded,
#endif
                isNonReflectingBC CUDA_STREAM(S(7)));

      // +-----------------------------------------------+
      // |        Debug field stats (validation)         |
      // +-----------------------------------------------+
      if (spit_out_solutions_each_time_step && rank == 0) {
#ifdef ENABLE_CUDA
        cudaDeviceSynchronize();
        thrust::host_vector<double> _hG_h(h_G.begin(), h_G.end());
        thrust::host_vector<double> _hsn_h(h_sn.begin(), h_sn.end());
        thrust::host_vector<double> _u_h(u.begin(), u.end());
        thrust::host_vector<double> _v_h(v.begin(), v.end());
#else
        const auto& _hG_h  = h_G;
        const auto& _hsn_h = h_sn;
        const auto& _u_h   = u;
        const auto& _v_h   = v;
#endif
        auto field_stats = [&](const auto& f, const auto& ids) {
          double mx = 0., l2 = 0.;
          for (auto id : ids) { double v = f[id]; if (std::abs(v)>mx) mx=std::abs(v); l2+=v*v; }
          return std::make_pair(mx, std::sqrt(l2));
        };
        auto [hG_mx, hG_l2]   = field_stats(_hG_h,  idBasinVect_excluded);
        auto [hsn_mx, hsn_l2] = field_stats(_hsn_h, idBasinVect_excluded);
        auto [u_mx,   u_l2]   = field_stats(_u_h,   idStaggeredInternalVectHorizontal_excluded);
        auto [v_mx,   v_l2]   = field_stats(_v_h,   idStaggeredInternalVectVertical_excluded);
        printf("[STATS t=%.4f] h_G: max=%.10e l2=%.10e | h_sn: max=%.10e l2=%.10e | u: max=%.10e l2=%.10e | v: max=%.10e l2=%.10e\n",
               time, hG_mx, hG_l2, hsn_mx, hsn_l2, u_mx, u_l2, v_mx, v_l2);
        fflush(stdout);
      }

      // +-----------------------------------------------+
      // |             Sediment Transport                |
      // +-----------------------------------------------+
/*
      if (is_sediment_transport) {
        // tic();

        dt_sed = compute_dt_sediment(alfa_coeff, beta_coeff, slope_x_max,
                                     slope_y_max, u, v, pixel_size, dt_DSV,
                                     &numberOfSteps);
        c1_sed = dt_sed / pixel_size;

        // vertical and horizontal residuals truncated for Sediment Transport
        computeResidualsTruncated(
            u, v, N_cols, N_rows, N, (dt_sed / pixel_size), slope_x, slope_y,
            alfa_coeff,  // alfa
            beta_coeff,  // beta
            gamma_coeff, // gamma
#ifdef ENABLE_CUDA
	    idHorizontalAll_excluded, idVerticalAll_excluded,
	    Gamma_vect_x_1, Gamma_vect_x_2, Gamma_vect_y_1, Gamma_vect_y_2,
	    slope_x_mod, slope_y_mod CUDA_STREAM(S(11))
#else
            idStaggeredInternalVectHorizontal_excluded,
            idStaggeredInternalVectVertical_excluded,
            idStaggeredBoundaryVectWest_excluded,
            idStaggeredBoundaryVectEast_excluded,
            idStaggeredBoundaryVectNorth_excluded,
            idStaggeredBoundaryVectSouth_excluded,
	    Gamma_vect_x, Gamma_vect_y
#endif
	    );

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

        for (const auto &k : idBasinVect_excluded) {
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
          for (const auto &Id :
               idStaggeredInternalVectHorizontal_excluded) {
            const unsigned int i = Id / (N_cols + 1), IDeast = Id - i,
                               IDwest = Id - i - 1;

            const double &h_left = h_sd[IDwest], &h_right = h_sd[IDeast];

            h_interface_x[Id] =
                Gamma_vect_x[Id][0] * h_right + Gamma_vect_x[Id][1] * h_left;
          }

          for (const auto &Id : idStaggeredBoundaryVectWest_excluded) {
            const unsigned int i = Id / (N_cols + 1);

            const double h_left = 0, // 0,
                h_right = h_sd[Id - i];

            h_interface_x[Id] =
                Gamma_vect_x[Id][0] * h_right + Gamma_vect_x[Id][1] * h_left;
          }

          for (const auto &Id : idStaggeredBoundaryVectEast_excluded) {

            const unsigned int i = Id / (N_cols + 1);

            const double h_left = h_sd[Id - i - 1], h_right = 0;

            h_interface_x[Id] =
                Gamma_vect_x[Id][0] * h_right + Gamma_vect_x[Id][1] * h_left;
          }


          // assemble vertical fluxes
          for (const auto &Id :
               idStaggeredInternalVectVertical_excluded) {
            const unsigned int IDsouth = Id, IDnorth = Id - N_cols;

            const double &h_left = h_sd[IDnorth], &h_right = h_sd[IDsouth];

            h_interface_y[Id] =
                Gamma_vect_y[Id][0] * h_right + Gamma_vect_y[Id][1] * h_left;
          }

          for (const auto &Id : idStaggeredBoundaryVectNorth_excluded) {

            const double h_left = 0, 
                h_right = h_sd[Id];

            h_interface_y[Id] =
                Gamma_vect_y[Id][0] * h_right + Gamma_vect_y[Id][1] * h_left;
          }

          for (const auto &Id : idStaggeredBoundaryVectSouth_excluded) {

            const double h_left = h_sd[Id - N_cols], h_right = 0;

            h_interface_y[Id] =
                Gamma_vect_y[Id][0] * h_right + Gamma_vect_y[Id][1] * h_left;
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
          for (const auto &Id : idBasinVect_excluded) {
            const unsigned int i = Id / N_cols;
	    const auto Res_x_cell = h_interface_x[Id + 1 + i] - h_interface_x[Id + i];
	    const auto Res_y_cell = h_interface_y[Id + N_cols] - h_interface_y[Id]; 	    
            h_sd[Id] += -(Res_x_cell + Res_y_cell) + W_Gav[Id] +
                        additional_source_term[Id];
          }
        }

      } // End of if(sediment_transport)
*/
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

#ifndef ENABLE_CUDA
      if (save_temporal_sequence) {

	// Copy back to Host now

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
          for (const auto &candidate : kk_gauges[number - 1]) {
            const unsigned int i = candidate / N_cols;
            const auto &cc = H[candidate];

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
      for (const auto &k : idBasinVect) {
        W_Gav_cum_pot[k] += W_Gav[k] * (dt_DSV / dt_sed);
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
                     W_Gav_cum_pot);
        saveSolution(output_dir + "hG_", " ", N_rows, N_cols, xllcorner,
                     yllcorner, pixel_size, NODATA_value, iter, u, v, h_G);
        saveSolution(output_dir + "hsn_", " ", N_rows, N_cols, xllcorner,
                     yllcorner, pixel_size, NODATA_value, iter, u, v, h_sn);

      } else {

        if (std::floor(time / (frequency_save * 3600)) >
            std::floor((time - dt_DSV) / (frequency_save * 3600))) {
          iter++;

          const auto &currentDay =
              iter; 

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
                       W_Gav_cum_pot);
          saveSolution(output_dir + "hG_", " ", N_rows, N_cols, xllcorner,
                       yllcorner, pixel_size, NODATA_value, currentDay, u, v,
                       h_G);
          saveSolution(output_dir + "hsn_", " ", N_rows, N_cols, xllcorner,
                       yllcorner, pixel_size, NODATA_value, currentDay, u, v,
                       h_sn);
        }
      }
#else

      if (spit_out_solutions_each_time_step) {

	      saveToNetCDF(iter++, time);

      } else {

	      if (std::floor(time / (frequency_save * 3600)) >
			      std::floor((time - dt_DSV) / (frequency_save * 3600))) {

		      if (rank == 0)
			      std::cout << "Saving solution ..., current saving "
				      << iter << std::endl;

		      saveToNetCDF(iter++, time);
	      }
      }


#endif
      // +-----------------------------------------------+
      // |               Update time step                |
      // +-----------------------------------------------+

      dt_DSV = maxdt(u, v, maxH, pixel_size, gravity);
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

      // +-----------------------------------------------+
      // |   BENCH: count accepted step, maybe stop      |
      // +-----------------------------------------------+
      ++bench_completed_steps;
      if (stop_after_first_step ||
          (max_bench_steps > 0 && bench_completed_steps >= max_bench_steps)) {
        if (rank == 0)
          std::cout << "[BENCH] stopping after " << bench_completed_steps
                    << " accepted step(s)" << std::endl;
        break;
      }

    } // End Time Loop

    // +-------------------------------------------------------------+
    // |   BENCH report: wall-clock of the time loop + field        |
    // |   fingerprints (compare CPU vs CUDA run for validation)    |
    // +-------------------------------------------------------------+
    {
#ifdef ENABLE_CUDA
      cudaDeviceSynchronize();
#endif
      const double bench_wall_s =
          std::chrono::duration<double>(
              std::chrono::steady_clock::now() - bench_t0)
              .count();

      auto bench_fp = [&](const char *nm, const auto &fld, const auto &ids) {
#ifdef ENABLE_CUDA
        const thrust::host_vector<double> h(fld.begin(), fld.end());
        const thrust::host_vector<unsigned int> ih(ids.begin(), ids.end());
#else
        const auto &h = fld;
        const auto &ih = ids;
#endif
        long double s = 0.0L, s2 = 0.0L, mx = 0.0L;
        for (auto id : ih) {
          const long double x = static_cast<long double>(h[id]);
          const long double a = x < 0.0L ? -x : x;
          s += x;
          s2 += x * x;
          if (a > mx)
            mx = a;
        }
        if (rank == 0)
          printf("[FP] %-4s n=%zu  sum=%+.10Le  l2=%.10Le  max=%.10Le\n", nm,
                 static_cast<size_t>(ih.size()), s, std::sqrt(s2), mx);
      };

      if (rank == 0) {
        printf("[BENCH] backend=%s  accepted_steps=%u  loop_iters=%u  "
               "wall_s=%.6f  per_step_ms=%.4f  per_iter_ms=%.4f\n",
#ifdef ENABLE_CUDA
               "CUDA",
#else
               "CPU",
#endif
               bench_completed_steps, bench_loop_iters, bench_wall_s,
               bench_completed_steps
                   ? 1.0e3 * bench_wall_s / bench_completed_steps
                   : 0.0,
               bench_loop_iters ? 1.0e3 * bench_wall_s / bench_loop_iters : 0.0);
      }
      bench_fp("H", H, idBasinVect_excluded);
      bench_fp("eta", eta, idBasinVect_excluded);
      bench_fp("hG", h_G, idBasinVect_excluded);
      bench_fp("hsn", h_sn, idBasinVect_excluded);
      bench_fp("hsd", h_sd, idBasinVect_excluded);
      bench_fp("u", u, idStaggeredInternalVectHorizontal_excluded);
      bench_fp("v", v, idStaggeredInternalVectVertical_excluded);
      if (rank == 0)
        fflush(stdout);
    }

#ifdef ENABLE_CUDA    

    // cleanup ──────────────────────
    cudaEventDestroy(e_ET);
    cudaEventDestroy(e_temperature);
    cudaEventDestroy(e_precip);
    cudaEventDestroy(e_H_interface_x);
    cudaEventDestroy(e_alfa_x);
    cudaEventDestroy(e_H_interface_y);
    cudaEventDestroy(e_alfa_y);

    // ### Free resources ###
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_B.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_X.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_R.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_R_aux.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_P.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_T.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_tmp.vec) )
    CHECK_CUSPARSE( cusparseDestroySpMat(matA) )
    CHECK_CUSPARSE( cusparseDestroySpMat(matL) )
    CHECK_CUSPARSE( cusparseDestroy(cusparseHandle) )
    CHECK_CUBLAS( cublasDestroy(cublasHandle) )
    CHECK_CUSPARSE( cusparseDestroyCsric02Info(infoM) )
    CHECK_CUSPARSE( cusparseDestroyMatDescr(descrM) )
    CHECK_CUDA( cudaFree(d_bufferIC) )
    CHECK_CUDA( cudaFree(d_bufferMV) )

    CHECK_CUSPARSE( cusparseSpSV_destroyDescr(spsvDescrL) )
    CHECK_CUSPARSE( cusparseSpSV_destroyDescr(spsvDescrLT) )
    CHECK_CUDA( cudaFree(d_bufferL) )
    CHECK_CUDA( cudaFree(d_bufferLT) )

    // add to cleanup:
    CHECK_CUDA( cudaFree(d_A_rows) )
    CHECK_CUDA( cudaFree(d_A_columns) )
    CHECK_CUDA( cudaFree(d_A_values) )
    CHECK_CUDA( cudaFree(d_L_values) )
    CHECK_CUDA( cudaFree(d_B.ptr) )
    CHECK_CUDA( cudaFree(d_X.ptr) )
    CHECK_CUDA( cudaFree(d_R.ptr) )
    CHECK_CUDA( cudaFree(d_R_aux.ptr) )
    CHECK_CUDA( cudaFree(d_P.ptr) )
    CHECK_CUDA( cudaFree(d_T.ptr) )
    CHECK_CUDA( cudaFree(d_tmp.ptr) )
#endif

  } // End Monte Carlo loop

#ifdef ENABLE_CUDA
  for (auto &s : streams)
    cudaStreamDestroy(s);
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    print_timing_report();
  }
  MPI_Finalize();

  return EXIT_SUCCESS;
}
