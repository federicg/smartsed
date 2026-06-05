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

  const int number_stations = dataFile("files/meteo_data/number_stations", 1);

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

    unsigned int N_rows, N_cols, N, numberOfSteps = 1, iter = 0;

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

    constexpr double alfa_coeff = 2.5, beta_coeff = 1.6, gamma_coeff = 1.;
    static constexpr double gravity = 9.81;
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

        orography_pot, h_G_pot, h_sd_pot, h_sn_pot, W_Gav_pot, W_Gav_cum_pot,
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

    thrust::device_vector<double> H = H_pot;
    thrust::device_vector<double> u = u_pot;
    thrust::device_vector<double> v = v_pot;
    thrust::device_vector<double> additional_source_term = additional_source_term_pot;
    thrust::device_vector<double> orography = orography_pot;
    thrust::device_vector<double> h_G = h_G_pot;
    thrust::device_vector<double> h_sd = h_sd_pot;
    thrust::device_vector<double> h_sn = h_sn_pot;
    thrust::device_vector<double> W_Gav = W_Gav_pot;
    thrust::device_vector<double> W_Gav_cum = W_Gav_cum_pot;
    thrust::device_vector<double> hydraulic_conductivity = hydraulic_conductivity_pot;
    thrust::device_vector<double> Z_Gav = Z_Gav_pot;
    thrust::device_vector<double> n_x = n_x_pot;
    thrust::device_vector<double> n_y = n_y_pot;
    thrust::device_vector<double> u_star = u_star_pot;
    thrust::device_vector<double> v_star = v_star_pot;
    thrust::device_vector<double> h_interface_x = h_interface_x_pot;
    thrust::device_vector<double> h_interface_y = h_interface_y_pot;
    thrust::device_vector<double> slope_x = slope_x_pot;
    thrust::device_vector<double> slope_y = slope_y_pot;
    thrust::device_vector<double> soilMoistureRetention = soilMoistureRetention_pot;
    thrust::device_vector<double> eta = eta_pot;

    thrust::device_vector<unsigned int> idStaggeredInternalVectHorizontal_excluded = idStaggeredInternalVectHorizontal_excluded_pot;
    thrust::device_vector<unsigned int> idStaggeredBoundaryVectWest_excluded = idStaggeredBoundaryVectWest_excluded_pot;
    thrust::device_vector<unsigned int> idStaggeredBoundaryVectEast_excluded = idStaggeredBoundaryVectEast_excluded_pot;
    thrust::device_vector<unsigned int> idStaggeredInternalVectVertical_excluded = idStaggeredInternalVectVertical_excluded_pot;
    thrust::device_vector<unsigned int> idStaggeredBoundaryVectNorth_excluded = idStaggeredBoundaryVectNorth_excluded_pot;
    thrust::device_vector<unsigned int> idStaggeredBoundaryVectSouth_excluded = idStaggeredBoundaryVectSouth_excluded_pot;
    thrust::device_vector<unsigned int> idBasinVect_excluded = idBasinVect_excluded_pot;
    thrust::device_vector<unsigned int> idBasinVectReIndex_excluded = idBasinVectReIndex_excluded_pot;

    thrust::device_vector<unsigned int> idStaggeredInternalVectHorizontal = idStaggeredInternalVectHorizontal_pot;
    thrust::device_vector<unsigned int> idStaggeredBoundaryVectWest = idStaggeredBoundaryVectWest_pot;
    thrust::device_vector<unsigned int> idStaggeredBoundaryVectEast = idStaggeredBoundaryVectEast_pot;
    thrust::device_vector<unsigned int> idStaggeredInternalVectVertical = idStaggeredInternalVectVertical_pot;
    thrust::device_vector<unsigned int> idStaggeredBoundaryVectNorth = idStaggeredBoundaryVectNorth_pot;
    thrust::device_vector<unsigned int> idStaggeredBoundaryVectSouth = idStaggeredBoundaryVectSouth_pot;
    thrust::device_vector<unsigned int> idBasinVect = idBasinVect_pot;
 
    // set initial quantities for the time step adaptation
    thrust::device_vector<double> H_old = H_pot;
    thrust::device_vector<double> H_oldold = H_pot;

    thrust::device_vector<double> Gamma_vect_x_1 = Gamma_vect_x_1_pot;
    thrust::device_vector<double> Gamma_vect_x_2 = Gamma_vect_x_2_pot;
    thrust::device_vector<double> Gamma_vect_y_1 = Gamma_vect_y_1_pot;
    thrust::device_vector<double> Gamma_vect_y_2 = Gamma_vect_y_2_pot;

    // build now the slope_x_mod_pot and slope_y_mod_pot to be used in the computeResidualsTruncated function
    thrust::device_vector<double> slope_x_mod = slope_x_mod_pot;
    thrust::device_vector<double> slope_y_mod = slope_y_mod_pot;

#else

    auto H = std::move(H_pot);
    auto u = std::move(u_pot);
    auto v = std::move(v_pot);
    auto additional_source_term = std::move(additional_source_term_pot);
    auto orography = std::move(orography_pot); 
    auto h_G = std::move(h_G_pot); 
    auto h_sd = std::move(h_sd_pot); 
    auto h_sn = std::move(h_sn_pot); 
    auto W_Gav = std::move(W_Gav_pot); 
    auto W_Gav_cum = std::move(W_Gav_cum_pot);
    auto hydraulic_conductivity = std::move(hydraulic_conductivity_pot);
    auto Z_Gav = std::move(Z_Gav_pot); 
    auto n_x = std::move(n_x_pot); 
    auto n_y = std::move(n_y_pot);
    auto u_star = std::move(u_star_pot); 
    auto v_star = std::move(v_star_pot); 
    auto h_interface_x = std::move(h_interface_x_pot);
    auto h_interface_y = std::move(h_interface_y_pot); 
    auto slope_x = std::move(slope_x_pot); 
    auto slope_y = std::move(slope_y_pot);
    auto soilMoistureRetention = std::move(soilMoistureRetention_pot); 
    auto eta = std::move(eta_pot); 
    
    auto idStaggeredInternalVectHorizontal_excluded = std::move(idStaggeredInternalVectHorizontal_excluded_pot);
    auto idStaggeredBoundaryVectWest_excluded = std::move(idStaggeredBoundaryVectWest_excluded_pot);
    auto idStaggeredBoundaryVectEast_excluded = std::move(idStaggeredBoundaryVectEast_excluded_pot);
    auto idStaggeredInternalVectVertical_excluded = std::move(idStaggeredInternalVectVertical_excluded_pot);
    auto idStaggeredBoundaryVectNorth_excluded = std::move(idStaggeredBoundaryVectNorth_excluded_pot);
    auto idStaggeredBoundaryVectSouth_excluded = std::move(idStaggeredBoundaryVectSouth_excluded_pot);
    auto idBasinVect_excluded = std::move(idBasinVect_excluded_pot);
    auto idBasinVectReIndex_excluded = std::move(idBasinVectReIndex_excluded_pot);

    auto idStaggeredInternalVectHorizontal = std::move(idStaggeredInternalVectHorizontal_pot);
    auto idStaggeredBoundaryVectWest = std::move(idStaggeredBoundaryVectWest_pot);
    auto idStaggeredBoundaryVectEast = std::move(idStaggeredBoundaryVectEast_pot);
    auto idStaggeredInternalVectVertical = std::move(idStaggeredInternalVectVertical_pot);
    auto idStaggeredBoundaryVectNorth = std::move(idStaggeredBoundaryVectNorth_pot);
    auto idStaggeredBoundaryVectSouth = std::move(idStaggeredBoundaryVectSouth_pot);
    auto idBasinVect = std::move(idBasinVect_pot);

    // set initial quantities for the time step adaptation
    auto H_old = H;
    auto H_oldold = H;

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

    Temperature temp(file_dir + temperature_file, N, max_Days, T_thr, 
                     std::round(max_Days * 24 / time_spacing_temp),
                     steps_per_hour, time_spacing_temp, height_thermometer,
                     format_temp, 
		     orography, idBasinVect); // here we pass the device vectors, because are not used during initialization 

    // +-----------------------------------------------+
    // |              Evapotranspiration               |
    // +-----------------------------------------------+

    evapoTranspiration ET(ET_model, N, temp.J, max_Days, phi_rad,
                          height_thermometer, 
			  idBasinVect, orography); // here we pass the device vectors, because are not used during initialization

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

    dt_DSV = maxdt(u, v, maxH, pixel_size, gravity, dt_DSV_given, t_final);

    c1_DSV_ = c1_DSV(dt_DSV, pixel_size);
    c2_DSV_ = c2_DSV(c1_DSV_);
    c3_DSV_ = c3_DSV(c1_DSV_);

    time = 0.; 
    timed = -dt_DSV;
    timedd = -2. * dt_DSV;

    // Up to this point, no computations are carried out on the device.
    // Just standard copies to Device are done (which are blocking)


    while (!is_last_step) {

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

      // Compute alfa coefficients
      alfa.f_x(CUDA_STREAM_ONLY(S(2)));
      alfa.f_y(CUDA_STREAM_ONLY(S(3)));

      // Update the temperature 
      temp.computeTemperature(time, dt_temp, dt_DSV CUDA_STREAM(S(4)));

      // ET varies daily
      if (std::floor(time / (24. * 3600)) >
          std::floor((time - dt_DSV) / (24. * 3600))) {
        ET.ET(temp.T_dailyMean, temp.T_dailyMin, temp.T_dailyMax,
              std::floor(time / (24 * 3600)) CUDA_STREAM(S(5)));
      }
 
      // +-----------------------------------------------+
      // |    Gravitational Layer, Snow Accumulation     |
      // +-----------------------------------------------+

      precipitation.computePrecipitation(time, soilMoistureRetention,
                                         temp.melt_mask, h_G, H, N_cols,
                                         idBasinVect CUDA_STREAM(S(6)));

      // update only if necessary
      if (std::floor(time / dt_min) > std::floor((time - dt_DSV) / dt_min)) {

        // vertical and horizontal residuals for Gravitational Layer
        computeResiduals(
            n_x, n_y, N_cols, hydraulic_conductivity,
            idStaggeredInternalVectHorizontal, idStaggeredInternalVectVertical,
            idStaggeredBoundaryVectWest, idStaggeredBoundaryVectEast,
            idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectSouth,
	    idBasinVect, temp.T_raster, temp.melt_mask,
	    h_sn, h_G, ET.ET_vec,
            precipitation.DP_infiltrated,
	    precipitation.DP_total,
	    c1_min,
	    dt_min, 
	    T_thr, 
            h_interface_x, h_interface_y CUDA_STREAM(S(7)));

      }

      // +-----------------------------------------------+
      // |               De Saint Venant                 |
      // +-----------------------------------------------+

      //  fill u_star and v_star with a Bilinear Interpolation
      bilinearInterpolation(u, v, u_star, v_star, N_rows, N_cols, dt_DSV,
                            pixel_size,
                            idStaggeredInternalVectHorizontal_excluded,
                            idStaggeredInternalVectVertical_excluded,
                            idStaggeredBoundaryVectWest_excluded,
                            idStaggeredBoundaryVectEast_excluded,
                            idStaggeredBoundaryVectNorth_excluded,
                            idStaggeredBoundaryVectSouth_excluded CUDA_STREAM(S(8)));

#ifndef ENABLE_CUDA
      coefficients.reserve(
          idBasinVect_excluded.size() +
          4 * idStaggeredInternalVectHorizontal_excluded.size() +
          4 * idStaggeredInternalVectVertical_excluded.size());

      additional_source_term.assign(N, 0.);
#endif

      buildMatrix(H_interface.horizontal, H_interface.vertical, orography,
                  u_star, v_star, u, v, H, N_cols, N_rows, N, c1_DSV_, c3_DSV_,
                  0, precipitation.DP_cumulative, dt_DSV, alfa.alfa_x,
                  alfa.alfa_y, idStaggeredInternalVectHorizontal_excluded,
                  idStaggeredInternalVectVertical_excluded,
                  idStaggeredBoundaryVectWest_excluded,
                  idStaggeredBoundaryVectEast_excluded,
                  idStaggeredBoundaryVectNorth_excluded,
                  idStaggeredBoundaryVectSouth_excluded, idBasinVect_excluded,
                  idBasinVect, idStaggeredInternalVectHorizontal,
                  idStaggeredInternalVectVertical, idBasinVectReIndex_excluded,
                  isNonReflectingBC, true,
		  additional_source_term CUDA_STREAM(S(9))
#ifndef ENABLE_CUDA
                  , excluded_ids, coefficients, rhs
#endif
		  );

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
/*
    //--------------------------------------------------------------------------
    // ### Run CG computation ###
    printf("CG loop:\n");
    gpu_CG(cublasHandle, cusparseHandle, m,
           matA, matL, d_B, d_X, d_R, d_R_aux, d_P, d_T,
           d_tmp, d_bufferMV, 
	   1000,  // max iterations 
	   1e-6); // tolerance
    //--------------------------------------------------------------------------
*/
#endif

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
                idStaggeredBoundaryVectSouth_excluded, isNonReflectingBC CUDA_STREAM(S(10)));
/*
      // +-----------------------------------------------+
      // |             Sediment Transport                |
      // +-----------------------------------------------+

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
            idStaggeredInternalVectHorizontal_excluded,
            idStaggeredInternalVectVertical_excluded,
            idStaggeredBoundaryVectWest_excluded,
            idStaggeredBoundaryVectEast_excluded,
            idStaggeredBoundaryVectNorth_excluded,
            idStaggeredBoundaryVectSouth_excluded, 
#ifdef ENABLE_CUDA
	    Gamma_vect_x_1, Gamma_vect_x_2, Gamma_vect_y_1, Gamma_vect_y_2,
	    slope_x_mod, slope_y_mod CUDA_STREAM(S(11))
#else
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
*/
    } // End Time Loop
      
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
