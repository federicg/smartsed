#ifndef CODE_INIT_H
#define CODE_INIT_H

//! std library
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

//! Eigen library
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

template <class T>
void resize_rasters(
    const T &dataFile,

    unsigned int &N_rows, unsigned int &N_cols, unsigned int &N,

    std::vector<unsigned int> &idStaggeredBoundaryVectSouth,
    std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    std::vector<unsigned int> &idStaggeredBoundaryVectEast,
    std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    std::vector<unsigned int> &idStaggeredInternalVectVertical,
    std::vector<unsigned int> &idBasinVect,
    std::vector<unsigned int> &idBasinVectReIndex,

    std::vector<unsigned int> &idStaggeredBoundaryVectSouth_excluded,
    std::vector<unsigned int> &idStaggeredBoundaryVectNorth_excluded,
    std::vector<unsigned int> &idStaggeredBoundaryVectWest_excluded,
    std::vector<unsigned int> &idStaggeredBoundaryVectEast_excluded,
    std::vector<unsigned int> &idStaggeredInternalVectHorizontal_excluded,
    std::vector<unsigned int> &idStaggeredInternalVectVertical_excluded,
    std::vector<unsigned int> &idBasinVect_excluded,
    std::vector<unsigned int> &idBasinVectReIndex_excluded,

    std::vector<std::array<double, 2>> &Gamma_vect_x,
    std::vector<std::array<double, 2>> &Gamma_vect_y,

    std::vector<std::tuple<bool, int>> &excluded_ids,
    std::vector<double> &additional_source_term,

    std::vector<double> &basin_mask_Vec, std::vector<double> &orography,
    std::vector<double> &h_G, std::vector<double> &h_sd,
    std::vector<double> &h_sn, std::vector<double> &S_coeff,
    std::vector<double> &W_Gav, std::vector<double> &W_Gav_cum,
    std::vector<double> &hydraulic_conductivity, std::vector<double> &Z_Gav,
    std::vector<double> &d_90, std::vector<double> &Res_x,
    std::vector<double> &Res_y, std::vector<double> &u, std::vector<double> &v,
    std::vector<double> &n_x, std::vector<double> &n_y,
    std::vector<double> &u_star, std::vector<double> &v_star,
    std::vector<double> &h_interface_x, std::vector<double> &h_interface_y,
    std::vector<double> &slope_x, std::vector<double> &slope_y,
    std::vector<double> &slope_cell, std::vector<double> &soilMoistureRetention,
    std::vector<double> &roughness_vect, std::vector<double> &eta,
    std::vector<double> &H,

    Eigen::VectorXd &H_basin, Eigen::VectorXd &rhs,

    double &pixel_size, // meter/pixel
    double &xllcorner, double &yllcorner, double &xllcorner_staggered_u,
    double &yllcorner_staggered_u, double &xllcorner_staggered_v,
    double &yllcorner_staggered_v, double &NODATA_value,

    const int currentSimNumber, const int rank, const int scaling_factor,
    const std::string friction_model, std::string &output_dir,
    const std::string file_dir, const bool save_temporal_sequence,
    const double max_Days, std::string &temperature_file, std::string &ET_model,
    std::string &infiltrationModel, double &phi_rad, double &dt_min,
    const double dt_temp, double &slope_x_max, double &slope_y_max,
    const int number_gauges, std::vector<std::vector<unsigned int>> &kk_gauges)

{

  output_dir = "../Outputs/" + std::to_string(currentSimNumber) + "/";

  const std::string orography_file =
      dataFile("files/orography_file", "DEM.txt");
  const std::string mask_file = dataFile("files/mask_file", "Mask_bin.txt");

  // precipitation files and temperature
  temperature_file =
      dataFile("files/meteo_data/temperature_file", "Temperature.txt");

  const bool restart_H = dataFile("files/initial_conditions/restart_H", false);
  const bool restart_vel =
      dataFile("files/initial_conditions/restart_vel", false);
  const bool restart_snow =
      dataFile("files/initial_conditions/restart_snow", false);
  const bool restart_sediment =
      dataFile("files/initial_conditions/restart_sediment", false);
  const bool restart_gravitational =
      dataFile("files/initial_conditions/restart_gravitational", false);
  const bool restart_soilMoisture =
      dataFile("files/initial_conditions/restart_soilMoisture", false);

  ET_model = dataFile("files/evapotranspiration/ET_model", "None");
  phi_rad =
      M_PI / 180. * dataFile("files/evapotranspiration/latitude_deg", 45.);

  infiltrationModel = dataFile("files/infiltration/infiltration_model", "None");

  const double delta_gauges = dataFile("discretization/delta_gauges", 0);

  const std::string corineCode_file =
      dataFile("files/infiltration/corineCode_file", "CLC_RASTER.txt");

  std::filesystem::create_directories(output_dir);

  // Resize the input Geometry
  {
    Raster orographyMat(file_dir + orography_file);
    Raster basin_mask(file_dir + mask_file);

    if (basin_mask.cellsize != orographyMat.cellsize) {
      if (rank == 0)
        std::cout << mask_file << " cellsize and " << orography_file
                  << " cellsize are not equal" << std::endl;
      exit(-1);
    }

    pixel_size = double(scaling_factor) * basin_mask.cellsize;

    if (rank == 0) {
      std::cout << "cell resolution for the current simulation = " << pixel_size
                << " meters" << std::endl;
      std::cout << "-------------------- " << std::endl;
    }

    xllcorner = basin_mask.xllcorner;
    yllcorner = basin_mask.yllcorner;
    NODATA_value = basin_mask.NODATA_value;
    if (basin_mask.cellsize <= pixel_size) {

      const std::string bashCommand =
          std::string("Rscript -e ") + "\"library(raster);" + "dem=raster('" +
          file_dir + orography_file + "');" + "basin=raster('" + file_dir +
          mask_file + "');" + "basin=aggregate(basin," +
          std::to_string(scaling_factor) + ");" +
          "values(basin)[values(basin)>0]=1;" +
          "dem=resample(dem,basin,method='bilinear');" +
          "writeRaster( dem, file=paste0('" + output_dir +
          "DEM.asc'), overwrite=TRUE );" + "writeRaster( basin, file=paste0('" +
          output_dir + "basin_mask.asc'), overwrite=TRUE )\"";
      std::system(bashCommand.c_str());
    } else {
      if (rank == 0)
        std::cout << "Basin mask greater than simulation resolution, i.e. "
                  << pixel_size << std::endl;
      exit(-1.);
    }

    if (dataFile("discretization/FillSinks", false)) {

      std::string bashCommand =
          std::string("python3 -c ") +
          "\"import os; import sys; import gdal; cwd = os.getcwd();" +
          "sys.path.append( cwd + '/../DeterministicProgram/include/py' ); "
          "import richdem as rd;" +
          "dem = rd.LoadGDAL( '" + output_dir + "DEM.asc' );" +
          "rd.FillDepressions( dem, in_place=True, "
          "epsilon=True,topology='D4' );" +
          "rd.SaveGDAL( '" + output_dir + "DEM.tif', dem )\"";

      std::system(bashCommand.c_str());

      bashCommand = std::string("Rscript -e ") + "\"library(raster);" +
                    "dem=raster('" + output_dir + "DEM.tif');" +
                    "writeRaster( dem, file=paste0('" + output_dir +
                    "DEM.asc'), overwrite=TRUE )\"";

      std::system(bashCommand.c_str());
    }
  }

  {

    Raster orographyMat(output_dir + "DEM.asc");
    Raster basin_mask(output_dir + "basin_mask.asc");

    N_rows = basin_mask.Coords.rows();
    N_cols = basin_mask.Coords.cols();
    N = N_rows * N_cols;

    H.resize(N);
    orography.resize(N);
    basin_mask_Vec.resize(N);
    eta.resize(N);
    h_G.resize(N);
    h_sd.resize(N);
    h_sn.resize(N);
    S_coeff.resize(N);
    W_Gav.resize(N);
    W_Gav_cum.resize(N);
    Res_x.resize(N);
    Res_y.resize(N);
    Z_Gav.resize(N);
    d_90.resize(N);
    soilMoistureRetention.resize(N);
    hydraulic_conductivity.resize(N);
    roughness_vect.resize(N);
    excluded_ids.resize(N);
    additional_source_term.resize(N);
    slope_cell.resize(N);

    u.resize((N_cols + 1) * N_rows);
    v.resize((N_rows + 1) * N_cols);
    n_x.resize(u.size());
    n_y.resize(v.size());
    u_star.resize(u.size());
    v_star.resize(v.size());

    Gamma_vect_x.resize(u.size());
    Gamma_vect_y.resize(v.size());

    h_interface_x.resize(u.size());
    h_interface_y.resize(v.size());

    slope_x.resize(u.size());
    slope_y.resize(v.size());

    for (unsigned int i = 0; i < N_rows; i++) {
      for (unsigned int j = 0; j < N_cols; j++) {
        const auto k = j + i * N_cols;
        basin_mask_Vec[k] = basin_mask.Coords.coeff(i, j) > 0;
        orography[k] = orographyMat.Coords.coeff(i, j);
      }
    }
  }

  if (restart_H) {
    {
      const std::string H_file =
          dataFile("files/initial_conditions/H_file", "H.txt");
      Raster HMat(file_dir + H_file);

      if (HMat.cellsize <= pixel_size) {
        const std::string bashCommand =
            std::string("Rscript -e ") + "\"library(raster);" +
            "basin=raster('" + output_dir + "basin_mask.asc" + "');" +
            "H=raster('" + file_dir + H_file + "');" +
            "H=resample(H,basin,method='bilinear');" +
            "values(H)[is.na(values(H))]=0;" + "writeRaster( H, file=paste0('" +
            output_dir + "H_0.asc'), overwrite=TRUE )\"";
        std::system(bashCommand.c_str());

      } else {
        if (rank == 0)
          std::cout << "Error! resolution of surface water is greater than "
                       "simulation resolution, i.e. "
                    << pixel_size << std::endl;
        exit(-1.);
      }
    }
    Raster HMat(output_dir + "H_0.asc");

    for (unsigned int i = 0; i < N_rows; i++) {
      for (unsigned int j = 0; j < N_cols; j++) {
        const auto k = j + i * N_cols;
        H[k] = HMat.Coords.coeff(i, j) * basin_mask_Vec[k];
        eta[k] = H[k] + orography[k];
      }
    }

  } else {
    for (unsigned int i = 0; i < N; i++) {
      H[i] = 0.;
      eta[i] = 0.;
    }

    saveSolution(output_dir + "H_0", " ", N_rows, N_cols, xllcorner, yllcorner,
                 pixel_size, NODATA_value, H);
  }

  if (restart_vel) {
    {
      const std::string restart_vel_u_file =
          dataFile("files/initial_conditions/vel_file_u", "restart_vel_u.txt");
      Raster vel_u_Mat(file_dir + restart_vel_u_file);

      if (vel_u_Mat.cellsize <= pixel_size) {
        const std::string bashCommand =
            std::string("Rscript -e ") + "\"library(raster);" + "dem=raster('" +
            output_dir + "DEM.asc" + "');" + "u=raster('" + file_dir +
            restart_vel_u_file + "');" + "dem_u=dem;" +
            "extent(dem_u)@xmin=extent(dem)@xmin-res(dem)[[1]]/2;" +
            "extent(dem_u)@xmax=extent(dem)@xmax+res(dem)[[1]]/2;" +
            "res(dem_u)=res(dem);" +
            "extent(u)@xmin=extent(dem)@xmin-res(u)[[1]]/2;" +
            "extent(u)@xmax=extent(dem)@xmax+res(u)[[1]]/2;" +
            "u=resample(u,dem_u,method='bilinear');" +
            "values(u)[is.na(values(u))]=0;" + "writeRaster( u, file=paste0('" +
            output_dir + "u_0.asc'), overwrite=TRUE )\"";
        std::system(bashCommand.c_str());

      } else {
        if (rank == 0)
          std::cout << "Error! resolution of horizontal velocity is greater "
                       "than simulation resolution i.e. "
                    << pixel_size << std::endl;
        exit(-1.);
      }
    }
    Raster vel_u_Mat(output_dir + "u_0.asc");
    for (unsigned int i = 0; i < N_rows; i++) {
      for (unsigned int j = 0; j <= N_cols; j++) {
        const auto Id = j + (N_cols + 1) * i;
        u[Id] = vel_u_Mat.Coords.coeff(i, j);
      }
    }

    xllcorner_staggered_u = vel_u_Mat.xllcorner;
    yllcorner_staggered_u = vel_u_Mat.yllcorner;

    {
      const std::string restart_vel_v_file =
          dataFile("files/initial_conditions/vel_file_v", "restart_vel_v.txt");
      Raster vel_v_Mat(file_dir + restart_vel_v_file);

      if (vel_v_Mat.cellsize <= pixel_size) {
        const std::string bashCommand =
            std::string("Rscript -e ") + "\"library(raster);" + "dem=raster('" +
            output_dir + "DEM.asc" + "');" + "v=raster('" + file_dir +
            restart_vel_v_file + "');" + "dem_v=dem;" +
            "extent(dem_v)@ymin=extent(dem)@ymin-res(dem)[[1]]/2;" +
            "extent(dem_v)@ymax=extent(dem)@ymax+res(dem)[[1]]/2;" +
            "res(dem_v)=res(dem);" +
            "extent(v)@ymin=extent(dem)@ymin-res(v)[[1]]/2;" +
            "extent(v)@ymax=extent(dem)@ymax+res(v)[[1]]/2;" +
            "v=resample(v,dem_v,method='bilinear');" +
            "values(v)[is.na(values(v))]=0;" + "writeRaster( v, file=paste0('" +
            output_dir + "v_0.asc'), overwrite=TRUE )\"";
        std::system(bashCommand.c_str());

      } else {
        if (rank == 0)
          std::cout << "Error! resolution of vertical velocity is greater "
                       "than simulation resolution, i.e. "
                    << pixel_size << std::endl;
        exit(-1.);
      }
    }
    Raster vel_v_Mat(output_dir + "v_0.asc");

    xllcorner_staggered_v = vel_v_Mat.xllcorner;
    yllcorner_staggered_v = vel_v_Mat.yllcorner;

    for (unsigned int i = 0; i <= N_rows; i++) {
      for (unsigned int j = 0; j < N_cols; j++) {
        const auto Id = j + i * N_cols;
        v[Id] = vel_v_Mat.Coords.coeff(i, j);
      }
    }

  } else {
    xllcorner_staggered_u = xllcorner - pixel_size / 2.;
    yllcorner_staggered_u = yllcorner;

    u.assign(u.size(), 0.0);

    saveSolution(output_dir + "u_0", "u", N_rows, N_cols, xllcorner, yllcorner,
                 pixel_size, NODATA_value, u);

    xllcorner_staggered_v = xllcorner;
    yllcorner_staggered_v = yllcorner - pixel_size / 2.;

    v.assign(v.size(), 0.0);

    saveSolution(output_dir + "v_0", "v", N_rows, N_cols, xllcorner, yllcorner,
                 pixel_size, NODATA_value, v);
  }

  if (restart_snow) {
    {
      const std::string restart_snow_file =
          dataFile("files/initial_conditions/snow_file", "h_sn.asc");
      Raster snow_Mat(file_dir + restart_snow_file);

      if (snow_Mat.cellsize <= pixel_size) {
        const std::string bashCommand =
            std::string("Rscript -e ") + "\"library(raster);" + "dem=raster('" +
            output_dir + "DEM.asc" + "');" + "hsn=raster('" + file_dir +
            restart_snow_file + "');" +
            "hsn=resample(hsn,dem,method='bilinear');" +
            "values(hsn)[is.na(values(hsn))]=0;" +
            "writeRaster( hsn, file=paste0('" + output_dir +
            "hsn_0.asc'), overwrite=TRUE )\"";
        std::system(bashCommand.c_str());

      } else {
        if (rank == 0)
          std::cout << "Error! resolution of snow file is greater than "
                       "simulation resolution, i.e. "
                    << pixel_size << std::endl;
        exit(-1.);
      }
    }
    Raster snow_Mat(output_dir + "hsn_0.asc");

    for (unsigned int i = 0; i < N_rows; i++) {
      for (unsigned int j = 0; j < N_cols; j++) {
        const auto k = j + i * N_cols;
        h_sn[k] = snow_Mat.Coords.coeff(i, j) * basin_mask_Vec[k];
      }
    }

  } else {
    h_sn.assign(h_sn.size(), 0.0);
    saveSolution(output_dir + "hsn_0", " ", N_rows, N_cols, xllcorner,
                 yllcorner, pixel_size, NODATA_value, h_sn);
  }

  if (restart_sediment) {
    {
      const std::string restart_sediment_file =
          dataFile("files/initial_conditions/sediment_file", "h_sd.asc");
      Raster sediment_Mat(file_dir + restart_sediment_file);

      if (sediment_Mat.cellsize <= pixel_size) {
        const std::string bashCommand =
            std::string("Rscript -e ") + "\"library(raster);" + "dem=raster('" +
            output_dir + "DEM.asc" + "');" + "hsd=raster('" + file_dir +
            restart_sediment_file + "');" +
            "hsd=resample(hsd,dem,method='bilinear');" +
            "values(hsd)[is.na(values(hsd))]=0;" +
            "writeRaster( hsd, file=paste0('" + output_dir +
            "hsd_0.asc'), overwrite=TRUE )\"";
        std::system(bashCommand.c_str());

      } else {
        if (rank == 0)
          std::cout << "Error! resolution of sediment file is greater than "
                       "simulation resolution, i.e. "
                    << pixel_size << std::endl;
        exit(-1.);
      }
    }
    Raster sediment_Mat(output_dir + "hsd_0.asc");

    for (unsigned int i = 0; i < N_rows; i++) {
      for (unsigned int j = 0; j < N_cols; j++) {
        const auto k = j + i * N_cols;
        h_sd[k] = sediment_Mat.Coords.coeff(i, j) * basin_mask_Vec[k];
      }
    }
  } else {
    h_sd.assign(h_sd.size(), 0.0);
    saveSolution(output_dir + "hsd_0", " ", N_rows, N_cols, xllcorner,
                 yllcorner, pixel_size, NODATA_value, h_sd);
  }

  if (restart_gravitational) {
    {
      const std::string restart_gravitational_file =
          dataFile("files/initial_conditions/gravitational_file", "h_G.asc");
      Raster gravitational_Mat(file_dir + restart_gravitational_file);

      if (gravitational_Mat.cellsize <= pixel_size) {
        const std::string bashCommand =
            std::string("Rscript -e ") + "\"library(raster);" + "dem=raster('" +
            output_dir + "DEM.asc" + "');" + "hG=raster('" + file_dir +
            restart_gravitational_file + "');" +
            "hG=resample(hG,dem,method='bilinear');" +
            "values(hG)[is.na(values(hG))]=0;" +
            "writeRaster( hG, file=paste0('" + output_dir +
            "hG_0.asc'), overwrite=TRUE )\"";
        std::system(bashCommand.c_str());

      } else {
        if (rank == 0)
          std::cout << "Error! resolution of gravitational file is greater "
                       "than simulation resolution, i.e. "
                    << pixel_size << std::endl;
        exit(-1.);
      }
    }
    Raster gravitational_Mat(output_dir + "hG_0.asc");

    for (unsigned int i = 0; i < N_rows; i++) {
      for (unsigned int j = 0; j < N_cols; j++) {
        const auto k = j + i * N_cols;
        h_G[k] = gravitational_Mat.Coords.coeff(i, j) * basin_mask_Vec[k];
      }
    }

  } else {
    h_G.assign(h_G.size(), 0.0);
    saveSolution(output_dir + "hG_0", " ", N_rows, N_cols, xllcorner, yllcorner,
                 pixel_size, NODATA_value, h_G);
  }

  // +-----------------------------------------------+
  // |    Construct soilMoistureRetention vector     |
  // +-----------------------------------------------+

  std::vector<int> corineCode_Vec(N);
  std::vector<double> X_Gav(N), Y_Gav(N);

  {
    const std::string check_presence_string = file_dir + corineCode_file;

    if (!is_file_exist(check_presence_string.c_str())) {
      if (rank == 0)
        std::cout << check_presence_string << " is not present!" << std::endl;
      exit(-1.);
    }

    // interpolate CLC to make sure to match correct dimensions
    {
      const std::string bashCommand =
          std::string("Rscript -e ") + "\"library(raster);" + "dem=raster('" +
          output_dir + "DEM.asc" + "');" + "clc=raster('" + file_dir +
          corineCode_file + "');" + "clc=resample(clc,dem,method='ngb');" +
          "values(clc)[is.na(values(clc))]=0;" +
          "writeRaster( clc, file=paste0('" + output_dir +
          "CLC.asc'), overwrite=TRUE )\"";
      std::system(bashCommand.c_str());
    }

    Raster corineCode(output_dir + "CLC.asc");

    if (corineCode.cellsize != pixel_size) {
      if (rank == 0)
        std::cout << "Please check that the " << corineCode_file
                  << " cellsize is consistent with " << mask_file << " and "
                  << orography_file << " ones" << std::endl;
      exit(-1.);
    }

    for (unsigned int i = 0; i < N_rows; i++) {
      for (unsigned int j = 0; j < N_cols; j++) {
        const auto k = j + i * N_cols;
        corineCode_Vec[k] = corineCode.Coords.coeff(i, j);
      }
    }
  }

  {
    std::vector<int> HSG(N);
    std::vector<double> clayPercentage_Vec(N), sandPercentage_Vec(N);

    std::string str1, str2;
    if (restart_soilMoisture) {
      const std::string restart_clay_file = dataFile(
                            "files/initial_conditions/clay_file", "clay.asc"),
                        restart_sand_file = dataFile(
                            "files/initial_conditions/sand_file", "sand.asc");

      str1 = file_dir + restart_clay_file;
      str2 = file_dir + restart_sand_file;
    } else {
      str1 =
          output_dir + "clay_sim_" + std::to_string(currentSimNumber) + ".asc";
      str2 =
          output_dir + "sand_sim_" + std::to_string(currentSimNumber) + ".asc";

      if (!is_file_exist(str1.c_str())) {
        if (rank == 0)
          std::cout << str1
                    << " is not present! Make sure you have put nsim>=0 or "
                       "in case you want to provide directly the particle "
                       "size fractions "
                       "make sure you have put restart_soilMoisture=true and "
                       "specified "
                       "the correct paths"
                    << std::endl;
        exit(-1.);
      }

      if (!is_file_exist(str2.c_str())) {
        if (rank == 0)
          std::cout << str2
                    << " is not present! Make sure you have put nsim>=0 or "
                       "in case you want to provide directly the particle "
                       "size fractions "
                       "make sure you have put restart_soilMoisture=true and "
                       "specified "
                       "the correct paths"
                    << std::endl;
        exit(-1.);
      }
    }

    double cellsize_psf = 0;
    if (infiltrationModel != "None" || friction_model == "Rickenmann") {

      // interpolate psfs to make sure to match correct dimensions
      if (restart_soilMoisture) {
        std::string bashCommand;

        bashCommand = std::string("Rscript -e ") + "\"library(raster);" +
                      "dem=raster('" + output_dir + "DEM.asc" + "');" +
                      "clay=raster('" + str1 + "');" +
                      "clay=resample(clay,dem,method='ngb');" +
                      "values(clay)[is.na(values(clay))]=0;" +
                      "writeRaster( clay, file=paste0('" + output_dir +
                      "clay.asc'), overwrite=TRUE )\"";
        std::system(bashCommand.c_str());

        bashCommand = std::string("Rscript -e ") + "\"library(raster);" +
                      "dem=raster('" + output_dir + "DEM.asc" + "');" +
                      "sand=raster('" + str2 + "');" +
                      "sand=resample(sand,dem,method='ngb');" +
                      "values(sand)[is.na(values(sand))]=0;" +
                      "writeRaster( sand, file=paste0('" + output_dir +
                      "sand.asc'), overwrite=TRUE )\"";
        std::system(bashCommand.c_str());
      }

      Raster clayPercentage(str1), sandPercentage(str2);

      cellsize_psf = clayPercentage.cellsize;
      if (cellsize_psf != sandPercentage.cellsize) {
        if (rank == 0)
          std::cout << "Please check that the soil texture files have the "
                       "same resolution"
                    << std::endl;
        exit(-1.);
      }

      for (unsigned int i = 0; i < N_rows; i++) {
        for (unsigned int j = 0; j < N_cols; j++) {
          const auto k = j + i * N_cols;

          clayPercentage_Vec[k] = clayPercentage.Coords.coeff(i, j);
          sandPercentage_Vec[k] = sandPercentage.Coords.coeff(i, j);
        }
      }
    }

    if (infiltrationModel != "None" || friction_model == "Rickenmann") {

      for (unsigned int i = 0; i < N_rows; i++) {
        for (unsigned int j = 0; j < N_cols; j++) {

          const auto k = j + i * N_cols;

          const auto &clay = clayPercentage_Vec[k],
                     &sand = sandPercentage_Vec[k];

          if (sand > .9 && sand <= 1 && clay < .1 && clay >= 0) // A
          {
            HSG[k] = 0;
          } else if (sand > .5 && sand < .9 && clay > .1 && clay < .2) // B
          {
            HSG[k] = 1;
          } else if (sand < .5 && sand >= 0 && clay > .2 && clay < .4) // C
          {
            HSG[k] = 2;
          } else if (sand < .5 && sand >= 0 && clay > .4 && clay <= 1) // D
          {
            HSG[k] = 3;
          } else if (sand >= 0 && sand <= 1 && clay >= 0 && clay <= 1) {

            Vector2D point(std::array<double, 2>{{clay, sand}});

            Vector2D point_A(std::array<double, 2>{{0, 1}});
            Vector2D point_B(std::array<double, 2>{{.1, 1}});
            Vector2D point_C(std::array<double, 2>{{.1, .9}});
            Vector2D point_D(std::array<double, 2>{{0, .9}});

            Vector2D point_E(std::array<double, 2>{{.1, .5}});
            Vector2D point_F(std::array<double, 2>{{.2, .5}});
            Vector2D point_G(std::array<double, 2>{{.2, .9}});

            Vector2D point_H(std::array<double, 2>{{.2, 0}});
            Vector2D point_I(std::array<double, 2>{{.4, 0}});
            Vector2D point_L(std::array<double, 2>{{.4, .5}});

            Vector2D point_M(std::array<double, 2>{{1, 0}});
            Vector2D point_N(std::array<double, 2>{{1, .5}});

            std::vector<Vector2D> vv = {point_A, point_D, point_C, point_B,

                                        point_C, point_E, point_F, point_G,

                                        point_F, point_H, point_I, point_L,

                                        point_L, point_I, point_M, point_N};

            std::pair<double, int> min = std::make_pair(1.e4, -1);

            for (auto ii = 0; ii < vv.size(); ii += 4) {

              const auto A = vv[ii], D = vv[ii + 1], C = vv[ii + 2],
                         B = vv[ii + 3];

              double d1 = 1.e4, d2 = 1.e4, d3 = 1.e4, d4 = 1.e4;

              Vector2D e1(std::array<double, 2>{{1, 0}}),
                  e2(std::array<double, 2>{{0, 1}});

              if (e1.dot(point - D) >= 0 && e1.dot(point - D) <= (C(0) - D(0)))
                d1 = std::abs((point - D).dot(e2));
              if (e2.dot(point - D) >= 0 && e2.dot(point - D) <= (A(1) - D(1)))
                d2 = std::abs((point - D).dot(e1));
              if (e1.dot(point - A) >= 0 && e1.dot(point - A) <= (B(0) - A(0)))
                d3 = std::abs((point - A).dot(e2));
              if (e2.dot(point - C) >= 0 && e2.dot(point - C) <= (B(1) - C(1)))
                d4 = std::abs((point - C).dot(e1));

              if (d1 < min.first)
                min = std::pair<double, int>(d1, ii / 4.);
              if (d2 < min.first)
                min = std::pair<double, int>(d2, ii / 4.);
              if (d3 < min.first)
                min = std::pair<double, int>(d3, ii / 4.);
              if (d4 < min.first)
                min = std::pair<double, int>(d4, ii / 4.);
            }

            const auto &id = min.second;
            if (id == 0) // A
            {
              HSG[k] = 0;
            } else if (id == 1) // B
            {
              HSG[k] = 1;
            } else if (id == 2) // C
            {
              HSG[k] = 2;
            } else if (id == 3) // D
            {
              HSG[k] = 3;
            } else {
              if (rank == 0)
                std::cout << "Something wrong in HSG classification"
                          << std::endl;
              exit(-1.);
            }

          } else {
            HSG[k] = -1;
          }
        }
      }
    }

    static constexpr double S_0 = .254; // 254 mm

    const auto CN_map = createCN_map();

    for (unsigned int k = 0; k < N; k++) {
      const auto key = std::array<int, 2>{{corineCode_Vec[k], HSG[k]}};

      const auto it = CN_map.find(key);

      if (it != CN_map.end() && infiltrationModel != "None") {
        soilMoistureRetention[k] =
            S_0 * (100. / double(it->second) - 1.) * basin_mask_Vec[k];
      } else {
        soilMoistureRetention[k] = 0.;
      }
    }

    // build X_Gav and Y_Gav
    const std::string Gavrilovic_file_input =
        dataFile("physics/Gavrilovic_txt", "Gav.txt");
    const auto CN_Gav_map = createCN_map_Gav(file_dir + Gavrilovic_file_input);

    for (unsigned int k = 0; k < N; k++) {
      const auto key = corineCode_Vec[k];

      const auto it = CN_Gav_map.find(key);

      if (it != CN_Gav_map.end()) // remains zero else
      {
        X_Gav[k] = it->second[0];
        Y_Gav[k] = it->second[1];
      }

      Z_Gav[k] = X_Gav[k] * Y_Gav[k];
    }

    if (clayPercentage_Vec[0] == 0 && sandPercentage_Vec[0] == 0 &&
        friction_model == "Rickenmann") {
      if (rank == 0)
        std::cout
            << "clay and sand are both zero (can't compute d90 for friction, "
               "maybe change friction_model in SMARTSED_input in Manning if "
               "you don't want to run the R script), probably you have not "
               "run correctly the Geostatistical preprocessor!, STOP!"
            << std::endl;
      exit(1.);
    }

    // build d_10 (for k_c) and d_90 (frictionClass)
    auto d_10 = d_90;
    if (infiltrationModel != "None") {
      d_10 = compute_d_perc(clayPercentage_Vec, sandPercentage_Vec, 10);

      // Equations for hydraulic conductivity estimation from particle size
      // distribution: A dimensional analysis Ji-Peng Wang1, Bertrand
      // François, and Pierre Lambert
      static constexpr double C_H = 6.54e-4;
      static constexpr double gravity = 9.81;
      static constexpr double kin_visc = 0.89e-6;
      for (unsigned int i = 0; i < N; i++) {
        hydraulic_conductivity[i] =
            C_H * gravity / kin_visc * std::pow(d_10[i], 2.);
      }
    }

    if (infiltrationModel != "None" || friction_model == "Rickenmann")
      d_90 = compute_d_perc(clayPercentage_Vec, sandPercentage_Vec, 90);

    saveSolution(output_dir + "soilMoistureRetention", " ", N_rows, N_cols,
                 xllcorner, yllcorner, pixel_size, NODATA_value,
                 soilMoistureRetention);
    saveSolution(output_dir + "d_10", " ", N_rows, N_cols, xllcorner, yllcorner,
                 pixel_size, NODATA_value, d_10);
    saveSolution(output_dir + "d_90", " ", N_rows, N_cols, xllcorner, yllcorner,
                 pixel_size, NODATA_value, d_90);
    saveSolution(output_dir + "k_c", " ", N_rows, N_cols, xllcorner, yllcorner,
                 pixel_size, NODATA_value, hydraulic_conductivity);
  }

  if (rank == 0)
    std::cout << "maximum and minimum hydraulic_conductivity  "
              << *std::max_element(hydraulic_conductivity.begin(),
                                   hydraulic_conductivity.end())
              << " "
              << *std::min_element(hydraulic_conductivity.begin(),
                                   hydraulic_conductivity.end())
              << std::endl;

  // +-----------------------------------------------+
  // |                Compute Slopes                 |
  // +-----------------------------------------------+

  for (unsigned int i = 0; i < N_rows; i++) {
    for (unsigned int j = 1; j < N_cols; j++) {
      const auto Id = j + i * (N_cols + 1);

      slope_x[Id] = (orography[Id - i] - orography[Id - 1 - i]) / pixel_size;
      n_x[Id] = -(orography[Id - i] - orography[Id - 1 - i]) /
                std::abs((orography[Id - i] - orography[Id - 1 - i]));

      if (std::isnan(n_x[Id]))
        n_x[Id] = 0;
    }
  }

  for (unsigned int i = 0, j = 0; i < N_rows; i++) {
    const auto Id = j + i * (N_cols + 1);

    slope_x[Id] = slope_x[Id + 1];
    n_x[Id] = n_x[Id + 1];
  }

  for (unsigned int i = 0, j = N_cols; i < N_rows; i++) {
    const auto Id = j + i * (N_cols + 1);

    slope_x[Id] = slope_x[Id - 1];
    n_x[Id] = n_x[Id - 1];
  }

  for (unsigned int i = 1; i < N_rows; i++) {
    for (unsigned int j = 0; j < N_cols; j++) {
      const auto Id = j + i * N_cols;

      slope_y[Id] = (orography[Id] - orography[Id - N_cols]) / pixel_size;
      n_y[Id] = -(orography[Id] - orography[Id - N_cols]) /
                std::abs((orography[Id] - orography[Id - N_cols]));

      if (std::isnan(n_y[Id]))
        n_y[Id] = 0;
    }
  }

  for (unsigned int j = 0, i = 0; j < N_cols; j++) {
    const auto Id = j + i * N_cols;

    slope_y[Id] = slope_y[Id + N_cols];
    n_y[Id] = n_y[Id + N_cols];
  }

  for (unsigned int j = 0, i = N_rows; j < N_cols; j++) {
    const auto Id = j + i * N_cols;

    slope_y[Id] = slope_y[Id - N_cols];
    n_y[Id] = n_y[Id - N_cols];
  }

  saveSolution(output_dir + "slope_x", "u", N_rows, N_cols, xllcorner,
               yllcorner, pixel_size, NODATA_value, slope_x);
  saveSolution(output_dir + "slope_y", "v", N_rows, N_cols, xllcorner,
               yllcorner, pixel_size, NODATA_value, slope_y);

  // +-----------------------------------------------+
  // |     Compute boundaries of basin domain        |
  // +-----------------------------------------------+

  computeAdjacencies(basin_mask_Vec, idStaggeredBoundaryVectSouth,
                     idStaggeredBoundaryVectNorth, idStaggeredBoundaryVectWest,
                     idStaggeredBoundaryVectEast,
                     idStaggeredInternalVectHorizontal,
                     idStaggeredInternalVectVertical, idBasinVect,
                     idBasinVectReIndex, N_rows, N_cols);

  // +-----------------------------------------------+
  // |                Gavrilovic Coeff.              |
  // +-----------------------------------------------+

  for (const auto &k : idBasinVect) {
    const unsigned int i = k / N_cols;
    slope_cell[k] =
        std::sqrt(std::pow(.5 * (slope_x[k + i] + slope_x[k + i + 1]), 2.) +
                  std::pow(.5 * (slope_y[k] + slope_y[k + N_cols]), 2.));

    Z_Gav[k] *= (.5 + std::sqrt(slope_cell[k]));
    Z_Gav[k] = std::pow(Z_Gav[k], 1.5);
  }

  // +-----------------------------------------------+
  // |                    Gauges i,j                 |
  // +-----------------------------------------------+

  kk_gauges;
  kk_gauges.resize(number_gauges);

  for (int number = 1; number <= number_gauges; number++) {
    std::string filename_x = "discretization/X_gauges_",
                filename_y = "discretization/Y_gauges_";

    filename_x += std::to_string(number);
    filename_y += std::to_string(number);

    const double X_gauges = dataFile(filename_x.c_str(), 0.);
    const double Y_gauges = dataFile(filename_y.c_str(), 0.);

    const Vector2D XX_gauges(std::array<double, 2>{{X_gauges, Y_gauges}});

    if (save_temporal_sequence) {

      const Vector2D XX_O =
          std::array<double, 2>{{xllcorner, yllcorner + N_rows * pixel_size}};

      auto XX = (XX_gauges - XX_O) / pixel_size; // coordinate in the matrix

      auto XX_east =
          XX + Vector2D(std::array<double, 2>{{delta_gauges / pixel_size, 0}});
      auto XX_west =
          XX - Vector2D(std::array<double, 2>{{delta_gauges / pixel_size, 0}});

      auto XX_south =
          XX + Vector2D(std::array<double, 2>{{0, delta_gauges / pixel_size}});
      auto XX_north =
          XX - Vector2D(std::array<double, 2>{{0, delta_gauges / pixel_size}});

      XX(1) = -std::round(XX(1));
      XX(0) = std::round(XX(0));
      if (XX(0) < 0 || XX(1) < 0 || XX(1) >= N_rows || XX(0) >= N_cols) {
        if (rank == 0)
          std::cout << "The gauges in the input file are not good" << std::endl;
        exit(1.);
      }

      int i_1 = -std::round(XX_north(1));
      int i_2 = -std::round(XX_south(1));

      int j_1 = std::round(XX_west(0));
      int j_2 = std::round(XX_east(0));

      i_1 = std::min(std::max(i_1, 0), int(N_rows - 1));
      j_1 = std::min(std::max(j_1, 0), int(N_cols - 1));

      i_2 = std::min(std::max(i_2, 0), int(N_rows - 1));
      j_2 = std::min(std::max(j_2, 0), int(N_cols - 1));

      for (int i = i_2; i <= i_1; i++) {
        for (int j = j_1; j <= j_2; j++) {
          kk_gauges[number - 1].push_back(i * N_cols + j);
        }
      }
    }
  }

  // +-----------------------------------------------+
  // |           Beginning of core Part              |
  // +-----------------------------------------------+

  for (int ii = 0; ii < N; ii++) {
    const auto r1 =
        dataFile("files/infiltration/roughness_scale_factor1", 100.);
    const auto r2 =
        dataFile("files/infiltration/roughness_scale_factor2", 100.);
    const auto r3 =
        dataFile("files/infiltration/roughness_scale_factor3", 100.);

    if (slope_cell[ii] <= 0.2) {
      roughness_vect[ii] = r1;
    } else if (slope_cell[ii] <= 0.6) {
      roughness_vect[ii] = r2;
    } else {
      roughness_vect[ii] = r3;
    }
  }

  slope_y_max = 0.;
  slope_x_max = 0.;
  for (const auto &k : idBasinVect) {
    const unsigned int i = k / N_cols;

    const auto slope_x_l = std::abs(slope_x[k + i]);
    const auto slope_x_r = std::abs(slope_x[k + i + 1]);
    const auto slope_y_l = std::abs(slope_y[k]);
    const auto slope_y_r = std::abs(slope_y[k + N_cols]);

    if (slope_x_l > slope_x_max)
      slope_x_max = slope_x_l;
    if (slope_x_r > slope_x_max)
      slope_x_max = slope_x_r;
    if (slope_y_l > slope_y_max)
      slope_y_max = slope_y_l;
    if (slope_y_r > slope_y_max)
      slope_y_max = slope_y_r;
  }

  // +-----------------------------------------------+
  // |              compute_sub_basins               |
  // +-----------------------------------------------+

  excluded_ids.assign(N, std::tuple<bool, int>(false, -1));
  additional_source_term.assign(N, 0.);

  const bool static_subbasin_approx =
      dataFile("discretization/static_subbasin_approx", false);
  if (static_subbasin_approx) {

    const std::set<unsigned int> idBasinVect_set(idBasinVect.begin(),
                                                 idBasinVect.end()),
        idStaggeredBoundaryVectSouth_set(idStaggeredBoundaryVectSouth.begin(),
                                         idStaggeredBoundaryVectSouth.end()),
        idStaggeredBoundaryVectNorth_set(idStaggeredBoundaryVectNorth.begin(),
                                         idStaggeredBoundaryVectNorth.end()),
        idStaggeredBoundaryVectWest_set(idStaggeredBoundaryVectWest.begin(),
                                        idStaggeredBoundaryVectWest.end()),
        idStaggeredBoundaryVectEast_set(idStaggeredBoundaryVectEast.begin(),
                                        idStaggeredBoundaryVectEast.end());

    const double slope_thr = dataFile("discretization/slope_thr", 1.);

    // exclude high slopes,
    for (const unsigned int &Id : idBasinVect) {
      const auto &current_slope_cell = slope_cell[Id];
      if (current_slope_cell > slope_thr) {
        std::get<0>(excluded_ids[Id]) = true;
      }
    }

    // exclude also isolated cells, see below,
    // maybe cycle on intefaces, build a list of bool for each id
    std::vector<unsigned int> counter_near_excl;
    counter_near_excl.resize(N);
    counter_near_excl.assign(N, 0);

    // cycle over interfaces, internal vertical and horizontal
    for (const auto &Id : idStaggeredInternalVectHorizontal) {
      const unsigned int i = Id / (N_cols + 1),

                         IDleft = Id - i - 1, // H
          IDright = Id - i;

      if (std::get<0>(excluded_ids[IDleft])) {
        counter_near_excl[IDright] += 1;
      }
      if (std::get<0>(excluded_ids[IDright])) {
        counter_near_excl[IDleft] += 1;
      }
    }

    for (const auto &Id : idStaggeredInternalVectVertical) {

      const unsigned int IDleft = Id - N_cols, // H
          IDright = Id;

      if (std::get<0>(excluded_ids[IDleft])) {
        counter_near_excl[IDright] += 1;
      }
      if (std::get<0>(excluded_ids[IDright])) {
        counter_near_excl[IDleft] += 1;
      }
    }

    for (const auto &Id : idBasinVect) {
      if (counter_near_excl[Id] == 4) // it is surely an isolated cell! (can be
                                      // also an excluded cell but no problem..)
      {
        std::get<0>(excluded_ids[Id]) = true;
      }
    }

    // compute pour points, local movement in the 8 directions! (also diagonal
    // ones)
    for (const unsigned int &Id : idBasinVect) {

      auto &current_tuple = excluded_ids[Id];

      const auto &current_is = std::get<0>(current_tuple);
      auto &current_pour_id = std::get<1>(current_tuple);

      if (current_is) {
        // start from Id and perform a local search (by means of the gradient)
        // to get the pour point
        int candidate_pour_id = Id;

        while (true) {

          candidate_pour_id = computePourCell(
              candidate_pour_id, N_cols, orography, idBasinVect_set,
              idStaggeredBoundaryVectSouth_set,
              idStaggeredBoundaryVectNorth_set, idStaggeredBoundaryVectWest_set,
              idStaggeredBoundaryVectEast_set);

          // if the gradient points outside the basin leave current_pour_id =
          // -1
          if (candidate_pour_id < 0)
            break;

          const auto &is_current_id_again_excluded =
              std::get<0>(excluded_ids[candidate_pour_id]);
          if (!is_current_id_again_excluded) {
            current_pour_id = candidate_pour_id;
            break;
          }
        }
      }
    }
  }

  computeAdjacencies(
      basin_mask_Vec, excluded_ids, idStaggeredBoundaryVectSouth_excluded,
      idStaggeredBoundaryVectNorth_excluded,
      idStaggeredBoundaryVectWest_excluded,
      idStaggeredBoundaryVectEast_excluded,
      idStaggeredInternalVectHorizontal_excluded,
      idStaggeredInternalVectVertical_excluded, idBasinVect_excluded,
      idBasinVectReIndex_excluded, N_rows, N_cols);

  saveSolution(output_dir + "excluded_ids", N_rows, N_cols, xllcorner,
               yllcorner, pixel_size, NODATA_value, excluded_ids);
}

#endif
