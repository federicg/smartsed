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

    @author      Federico     Gatti        MOX Politecnico di Milano               <federico.gatti@polimi.it>
    @mantainer

    @supervisors Luca         Bonaventura  MOX Politecnico di Milano               <luca.bonaventura@polimi.it>
                 Alessandra   Menafoglio   MOX Politecnico di Milano               <alessandra.menafoglio@polimi.it>
                 Laura        Longoni      Applied Geology Politecnico di Milano   <laura.longoni@polimi.it>

    @date 20-06-2020


 */


// i: row    index
// j: column index



/*

   _______________________  --> x, j
  |
  |
  |
  |
  |
  |

  y, i


*/

/*

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

int
current_start_chunk(const int& rank, const std::vector<int>& chunk_length_vec)
{
  int result;
  if (rank-1<0)
  {
    result = 0;
  }
  else
  {
    result = chunk_length_vec[rank-1] + current_start_chunk(rank-1, chunk_length_vec);
  }

  return(result);

}

int main()
{
    int size=1; 
    int totSimNumber = 3;
    
    int chunk_length = totSimNumber/size;
    int residual = std::round((double(totSimNumber)/size - chunk_length)*size);
    
    std::vector<int> chunk_sim_vec(size);
    chunk_sim_vec.assign(size, chunk_length);
    
    int rank=0;
    
    cout << chunk_length << " " << current_start_chunk(rank, chunk_sim_vec) << " " << current_start_chunk(rank+1, chunk_sim_vec) << endl;
    
    for (int i = 0; i < residual; i++)
    {
        chunk_sim_vec[i] += 1;
    }
    
    
    for (int ii = min(totSimNumber,current_start_chunk(rank, chunk_sim_vec)+1); 
    ii <= min(totSimNumber,current_start_chunk(rank+1, chunk_sim_vec)); ii++)
    {
        cout << ii << endl;
    }
  
    cout<<"Hello World";

    return 0;
}


*/

#include "utils_H.h"

//! Parse library
#include "GetPot.hpp"


//! for simple profiling
#include "timing.h"

#include <mpi.h>

int
main (int argc, char** argv)
{

  int rank, size;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  tic();
  // Reading parameters through GetPot
  GetPot command_line (argc, argv);
  const std::string dataFileName = command_line.follow ( "SMARTSED_input", 2, "-f", "--file" );
  GetPot dataFile ( dataFileName );


  const Int         totSimNumber           = command_line.follow ( 2, "-sim" );
  const std::string friction_model         = dataFile ( "physics/friction_model", "None" );
  const Real        n_manning              = dataFile ( "physics/n_manning", 0.01 );

  const Real        height_thermometer     = dataFile ( "files/meteo_data/height_thermometer", 200. );


  const UInt        steps_per_hour         = dataFile ( "discretization/steps_per_hour", 10 );
  const Real        max_Days               = dataFile ( "discretization/max_Days", 20 );
  const Real        starting_day           = dataFile ( "discretization/starting_day", 0 );
  const Real        H_min                  = dataFile ( "discretization/H_min", 0.001 );
  const Real        T_thr                  = dataFile ( "discretization/T_thr", 0 );

  const bool        direct_method          = dataFile ( "linear_solver/direct_method", true );

  const bool        save_temporal_sequence = dataFile ( "discretization/save_temporal_sequence", false );
  const bool        isNonReflectingBC      = dataFile ( "discretization/isNonReflectingBC", false );

  const bool        is_sediment_transport  = dataFile ( "physics/is_sediment_transport", true );

  const bool        spit_out_matrix        = dataFile ( "debug/spit_out_matrix", false );
  const std::string matrix_name            = dataFile ( "debug/matrix_name", "/tmp/matrix_" );
  const std::string vector_name            = dataFile ( "debug/vector_name", "/tmp/vector_" );
        std::string tmpname = "";

  const bool        spit_out_solutions_each_time_step = dataFile ( "debug/spit_out_solutions_each_time_step", false );
  
  const Real        frequency_save         = dataFile( "debug/frequency_save", 24. );

  const UInt nstep   = steps_per_hour * max_Days * 24;
  const Real t_final = max_Days * 24 * 3600;
  const Real dt_DSV_given  = t_final / Real ( nstep );
  Real dt_DSV = dt_DSV_given;


  if ((size>totSimNumber && totSimNumber>0) || (size!=1 && totSimNumber<=1))
  {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }


  const int chunk_length = totSimNumber/size;
  const int residual = std::round((double(totSimNumber)/size - chunk_length)*size);

  std::vector<int> chunk_sim_vec(size);
  chunk_sim_vec.assign(size, chunk_length);

  for (UInt i = 0; i < residual; i++)
  {
    chunk_sim_vec[i] += 1;
  }
  toc ("parse command line");

  // execute R script
  //Rscript ../Geostatistics/Downscaling_Simulation_SoilGrids/Downscaling/DownscalingAitchisonSmartSed_2020.R $nsim $res
  if ( totSimNumber >= 0 && rank==0 )
    {
      const std::string bashCommand = std::string ( "Rscript ../Geostatistics/Downscaling_Simulation_SoilGrids/Downscaling/DownscalingAitchisonSmartSed_2020.R " ) + std::to_string( totSimNumber ) + " " + std::to_string ( command_line.follow ( 2, "-scale" ) );
      std::system ( bashCommand.c_str() );
    }
  MPI_Barrier (MPI_COMM_WORLD); 

  if (rank==0) std::cout << "mean # of simulations per rank, " << chunk_length << std::endl;

  for (int currentSimNumber = std::min(totSimNumber, current_start_chunk(rank, chunk_sim_vec)+1); 
    currentSimNumber <= std::min(totSimNumber,current_start_chunk(rank+1, chunk_sim_vec)); currentSimNumber++)
  {

    if (rank==0)
    {
      std::cout << "------------------------ "                            << std::endl;
      std::cout << "friction_model         = " << friction_model          << std::endl;
      std::cout << "n_manning              = " << n_manning               << std::endl;
      std::cout << "steps_per_hour         = " << steps_per_hour          << std::endl;
      std::cout << "max_Days               = " << max_Days                << std::endl;
      std::cout << "t_final                = " << t_final    << " sec."   << std::endl;
      std::cout << "dt_DSV_given           = " << dt_DSV     << " sec."   << std::endl;
      std::cout << "H_min                  = " << H_min                   << std::endl;
      std::cout << "------------------------ "                            << std::endl;
    }



  // +-----------------------------------------------+
  // |                 Reading Input files           |
  // +-----------------------------------------------+

    tic();
    const std::string file_dir       = "../Inputs/";
    const std::string output_dir     = "../Outputs/" + std::to_string ( currentSimNumber ) + "/";

    const std::string orography_file = dataFile ( "files/orography_file", "DEM.txt" );
    const std::string mask_file      = dataFile ( "files/mask_file", "Mask_bin.txt" );


  // precipitation files and temperature
    const std::string temperature_file = dataFile ( "files/meteo_data/temperature_file", "Temperature.txt" );



    const bool restart_H             = dataFile ( "files/initial_conditions/restart_H",             false );
    const bool restart_vel           = dataFile ( "files/initial_conditions/restart_vel",           false );
    const bool restart_snow          = dataFile ( "files/initial_conditions/restart_snow",          false );
    const bool restart_sediment      = dataFile ( "files/initial_conditions/restart_sediment",      false );
    const bool restart_gravitational = dataFile ( "files/initial_conditions/restart_gravitational", false );
    const bool restart_soilMoisture  = dataFile ( "files/initial_conditions/restart_soilMoisture",  false );



    {
      const std::string cmd_str = "mkdir -p " + output_dir;
      char* chararray_cmd = new char[ cmd_str.length() + 1 ];
      const char* cmd_bash = strcpy ( chararray_cmd, cmd_str.c_str() );

      std::system ( cmd_bash );
    }

    const std::string ET_model = dataFile ( "files/evapotranspiration/ET_model", "None" );
    const Real phi_rad         = M_PI / 180. * dataFile ( "files/evapotranspiration/latitude_deg", 45. );

    const std::string infiltrationModel = dataFile ( "files/infiltration/infiltration_model", "None" );


  /* Static Variables */

    UInt N_rows,
    N_cols,
    N;

    std::vector<UInt> idStaggeredBoundaryVectSouth,
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
    idBasinVectReIndex_excluded;


    std::vector<std::array<Real, 2> > Gamma_vect_x,
    Gamma_vect_y;


    std::vector<std::tuple<bool, int> > excluded_ids;
    std::vector<Real> additional_source_term;

    std::vector<Real> basin_mask_Vec,
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
    H;

    Eigen::VectorXd H_basin, rhs;

  Real pixel_size, // meter/pixel
  xllcorner,
  yllcorner,
  xllcorner_staggered_u,
  yllcorner_staggered_u,
  xllcorner_staggered_v,
  yllcorner_staggered_v,
  NODATA_value;

  /*                   */



  

  //
  {
    Raster orographyMat ( file_dir + orography_file );
    Raster basin_mask   ( file_dir + mask_file      );

    if ( basin_mask.cellsize != orographyMat.cellsize )
    {
      std::cout << mask_file << " cellsize and " << orography_file << " cellsize are not equal" << std::endl;
      exit ( -1 );
    }

    pixel_size = Real ( command_line.follow ( 2, "-scale" ) ) * basin_mask.cellsize;

    std::cout << "cell resolution for the current simulation = " << pixel_size << " meters" << std::endl;
    std::cout << "-------------------- "                                                    << std::endl;

    xllcorner    = basin_mask.xllcorner;
    yllcorner    = basin_mask.yllcorner;
    NODATA_value = basin_mask.NODATA_value;
    if ( basin_mask.cellsize <= pixel_size )
    {

      const std::string bashCommand = std::string ( "Rscript -e " ) + "\"library(raster);" +
      "dem=raster('" + file_dir + orography_file + "');" +
      "basin=raster('" + file_dir + mask_file + "');" +
      "basin=aggregate(basin," + std::to_string ( command_line.follow ( 2, "-scale" ) ) + ");" +
      "values(basin)[values(basin)>0]=1;" +
      "dem=resample(dem,basin,method='bilinear');" +
                                        //"values(dem)[is.na(values(dem))]=0;" +
      "writeRaster( dem, file=paste0('" + output_dir + "DEM.asc'), overwrite=TRUE );" +
      "writeRaster( basin, file=paste0('" + output_dir + "basin_mask.asc'), overwrite=TRUE )\"";
      std::system ( bashCommand.c_str() );
    }
    else
    {
      std::cout << "Basin mask greater than simulation resolution, i.e. " << pixel_size << std::endl;
      exit ( -1. );
    }


    if ( dataFile( "discretization/FillSinks", false ) )
    {

      std::string bashCommand = std::string( "python3 -c " ) + "\"import os; import sys; import gdal; cwd = os.getcwd();" +
      "sys.path.append( cwd + '/../DeterministicProgram/include/py' ); import richdem as rd;" +
      "dem = rd.LoadGDAL( '" + output_dir + "DEM.asc' );" +
      "rd.FillDepressions( dem, in_place=True, epsilon=True,topology='D4' );" +
      "rd.SaveGDAL( '" + output_dir + "DEM.tif', dem )\"";

      std::system( bashCommand.c_str() );

      bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
      "dem=raster('" + output_dir + "DEM.tif');" + 
      "writeRaster( dem, file=paste0('" + output_dir + "DEM.asc'), overwrite=TRUE )\"";

      std::system( bashCommand.c_str() );

    }


  }

  {

    Raster orographyMat ( output_dir + "DEM.asc" );
    Raster basin_mask   ( output_dir + "basin_mask.asc" );


    N_rows = basin_mask.Coords.rows();
    N_cols = basin_mask.Coords.cols();
    N      = N_rows * N_cols;


    H                     .resize ( N );
    orography             .resize ( N );
    basin_mask_Vec        .resize ( N );
    eta                   .resize ( N );
    h_G                   .resize ( N );
    h_sd                  .resize ( N );
    h_sn                  .resize ( N );
    S_coeff               .resize ( N );
    W_Gav                 .resize ( N );
    W_Gav_cum             .resize ( N );
    Res_x                 .resize ( N );
    Res_y                 .resize ( N );
    Z_Gav                 .resize ( N );
    d_90                  .resize ( N );
    soilMoistureRetention .resize ( N );
    hydraulic_conductivity.resize ( N );
    roughness_vect        .resize ( N );
    excluded_ids          .resize ( N );
    additional_source_term.resize ( N );
    slope_cell            .resize ( N );

    u             .resize ( ( N_cols + 1 ) * N_rows );
    v             .resize ( ( N_rows + 1 ) * N_cols );
    n_x           .resize ( u.size( ) );
    n_y           .resize ( v.size( ) );
    u_star        .resize ( u.size( ) );
    v_star        .resize ( v.size( ) );

    Gamma_vect_x         .resize ( u.size( ) );
    Gamma_vect_y         .resize ( v.size( ) );

    h_interface_x .resize ( u.size( ) );
    h_interface_y .resize ( v.size( ) );

    slope_x   .resize ( u.size( ) );
    slope_y   .resize ( v.size( ) );




    for ( UInt i = 0; i < N_rows; i++ )
    {
      for ( UInt j = 0; j < N_cols; j++ )
      {
        const auto k = j + i * N_cols;
        basin_mask_Vec[ k ] = basin_mask.Coords.coeff ( i, j ) > 0;
        orography[ k ]      = orographyMat.Coords.coeff ( i, j );
      }
    }

  }

  if ( restart_H )
  {
    {
      const std::string H_file = dataFile ( "files/initial_conditions/H_file", "H.txt" );
      Raster HMat ( file_dir + H_file );

      if ( HMat.cellsize <= pixel_size )
      {
        const std::string bashCommand = std::string ( "Rscript -e " ) + "\"library(raster);" +
        "basin=raster('" + output_dir + "basin_mask.asc" + "');" +
        "H=raster('" + file_dir + H_file + "');" +
        "H=resample(H,basin,method='bilinear');" +
        "values(H)[is.na(values(H))]=0;" +
        "writeRaster( H, file=paste0('" + output_dir + "H_0.asc'), overwrite=TRUE )\"";
        std::system ( bashCommand.c_str() );

      }
      else
      {
        std::cout << "Error! resolution of surface water is greater than simulation resolution, i.e. " << pixel_size << std::endl;
        exit ( -1. );
      }
    }
    Raster HMat ( output_dir + "H_0.asc" );


    for ( UInt i = 0; i < N_rows; i++ )
    {
      for ( UInt j = 0; j < N_cols; j++ )
      {
        const auto k = j + i * N_cols;
        H [ k ]   = HMat.Coords.coeff ( i, j ) * basin_mask_Vec[ k ];
        eta [ k ] = H [ k ] + orography[ k ];
      }
    }

  }
  else
  {
    for ( UInt i = 0; i < N; i++ )
    {
      H   [ i ] = 0.;
      eta [ i ] = 0.;
    }

    saveSolution ( output_dir + "H_0", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, H );
  }


  if ( restart_vel )
  {
    {
      const std::string restart_vel_u_file = dataFile ( "files/initial_conditions/vel_file_u", "restart_vel_u.txt" );
      Raster vel_u_Mat ( file_dir + restart_vel_u_file );

      if ( vel_u_Mat.cellsize <= pixel_size )
      {
        const std::string bashCommand = std::string ( "Rscript -e " ) + "\"library(raster);" +
        "dem=raster('" + output_dir + "DEM.asc" + "');" +
        "u=raster('" + file_dir + restart_vel_u_file + "');" +
        "dem_u=dem;" +
        "extent(dem_u)@xmin=extent(dem)@xmin-res(dem)[[1]]/2;" +
        "extent(dem_u)@xmax=extent(dem)@xmax+res(dem)[[1]]/2;" +
        "res(dem_u)=res(dem);" +
        "extent(u)@xmin=extent(dem)@xmin-res(u)[[1]]/2;" +
        "extent(u)@xmax=extent(dem)@xmax+res(u)[[1]]/2;" +
        "u=resample(u,dem_u,method='bilinear');" +
        "values(u)[is.na(values(u))]=0;" +
        "writeRaster( u, file=paste0('" + output_dir + "u_0.asc'), overwrite=TRUE )\"";
        std::system ( bashCommand.c_str() );

      }
      else
      {
        std::cout << "Error! resolution of horizontal velocity is greater than simulation resolution i.e. " << pixel_size << std::endl;
        exit ( -1. );
      }
    }
    Raster vel_u_Mat ( output_dir + "u_0.asc" );
    for ( UInt i = 0; i < N_rows; i++ )
    {
      for ( UInt j = 0; j <= N_cols; j++ )
      {
              const auto Id = j + ( N_cols + 1 ) * i; // u
              u[ Id ] = vel_u_Mat.Coords.coeff ( i, j );
            }
          }

          xllcorner_staggered_u = vel_u_Mat.xllcorner;
          yllcorner_staggered_u = vel_u_Mat.yllcorner;


          {
            const std::string restart_vel_v_file = dataFile ( "files/initial_conditions/vel_file_v", "restart_vel_v.txt" );
            Raster vel_v_Mat ( file_dir + restart_vel_v_file );

            if ( vel_v_Mat.cellsize <= pixel_size )
            {
              const std::string bashCommand = std::string ( "Rscript -e " ) + "\"library(raster);" +
              "dem=raster('" + output_dir + "DEM.asc" + "');" +
              "v=raster('" + file_dir + restart_vel_v_file + "');" +
              "dem_v=dem;" +
              "extent(dem_v)@ymin=extent(dem)@ymin-res(dem)[[1]]/2;" +
              "extent(dem_v)@ymax=extent(dem)@ymax+res(dem)[[1]]/2;" +
              "res(dem_v)=res(dem);" +
              "extent(v)@ymin=extent(dem)@ymin-res(v)[[1]]/2;" +
              "extent(v)@ymax=extent(dem)@ymax+res(v)[[1]]/2;" +
              "v=resample(v,dem_v,method='bilinear');" +
              "values(v)[is.na(values(v))]=0;" +
              "writeRaster( v, file=paste0('" + output_dir + "v_0.asc'), overwrite=TRUE )\"";
              std::system ( bashCommand.c_str() );

            }
            else
            {
              std::cout << "Error! resolution of vertical velocity is greater than simulation resolution, i.e. " << pixel_size << std::endl;
              exit ( -1. );
            }
          }
          Raster vel_v_Mat ( output_dir + "v_0.asc" );

          xllcorner_staggered_v = vel_v_Mat.xllcorner;
          yllcorner_staggered_v = vel_v_Mat.yllcorner;

          for ( UInt i = 0; i <= N_rows; i++ )
          {
            for ( UInt j = 0; j < N_cols; j++ )
            {
              const auto Id = j + i * N_cols;
              v[ Id ] = vel_v_Mat.Coords.coeff ( i, j );
            }
          }

        }
        else
        {
          xllcorner_staggered_u = xllcorner - pixel_size / 2.;
          yllcorner_staggered_u = yllcorner;

          u.assign (u.size (), 0.0);

          saveSolution ( output_dir + "u_0", "u", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, u );


          xllcorner_staggered_v = xllcorner;
          yllcorner_staggered_v = yllcorner - pixel_size / 2.;

          v.assign (v.size (), 0.0);

          saveSolution ( output_dir + "v_0", "v", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, v );
        }


        if ( restart_snow )
        {
          {
            const std::string restart_snow_file = dataFile ( "files/initial_conditions/snow_file", "h_sn.asc" );
            Raster snow_Mat ( file_dir + restart_snow_file );

            if ( snow_Mat.cellsize <= pixel_size )
            {
              const std::string bashCommand = std::string ( "Rscript -e " ) + "\"library(raster);" +
              "dem=raster('" + output_dir + "DEM.asc" + "');" +
              "hsn=raster('" + file_dir + restart_snow_file + "');" +
              "hsn=resample(hsn,dem,method='bilinear');" +
              "values(hsn)[is.na(values(hsn))]=0;" +
              "writeRaster( hsn, file=paste0('" + output_dir + "hsn_0.asc'), overwrite=TRUE )\"";
              std::system ( bashCommand.c_str() );

            }
            else
            {
              std::cout << "Error! resolution of snow file is greater than simulation resolution, i.e. " << pixel_size << std::endl;
              exit ( -1. );
            }
          }
          Raster snow_Mat ( output_dir + "hsn_0.asc" );

          for ( UInt i = 0; i < N_rows; i++ )
          {
            for ( UInt j = 0; j < N_cols; j++ )
            {
              const auto k = j + i * N_cols;
              h_sn[ k ] = snow_Mat.Coords.coeff ( i, j ) * basin_mask_Vec[ k ];
            }
          }


        }
        else
        {
          h_sn.assign (h_sn.size (), 0.0);
          saveSolution ( output_dir + "hsn_0", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, h_sn );
        }

        if ( restart_sediment )
        {
          {
            const std::string restart_sediment_file = dataFile ( "files/initial_conditions/sediment_file", "h_sd.asc" );
            Raster sediment_Mat ( file_dir + restart_sediment_file );

            if ( sediment_Mat.cellsize <= pixel_size )
            {
              const std::string bashCommand = std::string ( "Rscript -e " ) + "\"library(raster);" +
              "dem=raster('" + output_dir + "DEM.asc" + "');" +
              "hsd=raster('" + file_dir + restart_sediment_file + "');" +
              "hsd=resample(hsd,dem,method='bilinear');" +
              "values(hsd)[is.na(values(hsd))]=0;" +
              "writeRaster( hsd, file=paste0('" + output_dir + "hsd_0.asc'), overwrite=TRUE )\"";
              std::system ( bashCommand.c_str() );

            }
            else
            {
              std::cout << "Error! resolution of sediment file is greater than simulation resolution, i.e. " << pixel_size << std::endl;
              exit ( -1. );
            }
          }
          Raster sediment_Mat ( output_dir + "hsd_0.asc" );


          for ( UInt i = 0; i < N_rows; i++ )
          {
            for ( UInt j = 0; j < N_cols; j++ )
            {
              const auto k = j + i * N_cols;
              h_sd[ k ] = sediment_Mat.Coords.coeff ( i, j ) * basin_mask_Vec[ k ];
            }
          }
        }
        else
        {
          h_sd.assign (h_sd.size (), 0.0);
          saveSolution ( output_dir + "hsd_0", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, h_sd );

        }


        if ( restart_gravitational )
        {
          {
            const std::string restart_gravitational_file = dataFile ( "files/initial_conditions/gravitational_file", "h_G.asc" );
            Raster gravitational_Mat ( file_dir + restart_gravitational_file );

            if ( gravitational_Mat.cellsize <= pixel_size )
            {
              const std::string bashCommand = std::string ( "Rscript -e " ) + "\"library(raster);" +
              "dem=raster('" + output_dir + "DEM.asc" + "');" +
              "hG=raster('" + file_dir + restart_gravitational_file + "');" +
              "hG=resample(hG,dem,method='bilinear');" +
              "values(hG)[is.na(values(hG))]=0;" +
              "writeRaster( hG, file=paste0('" + output_dir + "hG_0.asc'), overwrite=TRUE )\"";
              std::system ( bashCommand.c_str() );

            }
            else
            {
              std::cout << "Error! resolution of gravitational file is greater than simulation resolution, i.e. " << pixel_size << std::endl;
              exit ( -1. );
            }
          }
          Raster gravitational_Mat ( output_dir + "hG_0.asc" );

          for ( UInt i = 0; i < N_rows; i++ )
          {
            for ( UInt j = 0; j < N_cols; j++ )
            {
              const auto k = j + i * N_cols;
              h_G[ k ] = gravitational_Mat.Coords.coeff ( i, j ) * basin_mask_Vec[ k ];
            }
          }

        }
        else
        {
          h_G.assign (h_G.size (), 0.0);
          saveSolution ( output_dir + "hG_0", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, h_G );
        }
        toc ("read input files");


  // +-----------------------------------------------+
  // |    Construct soilMoistureRetention vector     |
  // +-----------------------------------------------+

        tic();

        std::vector<Int> corineCode_Vec ( N ); 
        std::vector<Real> X_Gav ( N ), Y_Gav ( N );

        {
          const std::string corineCode_file = dataFile ( "files/infiltration/corineCode_file", "CLC_RASTER.txt" );

          const std::string check_presence_string = file_dir + corineCode_file; 

          if (!is_file_exist(check_presence_string.c_str()))
          {
            std::cout << check_presence_string << " is not present!" << std::endl;
            exit ( -1. );
          }

    // interpolate CLC to make sure to match correct dimensions
          {
            const std::string bashCommand = std::string ( "Rscript -e " ) + "\"library(raster);" +
            "dem=raster('" + output_dir + "DEM.asc" + "');" +
            "clc=raster('" + file_dir + corineCode_file + "');" +
            "clc=resample(clc,dem,method='ngb');" +
            "values(clc)[is.na(values(clc))]=0;" +
            "writeRaster( clc, file=paste0('" + output_dir + "CLC.asc'), overwrite=TRUE )\"";
            std::system ( bashCommand.c_str() );

          }

          Raster corineCode ( output_dir + "CLC.asc" );


          if ( corineCode.cellsize != pixel_size )
          {
            std::cout << "Please check that the " << corineCode_file << " cellsize is consistent with " << mask_file << " and "     << orography_file  << " ones" << std::endl;
            exit ( -1. );
          }



          for ( UInt i = 0; i < N_rows; i++ )
          {
            for ( UInt j = 0; j < N_cols; j++ )
            {
              const auto k = j + i * N_cols;
              corineCode_Vec[ k ] = corineCode.Coords.coeff ( i, j );
            }
          }

    //saveSolution ( output_dir + "CLC", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, corineCode_Vec );
        }



        {

          std::vector<Int>  HSG ( N );
          std::vector<Real> clayPercentage_Vec ( N ), sandPercentage_Vec ( N );


          std::string str1, str2;
          if ( restart_soilMoisture )
          {
            const std::string restart_clay_file = dataFile ( "files/initial_conditions/clay_file", "clay.asc" ),
            restart_sand_file = dataFile ( "files/initial_conditions/sand_file", "sand.asc" );

            str1 = file_dir + restart_clay_file;
            str2 = file_dir + restart_sand_file;
          }
          else
          {
            str1 = output_dir + "clay_sim_" + std::to_string ( currentSimNumber ) + ".asc";
            str2 = output_dir + "sand_sim_" + std::to_string ( currentSimNumber ) + ".asc";

            if (!is_file_exist(str1.c_str()))
            {
              std::cout << str1 << " is not present! Make sure you have put nsim>0 or in case you want to provide directly the particle size fractions make sure you have put restart_soilMoisture=true and specified the correct paths" << std::endl;
              exit ( -1. );
            }

            if (!is_file_exist(str2.c_str()))
            {
              std::cout << str2 << " is not present! Make sure you have put nsim>0 or in case you want to provide directly the particle size fractions make sure you have put restart_soilMoisture=true and specified the correct paths" << std::endl;
              exit ( -1. );
            }      
          }

          double cellsize_psf = 0;
          if (infiltrationModel != "None" || friction_model == "Rickenmann")
          {

      // interpolate psfs to make sure to match correct dimensions
            if ( restart_soilMoisture )
            {
              std::string bashCommand;

              bashCommand = std::string ( "Rscript -e " ) + "\"library(raster);" +
              "dem=raster('" + output_dir + "DEM.asc" + "');" +
              "clay=raster('" + str1 + "');" +
              "clay=resample(clay,dem,method='ngb');" +
              "values(clay)[is.na(values(clay))]=0;" +
              "writeRaster( clay, file=paste0('" + output_dir + "clay.asc'), overwrite=TRUE )\"";
              std::system ( bashCommand.c_str() );

              bashCommand = std::string ( "Rscript -e " ) + "\"library(raster);" +
              "dem=raster('" + output_dir + "DEM.asc" + "');" +
              "sand=raster('" + str1 + "');" +
              "sand=resample(sand,dem,method='ngb');" +
              "values(sand)[is.na(values(sand))]=0;" +
              "writeRaster( sand, file=paste0('" + output_dir + "sand.asc'), overwrite=TRUE )\"";
              std::system ( bashCommand.c_str() );
            }


            Raster clayPercentage ( str1 ),
            sandPercentage ( str2 );

            cellsize_psf = clayPercentage.cellsize;
            if ( cellsize_psf != sandPercentage.cellsize )
            {
              std::cout << "Please check that the soil texture files have the same resolution" << std::endl;
              exit ( -1. );
            }

            for ( UInt i = 0; i < N_rows; i++ )
            {
              for ( UInt j = 0; j < N_cols; j++ )
              {
                const auto k = j + i*N_cols;

                clayPercentage_Vec[ k ] = clayPercentage.Coords.coeff ( i, j );
                sandPercentage_Vec[ k ] = sandPercentage.Coords.coeff ( i, j );
              }
            }

          }


          if (infiltrationModel != "None" || friction_model == "Rickenmann")
          {


            for ( UInt i = 0; i < N_rows; i++ )
            {
              for ( UInt j = 0; j < N_cols; j++ )
              {


                const auto k = j + i * N_cols;

                const auto & clay = clayPercentage_Vec[ k ],
                & sand = sandPercentage_Vec[ k ];


            if ( sand > .9 && sand <= 1 && clay < .1 && clay >= 0 ) // A
            {
              HSG[ k ] = 0;
            }
            else if ( sand > .5 && sand < .9 && clay > .1 && clay < .2 ) // B
            {
              HSG[ k ] = 1;
            }
            else if ( sand < .5 && sand >= 0 && clay > .2 && clay < .4 ) // C
            {
              HSG[ k ] = 2;
            }
            else if ( sand < .5 && sand >= 0 && clay > .4 && clay <= 1 ) // D
            {
              HSG[ k ] = 3;
            }
            else if ( sand >= 0 && sand <= 1 && clay >= 0 && clay <= 1 )
            {

              Vector2D point ( std::array<Real, 2> {{ clay, sand }} );
              
              
              Vector2D point_A ( std::array<Real, 2> {{ 0,  1  }} );
              Vector2D point_B ( std::array<Real, 2> {{ .1, 1  }} );
              Vector2D point_C ( std::array<Real, 2> {{ .1, .9 }} );
              Vector2D point_D ( std::array<Real, 2> {{ 0, .9  }} );
              
              Vector2D point_E ( std::array<Real, 2> {{ .1, .5 }} );
              Vector2D point_F ( std::array<Real, 2> {{ .2, .5 }} );
              Vector2D point_G ( std::array<Real, 2> {{ .2, .9 }} );
              
              Vector2D point_H ( std::array<Real, 2> {{ .2, 0  }} );
              Vector2D point_I ( std::array<Real, 2> {{ .4, 0  }} );
              Vector2D point_L ( std::array<Real, 2> {{ .4, .5 }} );
              
              Vector2D point_M ( std::array<Real, 2> {{ 1, 0   }} );
              Vector2D point_N ( std::array<Real, 2> {{ 1, .5  }} );
              
              
              std::vector<Vector2D> vv = { point_A,
                point_D,
                point_C,
                point_B,
                
                point_C,
                point_E,
                point_F,
                point_G,
                
                point_F,
                point_H,
                point_I,
                point_L,
                
                point_L,
                point_I,
                point_M,
                point_N
              };
              
              
              std::pair<Real, Int> min = std::make_pair ( 1.e4, -1 );
              
              for ( auto ii = 0; ii < vv.size(); ii += 4 )
              {

                const auto A = vv[ ii ],
                D = vv[ ii + 1 ],
                C = vv[ ii + 2 ],
                B = vv[ ii + 3 ];
                
                Real d1 = 1.e4,
                d2 = 1.e4,
                d3 = 1.e4,
                d4 = 1.e4;
                
                
                Vector2D e1 ( std::array<Real, 2> {{ 1, 0 }} ),
                e2 ( std::array<Real, 2> {{ 0, 1 }} );
                
                
                if ( e1.dot ( point - D ) >= 0 && e1.dot ( point - D ) <= ( C ( 0 ) - D ( 0 ) ) ) d1 = std::abs ( ( point - D ).dot ( e2 ) );
                if ( e2.dot ( point - D ) >= 0 && e2.dot ( point - D ) <= ( A ( 1 ) - D ( 1 ) ) ) d2 = std::abs ( ( point - D ).dot ( e1 ) );
                if ( e1.dot ( point - A ) >= 0 && e1.dot ( point - A ) <= ( B ( 0 ) - A ( 0 ) ) ) d3 = std::abs ( ( point - A ).dot ( e2 ) );
                if ( e2.dot ( point - C ) >= 0 && e2.dot ( point - C ) <= ( B ( 1 ) - C ( 1 ) ) ) d4 = std::abs ( ( point - C ).dot ( e1 ) );
                
                if ( d1 < min.first ) min = std::pair<Real, Int> ( d1, ii / 4. );
                if ( d2 < min.first ) min = std::pair<Real, Int> ( d2, ii / 4. );
                if ( d3 < min.first ) min = std::pair<Real, Int> ( d3, ii / 4. );
                if ( d4 < min.first ) min = std::pair<Real, Int> ( d4, ii / 4. );
                
              }
              
              const auto& id = min.second;
              if ( id == 0 ) // A
              {
                HSG[ k ] = 0;
              }
              else if ( id == 1 ) // B
              {
                HSG[ k ] = 1;
              }
              else if ( id == 2 ) // C
              {
                HSG[ k ] = 2;
              }
              else if ( id == 3 ) // D
              {
                HSG[ k ] = 3;
              }
              else
              {
                std::cout << "Something wrong in HSG classification" << std::endl;
                exit ( -1. );
              }
              
            }
            else
            {
              HSG[ k ] = -1;
            }
            
          }
        }
      }


    const Real S_0 = .254;  // 254 mm
    
    const auto CN_map = createCN_map( );
    
    for ( UInt k = 0; k < N; k++ )
    {
      const auto key = std::array<Int, 2> {{ corineCode_Vec[ k ], HSG[ k ] }};
      
      const auto it = CN_map.find ( key );
      
      
      if ( it != CN_map.end( ) && infiltrationModel != "None" )
      {
        soilMoistureRetention[ k ] = S_0 * ( 100. / Real ( it->second ) - 1. ) * basin_mask_Vec[ k ];
      }
      else
      {
        soilMoistureRetention[ k ] = 0.;
      }
    }
    
    
    // build X_Gav and Y_Gav
    const auto CN_Gav_map = createCN_map_Gav( );
    
    for ( UInt k = 0; k < N; k++ )
    {
      const auto key = corineCode_Vec[ k ];
      
      const auto it = CN_Gav_map.find ( key );
      
      
      if ( it != CN_Gav_map.end( ) ) // remains zero else
      {
        X_Gav[ k ] = it->second[ 0 ];
        Y_Gav[ k ] = it->second[ 1 ];
      }
      
      Z_Gav[ k ] = X_Gav[ k ] * Y_Gav[ k ];
      
    }
    
    
    if ( clayPercentage_Vec[0] == 0 && sandPercentage_Vec[0] == 0 && friction_model == "Rickenmann" )
    {
      std::cout << "clay and sand are both zero (can't compute d90 for friction, maybe change friction_model in SMARTSED_input in Manning if you don't want to run the R script), probably you have not run correctly the Geostatistical preprocessor!, STOP!" << std::endl;
      exit ( 1. );
    }
    
    // build d_10 (for k_c) and d_90 (frictionClass)
    auto d_10 = d_90;
    if (infiltrationModel != "None")
    {
      d_10 = compute_d_perc ( clayPercentage_Vec, sandPercentage_Vec, 10 );
      
      // Equations for hydraulic conductivity estimation from particle size distribution: A dimensional analysis
      // Ji-Peng Wang1, Bertrand François, and Pierre Lambert
      const Real C_H = 6.54e-4;
      const Real gravity = 9.81;
      const Real kin_visc = 0.89e-6;
      for ( UInt i = 0; i < N; i++ )
      {
        hydraulic_conductivity[ i ] = C_H * gravity / kin_visc * std::pow ( d_10[ i ], 2. );
      }
    }
    
    if (infiltrationModel != "None" || friction_model == "Rickenmann") d_90 = compute_d_perc ( clayPercentage_Vec, sandPercentage_Vec, 90 );
    
    saveSolution ( output_dir + "soilMoistureRetention",  " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, soilMoistureRetention );
    saveSolution ( output_dir + "d_10",                   " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, d_10 );
    saveSolution ( output_dir + "d_90",                   " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, d_90 );
    saveSolution ( output_dir + "k_c",                    " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, hydraulic_conductivity );
    
  }

  
  std::cout << "maximum and minimum hydraulic_conductivity  " << *std::max_element ( hydraulic_conductivity.begin(), hydraulic_conductivity.end() ) << " " << *std::min_element ( hydraulic_conductivity.begin(), hydraulic_conductivity.end() ) << std::endl;
  toc ("build soil moisture vector");

  // +-----------------------------------------------+
  // |                Compute Slopes                 |
  // +-----------------------------------------------+


  tic();
  for ( UInt i = 0; i < N_rows; i++ )
  {
    for ( UInt j = 1; j < N_cols; j++ )
    {
      const auto Id = j + i * ( N_cols + 1 );

      slope_x[ Id ] = ( orography[ Id - i ] - orography[ Id - 1 - i ] ) / pixel_size;
      n_x    [ Id ] = - ( orography[ Id - i ] - orography[ Id - 1 - i ] ) / std::abs ( ( orography[ Id - i ] - orography[ Id - 1 - i ] ) );

      if ( std::isnan ( n_x[ Id ] ) ) n_x[ Id ] = 0;
    }
  }

  for ( UInt i = 0, j = 0; i < N_rows; i++ )
  {
    const auto Id = j + i * ( N_cols + 1 );

    slope_x[ Id ] = slope_x[ Id + 1 ];
    n_x    [ Id ] = n_x    [ Id + 1 ];
  }

  for ( UInt i = 0, j = N_cols; i < N_rows; i++ )
  {
    const auto Id = j + i * ( N_cols + 1 );

    slope_x[ Id ] = slope_x[ Id - 1 ];
    n_x    [ Id ] = n_x    [ Id - 1 ];
  }




  for ( UInt i = 1; i < N_rows; i++ )
  {
    for ( UInt j = 0; j < N_cols; j++ )
    {
      const auto Id = j + i * N_cols;

      slope_y[ Id ] = ( orography[ Id ] - orography[ Id - N_cols ] ) / pixel_size;
      n_y    [ Id ] = - ( orography[ Id ] - orography[ Id - N_cols ] ) / std::abs ( ( orography[ Id ] - orography[ Id - N_cols ] ) );

      if ( std::isnan ( n_y[ Id ] ) ) n_y[ Id ] = 0;
    }
  }

  for ( UInt j = 0, i = 0; j < N_cols; j++ )
  {
    const auto Id = j + i * N_cols;

    slope_y[ Id ] = slope_y[ Id + N_cols ];
    n_y    [ Id ] = n_y    [ Id + N_cols ];
  }

  for ( UInt j = 0, i = N_rows; j < N_cols; j++ )
  {
    const auto Id = j + i * N_cols;

    slope_y[ Id ] = slope_y[ Id - N_cols ];
    n_y    [ Id ] = n_y    [ Id - N_cols ];
  }


  saveSolution ( output_dir + "slope_x", "u", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, slope_x );
  saveSolution ( output_dir + "slope_y", "v", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, slope_y );
  toc ("compute slopes");

  // +-----------------------------------------------+
  // |     Compute boundaries of basin domain        |
  // +-----------------------------------------------+


  tic();
  computeAdjacencies ( basin_mask_Vec,
   idStaggeredBoundaryVectSouth,
   idStaggeredBoundaryVectNorth,
   idStaggeredBoundaryVectWest,
   idStaggeredBoundaryVectEast,
   idStaggeredInternalVectHorizontal,
   idStaggeredInternalVectVertical,
   idBasinVect,
   idBasinVectReIndex,
   N_rows,
   N_cols );
  
  
  toc ("compute basin boundaries");

  // +-----------------------------------------------+
  // |                Gavrilovic Coeff.              |
  // +-----------------------------------------------+

  tic();
  for ( const auto & k : idBasinVect )
  {
    const UInt i = k / N_cols;
    slope_cell[ k ] = std::sqrt ( std::pow ( .5 * ( slope_x[ k + i ] + slope_x[ k + i + 1 ] ), 2. ) + std::pow ( .5 * ( slope_y[ k ] + slope_y[ k + N_cols ] ), 2. ) );

    Z_Gav[ k ] *= ( .5 + std::sqrt ( slope_cell[ k ] ) );
    Z_Gav[ k ] = std::pow ( Z_Gav[ k ], 1.5 );
  }
  toc ("gavrilovic coeff");

  // +-----------------------------------------------+
  // |                    Gauges i,j                 |
  // +-----------------------------------------------+

  tic();

  std::vector<std::vector<UInt> > kk_gauges;

  const Int    number_gauges = dataFile ( "discretization/number_gauges", 1 );
  const double delta_gauges  = dataFile ( "discretization/delta_gauges",  0 );


  kk_gauges.resize(number_gauges);

  for (Int number=1; number<=number_gauges; number++)
  {
    std::string filename_x = "discretization/X_gauges_", 
    filename_y = "discretization/Y_gauges_";

    filename_x += std::to_string(number); 
    filename_y += std::to_string(number);

    const Real X_gauges = dataFile ( filename_x.c_str(), 0.);
    const Real Y_gauges = dataFile ( filename_y.c_str(), 0.);

    const Vector2D XX_gauges ( std::array<Real, 2> {{ X_gauges, Y_gauges }} );


    if ( save_temporal_sequence )
    {

      const Vector2D XX_O = std::array<Real, 2> {{ xllcorner, yllcorner + N_rows * pixel_size }};

      auto XX = ( XX_gauges - XX_O )/pixel_size; // coordinate in the matrix

      auto XX_east  = XX + Vector2D(std::array<Real,2>{{delta_gauges/pixel_size,0}});
      auto XX_west  = XX - Vector2D(std::array<Real,2>{{delta_gauges/pixel_size,0}});

      auto XX_south = XX + Vector2D(std::array<Real,2>{{0,delta_gauges/pixel_size}});
      auto XX_north = XX - Vector2D(std::array<Real,2>{{0,delta_gauges/pixel_size}});

      XX(1) = -std::round(XX(1));
      XX(0) =  std::round(XX(0));
      if ( XX(0) < 0 || XX(1) < 0 || XX(1) >= N_rows || XX(0) >= N_cols )
      {
        std::cout << "The gauges in the input file are not good" << std::endl;
        exit ( 1. );
      }


      int i_1 = -std::round(XX_north(1));
      int i_2 = -std::round(XX_south(1));

      int j_1 = std::round(XX_west(0));
      int j_2 = std::round(XX_east(0));

      i_1 = std::min(std::max(i_1, 0), int(N_rows-1));
      j_1 = std::min(std::max(j_1, 0), int(N_cols-1));

      i_2 = std::min(std::max(i_2, 0), int(N_rows-1));
      j_2 = std::min(std::max(j_2, 0), int(N_cols-1));


      for (int i = i_2; i <= i_1; i++)
      {
        for (int j = j_1; j <= j_2; j++)
        {
          kk_gauges[number-1].push_back(i * N_cols + j);
        }
      }

    }

  }
  toc ("gauges i,j ");
  

  // +-----------------------------------------------+
  // |                    Rain                       |
  // +-----------------------------------------------+

  tic();
  const bool is_precipitation = dataFile ( "files/meteo_data/precipitation", true );
  const bool constant_precipitation = dataFile ( "files/meteo_data/constant_precipitation", true );


  Real dt_rain = 0;

  Rain precipitation ( infiltrationModel, N, dataFile ( "files/infiltration/isInitialLoss", false ), dataFile ( "files/infiltration/perc_initialLoss", 0.05 ) );

  if ( constant_precipitation || !is_precipitation )
    {
      const std::string precipitation_file = dataFile ( "files/meteo_data/rain_file", " " );

      const Real time_spacing_rain  = dataFile ( "files/meteo_data/time_spacing_rain", 1. );

      dt_rain = time_spacing_rain * 3600;

      const auto ndata_rain = std::round ( max_Days * 24 / time_spacing_rain );

      precipitation.constant_precipitation ( file_dir + precipitation_file, ndata_rain, is_precipitation, time_spacing_rain );
    }
  else // IDW
    {
      std::vector<std::string> precipitation_file;
      std::vector<Real>        time_spacing_rain, X, Y;
      std::vector<UInt>        ndata_rain;

      const Int number_stations = dataFile ( "files/meteo_data/number_stations", 1 );
      for (Int number=1; number<=number_stations; number++)
      {
        std::string filename = "files/meteo_data/rain_file_";
        filename += std::to_string(number);
        const std::string precipitation_file_current = dataFile ( filename.c_str(), " " );
        
        filename = "files/meteo_data/time_spacing_rain_";
        filename += std::to_string(number);
        const Real time_spacing_rain_current  = dataFile ( filename.c_str(),  1. );

        filename = "files/meteo_data/X_";
        filename += std::to_string(number); 
        const Real X_current = dataFile ( filename.c_str(), 1. );
        
        filename = "files/meteo_data/Y_";
        filename += std::to_string(number);
        const Real Y_current = dataFile ( filename.c_str(), 1. );

        const UInt ndata_rain_current = std::round ( max_Days * 24 / time_spacing_rain_current );

        precipitation_file.push_back(file_dir + precipitation_file_current);
        time_spacing_rain .push_back(time_spacing_rain_current);
        X                 .push_back(X_current);
        Y                 .push_back(Y_current);
        ndata_rain        .push_back(ndata_rain_current);
      }


      dt_rain = *std::min_element ( time_spacing_rain.begin(), time_spacing_rain.end() ) * 3600;

      precipitation.IDW_precipitation ( precipitation_file, ndata_rain, time_spacing_rain, X, Y, xllcorner, yllcorner, pixel_size, N_rows, N_cols, idBasinVect );
    }
  toc ("rain");

  // +-----------------------------------------------+
  // |                 Temperature                   |
  // +-----------------------------------------------+

  tic();
  const Real        time_spacing_temp = dataFile ( "files/meteo_data/time_spacing_temp", 1. );
  const std::string format_temp       = dataFile ( "files/meteo_data/format_temp", "arpa" );
  const Real        dt_temp           = time_spacing_temp * 3600;

  Temperature temp ( file_dir + temperature_file,
                     N,
                     max_Days,
                     T_thr,
                     orography,
                     std::round ( max_Days * 24 / time_spacing_temp ),
                     steps_per_hour,
                     time_spacing_temp,
                     height_thermometer,
                     format_temp );

  toc ("temperature");

  // +-----------------------------------------------+
  // |              Evapotranspiration               |
  // +-----------------------------------------------+

  tic();
  evapoTranspiration ET ( ET_model, N, orography, temp.J, max_Days, phi_rad, height_thermometer );
  toc ("evapotranspiration");

  // +-----------------------------------------------+
  // |                   Core Part                   |
  // +-----------------------------------------------+


  tic();
  
  for (int ii = 0; ii < N; ii++)
  {
    const auto r1 =  dataFile ( "files/infiltration/roughness_scale_factor1", 100. );
    const auto r2 =  dataFile ( "files/infiltration/roughness_scale_factor2", 100. );
    const auto r3 =  dataFile ( "files/infiltration/roughness_scale_factor3", 100. );
    
    
    if (slope_cell[ ii ] <= 0.2)
    {
      roughness_vect[ ii ] = r1;
    }
    else if (slope_cell[ ii ] <= 0.6)
    {
      roughness_vect[ ii ] = r2;
    }
    else
    {
      roughness_vect[ ii ] = r3;
    }
  }



  const Real dt_min = std::min ( dt_rain, dt_temp );
  for ( const auto& kk : hydraulic_conductivity )
    {
      if ( kk * ( dt_min / pixel_size ) > 1. )
        {
          std::cout << "Error! Courant number for gravitational layer is greater than 1!" << std::endl;
          exit ( -1. );
        }
    }

  const Real c1_min = dt_min / pixel_size;




  const Real g = 9.81;
  std::function<Real(Real const&, Real const&)> c1_DSV = [](Real const& dt_DSV, Real const& pixel_size){return dt_DSV / pixel_size;};
  std::function<Real(Real const&, Real const&)> c2_DSV = [](Real const& g,      Real const& c1        ){return g * c1;             };
  std::function<Real(Real const&, Real const&)> c3_DSV = [](Real const& g,      Real const& c1        ){return g * c1 * c1;        };
  

  const Real area = std::pow ( pixel_size, 2 ) * 1.e-6; // km^2

//  std::cout << "c1_DSV " << c1_DSV(dt_DSV, pixel_size) << "   c2_DSV " << c2_DSV(g, c1_DSV(dt_DSV, pixel_size)) << "   c3_DSV " << c3_DSV(g, c1_DSV(dt_DSV, pixel_size)) << std::endl;


  Real slope_y_max = 0.;
  Real slope_x_max = 0.;  
  for ( const auto& k : idBasinVect )
    {
      const UInt i = k / N_cols;

      const auto slope_x_l = std::abs(slope_x[ k + i      ]);
      const auto slope_x_r = std::abs(slope_x[ k + i + 1  ]);
      const auto slope_y_l = std::abs(slope_y[ k          ]);
      const auto slope_y_r = std::abs(slope_y[ k + N_cols ]);

      if (slope_x_l > slope_x_max) slope_x_max = slope_x_l;
      if (slope_x_r > slope_x_max) slope_x_max = slope_x_r;
      if (slope_y_l > slope_y_max) slope_y_max = slope_y_l;
      if (slope_y_r > slope_y_max) slope_y_max = slope_y_r;
    }


  Real dt_sed,
       c1_sed;

  UInt numberOfSteps = 1;

  bool isHNegative =  false;
  
  

  toc ("prepare loop");
  
  
  
  // +-----------------------------------------------+
  // |              compute_sub_basins               |
  // +-----------------------------------------------+
  
  excluded_ids.assign(N, std::tuple<bool,int>( false, -1 ));
  additional_source_term.assign(N, 0.);
  
  const bool static_subbasin_approx = dataFile ( "discretization/static_subbasin_approx", false );
  if (static_subbasin_approx)
  {
  
    tic();
    
    const std::set<UInt> idBasinVect_set                 (idBasinVect.begin(),                  idBasinVect.end()),
                         idStaggeredBoundaryVectSouth_set(idStaggeredBoundaryVectSouth.begin(), idStaggeredBoundaryVectSouth.end()),
                         idStaggeredBoundaryVectNorth_set(idStaggeredBoundaryVectNorth.begin(), idStaggeredBoundaryVectNorth.end()),
                         idStaggeredBoundaryVectWest_set (idStaggeredBoundaryVectWest.begin(),  idStaggeredBoundaryVectWest.end()),
                         idStaggeredBoundaryVectEast_set (idStaggeredBoundaryVectEast.begin(),  idStaggeredBoundaryVectEast.end());
    
    const Real slope_thr = dataFile("discretization/slope_thr", 1.);
    
    
    
    // exclude high slopes,
    for ( const UInt& Id : idBasinVect )
    {
      
      const auto & current_slope_cell = slope_cell[ Id ];
      if (current_slope_cell > slope_thr)
      {
        std::get<0>( excluded_ids[ Id ] ) = true;
      }
    }
    
    // exclude also isolated cells, see below,
    // maybe cycle on intefaces, build a list of bool for each id
    std::vector<UInt> counter_near_excl;
    counter_near_excl.resize(N);
    counter_near_excl.assign(N, 0);
    
    // cycle over interfaces, internal vertical and horizontal
    for ( const auto & Id : idStaggeredInternalVectHorizontal )
    {
      const UInt i = Id / ( N_cols + 1 ),
      
      IDleft  = Id - i - 1, // H
      IDright = Id - i;
      
      if ( std::get<0>(excluded_ids[ IDleft ]) )
      {
        counter_near_excl[ IDright ] += 1;
      }
      if ( std::get<0>(excluded_ids[ IDright ]) )
      {
        counter_near_excl[ IDleft ] += 1;
      }
      
    }
    
    
    for ( const auto & Id : idStaggeredInternalVectVertical )
    {
      
      const UInt IDleft  = Id - N_cols, // H
      IDright = Id;
      
      if ( std::get<0>(excluded_ids[ IDleft ]) )
      {
        counter_near_excl[ IDright ] += 1;
      }
      if ( std::get<0>(excluded_ids[ IDright ]) )
      {
        counter_near_excl[ IDleft ] += 1;
      }
      
    }

    for ( const auto & Id : idBasinVect )
    {
      if ( counter_near_excl[ Id ] == 4 ) // it is surely an isolated cell! (can be also an excluded cell but no problem..)
      {
        std::get<0>( excluded_ids[ Id ] ) = true;
      }
    }
    

    // compute pour points, local movement in the 8 directions! (also diagonal ones)
    for ( const UInt& Id : idBasinVect )
    {

      auto & current_tuple = excluded_ids[ Id ];

      const auto & current_is      = std::get<0>( current_tuple );
            auto & current_pour_id = std::get<1>( current_tuple );

      if ( current_is )
      {
        // start from Id and perform a local search (by means of the gradient) to get the pour point
        int candidate_pour_id = Id;

        while (true)
        {

          candidate_pour_id = computePourCell(candidate_pour_id,
                                              N_cols,
                                              orography,
                                              idBasinVect_set,
                                              idStaggeredBoundaryVectSouth_set,
                                              idStaggeredBoundaryVectNorth_set,
                                              idStaggeredBoundaryVectWest_set,
                                              idStaggeredBoundaryVectEast_set);



          // if the gradient points outside the basin leave current_pour_id = -1
          if (candidate_pour_id < 0) break;


          const auto& is_current_id_again_excluded = std::get<0> ( excluded_ids[ candidate_pour_id ] );
          if ( !is_current_id_again_excluded )
          {
            current_pour_id = candidate_pour_id;
            break;
          }


        }

      }
    }
    

    toc ("get sub-basins");
    
  }
  
  
  computeAdjacencies ( basin_mask_Vec,
                      excluded_ids,
                      idStaggeredBoundaryVectSouth_excluded,
                      idStaggeredBoundaryVectNorth_excluded,
                      idStaggeredBoundaryVectWest_excluded,
                      idStaggeredBoundaryVectEast_excluded,
                      idStaggeredInternalVectHorizontal_excluded,
                      idStaggeredInternalVectVertical_excluded,
                      idBasinVect_excluded,
                      idBasinVectReIndex_excluded,
                      N_rows,
                      N_cols );

  upwind H_interface ( H, u, v, idStaggeredInternalVectHorizontal_excluded,
                                idStaggeredBoundaryVectWest_excluded,
                                idStaggeredBoundaryVectEast_excluded, 
                                idStaggeredInternalVectVertical_excluded,
                                idStaggeredBoundaryVectNorth_excluded,
                                idStaggeredBoundaryVectSouth_excluded, N_rows, N_cols );

  frictionClass alfa ( H_interface.horizontal, H_interface.vertical, u, v, 
                       idStaggeredInternalVectHorizontal_excluded,
                       idStaggeredBoundaryVectWest_excluded,
                       idStaggeredBoundaryVectEast_excluded,
                       idStaggeredInternalVectVertical_excluded,
                       idStaggeredBoundaryVectNorth_excluded,
                       idStaggeredBoundaryVectSouth_excluded, friction_model, n_manning, dt_DSV, d_90, roughness_vect, 0., N_rows, N_cols, slope_x, slope_y );
  
  
  H_basin.resize ( idBasinVect_excluded.size() );
  rhs    .resize ( idBasinVect_excluded.size() );
  

  for (int i = 0; i < H_basin.size(); i++) H_basin( i ) = 0.;
  for (int i = 0; i < rhs.size(); i++) rhs( i ) = 0.;
   
  saveSolution (output_dir + "excluded_ids", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, excluded_ids);
  
  
  SpMat A ( idBasinVect_excluded.size(), idBasinVect_excluded.size() );
  
  // row, column and value in the Triplet
  std::vector<Eigen::Triplet<Real> > coefficients;
  
  /*
  coefficients.reserve ( idBasinVect_excluded.size() +
                        4 * idStaggeredInternalVectHorizontal_excluded.size() +
                        4 * idStaggeredInternalVectVertical_excluded.size() );*/
  
  
  int iter = 0;
  
  
  double c1_DSV_, c2_DSV_, c3_DSV_, minH, maxH;

  maxH = *std::max_element( H.begin(), H.end() );
   
  dt_DSV = maxdt(u, v, g, maxH, pixel_size);
  dt_DSV = dt_DSV < dt_DSV_given ? dt_DSV : dt_DSV_given;
  dt_DSV = dt_DSV < t_final ? dt_DSV : t_final;
      
  c1_DSV_ = c1_DSV (dt_DSV, pixel_size);
  c2_DSV_ = c2_DSV (g, c1_DSV_);
  c3_DSV_ = c3_DSV (g, c1_DSV_);

  //omp_set_num_threads(omp_get_num_procs());

  //std::cout << "# of available threads, " << omp_get_num_procs() << std::endl; 

  auto H_old    = H;
  auto H_oldold = H;
  
  
  double time = 0., timed = -dt_DSV, timedd = -2.*dt_DSV; 
  bool is_last_step = false, check_last = false;
  while ( !is_last_step )
    {

      std::cout << "Simulation progress: " << time / t_final * 100 << " %" << " max surface run-off vel. based Courant: " << maxCourant ( u, v, c1_DSV_ ) << " max surface run-off cel. based Courant: " << maxCourant ( H, g, c1_DSV_ ) << std::endl;
      std::cout << "Current dt, " << dt_DSV << ", given dt, " << dt_DSV_given << std::endl;

      tic();
      // Compute interface fluxes via upwind method
      H_interface.computeHorizontal ( );
      H_interface.computeVertical   ( );
      toc ("H_interface");

      tic();
      // Compute alfa coefficients
      alfa.f_x ( );
      alfa.f_y ( );
      toc ("alfa");

      
      // update only if necessary  --> governed by temperature dynamics, i.e. time_spacing_temp
      if ( std::floor ( time / dt_temp ) > std::floor ( (time - dt_DSV) / dt_temp ) )
        {
          tic();
          // Compute temperature map
          temp.computeTemperature ( std::floor( time / dt_temp ), orography, idBasinVect );
          toc ("temperature");
        }
      

      
      // ET varies daily
      if ( std::floor ( time / (24. * 3600) ) > std::floor ( (time - dt_DSV) / (24. * 3600) ) )
        {
          tic();
          // Get ET rate at the current time
          ET.ET ( temp.T_dailyMean, temp.T_dailyMin, temp.T_dailyMax, std::floor( time / ( 24 * 3600 ) ), idBasinVect, orography );
          toc ("ET");
        }
      

      
      // update only if necessary
      if ( std::floor ( time / dt_min ) > std::floor ( (time - dt_DSV) / dt_min ) )
        {
          tic();
          // +-----------------------------------------------+
          // |            Gravitational Layer                |
          // +-----------------------------------------------+

          // h_G
          // vertical and horizontal residuals for Gravitational Layer
          computeResiduals ( n_x,
                             n_y,
                             N_cols,
                             N_rows,
                             h_G,
                             hydraulic_conductivity,
                             idStaggeredInternalVectHorizontal,
                             idStaggeredInternalVectVertical,
                             idStaggeredBoundaryVectWest,
                             idStaggeredBoundaryVectEast,
                             idStaggeredBoundaryVectNorth,
                             idStaggeredBoundaryVectSouth,
                             idBasinVect,
                             h_interface_x,
                             h_interface_y,
                             Res_x,
                             Res_y );
          toc ("gravit") ;
        }
      

      tic();
      //
      precipitation.computePrecipitation (time,
                                          soilMoistureRetention,
                                          temp.melt_mask,
                                          h_G,
                                          H,
                                          N_rows,
                                          N_cols,
                                          idBasinVect);
      toc ("precipitation");

      
      // update only if necessary
      if ( std::floor ( time / dt_min ) > std::floor ( (time - dt_DSV) / dt_min ) )
        {
          tic();
          for ( const UInt& k : idBasinVect )
            {
              S_coeff[ k ] = 4.62e-10 * h_sn[ k ] * ( temp.T_raster[ k ] - T_thr ) * temp.melt_mask[ k ];
            }

          for ( const auto& k : idBasinVect )
            {
              h_G[ k ] += ( S_coeff[ k ] - ET.ET_vec[ k ] ) * dt_min + precipitation.DP_infiltrated[ k ] * dt_min - c1_min * ( Res_x[ k ] + Res_y[ k ] );
              h_G[ k ] *= ( h_G[ k ] >= 0 ); // to account for evapotranspiration
            }

          // +-----------------------------------------------+
          // |              Snow Accumulation                |
          // +-----------------------------------------------+


          for ( const auto& k : idBasinVect )
            {
              const auto snow_acc = precipitation.DP_total[ k ] * ( 1. - temp.melt_mask[ k ] ) * dt_min - S_coeff[ k ] * dt_min;
              h_sn[ k ] += snow_acc;
            }
          toc ("snow");
        }
      


      // +-----------------------------------------------+
      // |               De Saint Venant                 |
      // +-----------------------------------------------+


      tic();
      // fill u_star and v_star with a Bilinear Interpolation
      bilinearInterpolation ( u,
                              v,
                              u_star,
                              v_star,
                              N_rows,
                              N_cols,
                              dt_DSV,
                              pixel_size,
                              idStaggeredInternalVectHorizontal_excluded,
                              idStaggeredInternalVectVertical_excluded,
                              idStaggeredBoundaryVectWest_excluded,
                              idStaggeredBoundaryVectEast_excluded,
                              idStaggeredBoundaryVectNorth_excluded,
                              idStaggeredBoundaryVectSouth_excluded );

      
      toc ("bilinearInterpolation");



      tic();

      coefficients.reserve ( idBasinVect_excluded.size() +
                        4 * idStaggeredInternalVectHorizontal_excluded.size() +
                        4 * idStaggeredInternalVectVertical_excluded.size() );

      additional_source_term.assign(N, 0.);
      
      buildMatrix ( H_interface.horizontal,
                   H_interface.vertical,
                   orography,
                   u_star,
                   v_star,
                   u,
                   v,
                   H,
                   N_cols,
                   N_rows,
                   N,
                   c1_DSV_,
                   c3_DSV_,
                   0, // 0 
                   precipitation.DP_cumulative,
                   dt_DSV,  
                   alfa.alfa_x,
                   alfa.alfa_y,
                   idStaggeredInternalVectHorizontal_excluded,
                   idStaggeredInternalVectVertical_excluded,
                   idStaggeredBoundaryVectWest_excluded,
                   idStaggeredBoundaryVectEast_excluded,
                   idStaggeredBoundaryVectNorth_excluded,
                   idStaggeredBoundaryVectSouth_excluded,
                   idBasinVect_excluded,
                   idBasinVect,
                   idStaggeredInternalVectHorizontal,
                   idStaggeredInternalVectVertical,
                   idBasinVectReIndex_excluded,
                   isNonReflectingBC,
                   true,
                   excluded_ids,
                   additional_source_term,
                   coefficients,
                   rhs );
      
      
      



      A.setFromTriplets ( coefficients.begin(), coefficients.end() );
      A.makeCompressed();
      coefficients.clear();

      toc ("assemble matrix");

      tic();
      if (spit_out_matrix)
        {
          tmpname = matrix_name + std::to_string (iter);
          std::cout << "saving " << tmpname << std::endl;
          Eigen::saveMarket (A, tmpname, false);
          tmpname = vector_name + std::to_string (iter);
          std::cout << "saving " << tmpname << std::endl;
          Eigen::saveMarketVector_lis (rhs, tmpname);
        }
      toc ("spit out matrix for debug");


      tic();
      if ( direct_method )  // Direct Sparse method: Cholesky being A spd
        {
          Eigen::SimplicialLDLT<SpMat, Eigen::Upper> solver;   


          solver.compute ( A );

          if ( solver.info() != Eigen::Success )
            {
              std::cout << "Decomposition Failed" << std::endl;
              exit ( -1. );
            }

          H_basin = solver.solve ( rhs );

        }
      else
        {

          double tol = 1.e-6;
          int result, maxit = 15000;
          
          /*
          Eigen::ConjugateGradient < SpMat, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double> > cg;
          cg.compute ( A );


          H_basin = cg.solveWithGuess ( rhs, H_basin );

          std::cout << "# threads used by Eigen, " << Eigen::nbThreads( ) << std::endl;
          std::cout << "# iterations:    "         << cg.iterations()     << std::endl;
          std::cout << "estimated error: "         << cg.error()          << std::endl;
          */
          
          
          // Now with IML++
          Eigen::IncompleteCholesky<double> IC(A);// Create I cholesky preconditioner
          //Eigen::IncompleteLUT<double> ILU(A); // create ILU preconditioner
          result = LinearAlgebra::CG(A, H_basin, rhs, IC, maxit, tol);   // Solve system
          
          std::cout <<" IML++ CG "<< std::endl;
          std::cout << "CG flag = " << result << std::endl;
          std::cout << "iterations performed: " << maxit << std::endl;
          std::cout << "tolerance achieved  : " << tol << std::endl;
           

        }
    
      
      
      minH = H_basin.minCoeff(); maxH = H_basin.maxCoeff();
      std::cout << "min H: " << minH << " max H: " << maxH << std::endl;
      
      
      
      if (minH < -H_min) 
      {
        
        dt_DSV = dt_DSV/10.;
        
        if (dt_DSV == 0)
        {
          std::cout << "dt has gone to zero, sorry, STOP!" << std::endl;
          exit( -1. );
        }

        c1_DSV_ = c1_DSV (dt_DSV, pixel_size);
        c2_DSV_ = c2_DSV (g, c1_DSV_);
        c3_DSV_ = c3_DSV (g, c1_DSV_);

        continue;

      }
      

      for ( const UInt& Id : idBasinVect_excluded )
        {
          const UInt IDreIndex = idBasinVectReIndex_excluded[ Id ];

          H_oldold[ Id ] = H_old[ Id ];
          H_old[ Id ] = H[ Id ];

          H [ Id ] = std::abs(H_basin ( IDreIndex ));
          eta [ Id ] = H [ Id ] + orography[ Id ];
        }



      // +-----------------------------------------------+
      // |                  Update u, v                  |
      // +-----------------------------------------------+
      
      
      updateVel ( u,
                  v,
                  u_star,
                  v_star,
                  alfa.alfa_x,
                  alfa.alfa_y,
                  N_rows,
                  N_cols,
                  c2_DSV_,
                  0, // 0
                  eta,
                  H,
                  orography,
                  idStaggeredInternalVectHorizontal_excluded,
                  idStaggeredInternalVectVertical_excluded,
                  idStaggeredBoundaryVectWest_excluded,
                  idStaggeredBoundaryVectEast_excluded,
                  idStaggeredBoundaryVectNorth_excluded,
                  idStaggeredBoundaryVectSouth_excluded,
                  isNonReflectingBC );
      
      
      
      /*
      for ( const UInt& k : idBasinVect )
        {
          if ( H ( k ) < 0 )
            {
              H ( k ) = 0.;
              isHNegative = true;
              
              eta ( k ) = H ( k ) + orography[ k ];
            }
        }*/
      toc ("solve");

      


      // +-----------------------------------------------+
      // |             Sediment Transport                |
      // +-----------------------------------------------+

      if (is_sediment_transport)
      {
        tic();
        const double alfa_coeff = 2.5, beta_coeff = 1.6, gamma_coeff = 1.;

        dt_sed = compute_dt_sediment ( alfa_coeff, beta_coeff, slope_x_max, slope_y_max, u, v, pixel_size, dt_DSV, numberOfSteps );
        c1_sed = dt_sed / pixel_size;



        // vertical and horizontal residuals truncated for Sediment Transport
        computeResidualsTruncated ( u,
                                    v,
                                    N_cols,
                                    N_rows,
                                    N,
                                    c1_sed,
                                    slope_x,
                                    slope_y,
                                    alfa_coeff,   // alfa
                                    beta_coeff,   // beta
                                    gamma_coeff,    // gamma
                                    idStaggeredInternalVectHorizontal_excluded,
                                    idStaggeredInternalVectVertical_excluded,
                                    idStaggeredBoundaryVectWest_excluded,
                                    idStaggeredBoundaryVectEast_excluded,
                                    idStaggeredBoundaryVectNorth_excluded,
                                    idStaggeredBoundaryVectSouth_excluded,
                                    Gamma_vect_x,
                                    Gamma_vect_y
                                    );


        additional_source_term.assign(N, 0.);
        
        for ( const auto & k : idBasinVect )
        {
          const auto & current_tuple = excluded_ids[ k ];
          if ( std::get<0>( current_tuple ) )
          {
            const auto & k_pour = std::get<1>( current_tuple );
            if ( k_pour >= 0 )
            {
              additional_source_term[ k_pour ] += 1.e-3 * M_PI * Z_Gav[ k ] * std::sqrt ( std::abs ( ( .1 + .1 * temp.T_raster[ k ] ) * temp.melt_mask[ k ] ) ) * precipitation.DP_total[ k ] * dt_sed;
            }
          }
        }


        for ( const UInt& k : idBasinVect_excluded )
        {
            // 1.e-3 is the conversion factor, look at EPM theory
          W_Gav[ k ] = 1.e-3 * M_PI * Z_Gav[ k ] * std::sqrt ( std::abs ( ( .1 + .1 * temp.T_raster[ k ] ) * temp.melt_mask[ k ] ) ) * precipitation.DP_total[ k ] * dt_sed + additional_source_term[ k ];
        }

        std::cout << "# steps for solid transport, " << numberOfSteps << std::endl;

        for ( UInt kk = 0; kk < numberOfSteps; kk++ )
        {

          additional_source_term.assign(N, 0.);

          // horizontal
          for ( const UInt& Id : idStaggeredInternalVectHorizontal_excluded )
          {
            const UInt i       = Id / ( N_cols + 1 ),
            IDeast  = Id - i,
            IDwest  = Id - i - 1;

            const Real& h_left  = h_sd[ IDwest ],
            & h_right = h_sd[ IDeast ];


            h_interface_x[ Id ] = Gamma_vect_x[ Id ][ 0 ] * h_right +
                                  Gamma_vect_x[ Id ][ 1 ] * h_left;
          }


          for ( const UInt& Id : idStaggeredBoundaryVectWest_excluded )
          {
            const UInt i = Id / ( N_cols + 1 );

            const Real h_left  = 0, //0,
            h_right = h_sd[ Id - i ];


            h_interface_x[ Id ] = Gamma_vect_x[ Id ][ 0 ] * h_right +
                                  Gamma_vect_x[ Id ][ 1 ] * h_left;
          }


          for ( const UInt& Id : idStaggeredBoundaryVectEast_excluded )
          {

            const UInt i = Id / ( N_cols + 1 );


            const Real h_left  = h_sd[ Id - i - 1 ],
            h_right = 0;


            h_interface_x[ Id ] = Gamma_vect_x[ Id ][ 0 ] * h_right +
                                  Gamma_vect_x[ Id ][ 1 ] * h_left;
          }

          for ( const UInt& Id : idBasinVect_excluded )
          {
            const UInt i = Id / N_cols;

            Res_x[ Id ] = h_interface_x[ Id + 1 + i ] - h_interface_x[ Id + i ];
          }






          // vertical
          for ( const UInt& Id : idStaggeredInternalVectVertical_excluded )
          {
            const UInt IDsouth = Id,
            IDnorth = Id - N_cols;

            const Real& h_left  = h_sd[ IDnorth ],
            & h_right = h_sd[ IDsouth ];


            h_interface_y[ Id ] = Gamma_vect_y[ Id ][ 0 ] * h_right +
                                  Gamma_vect_y[ Id ][ 1 ] * h_left;
          }


          for ( const UInt& Id : idStaggeredBoundaryVectNorth_excluded )
          {

            const Real h_left  = 0, // 0
            h_right = h_sd[ Id ];


            h_interface_y[ Id ] = Gamma_vect_y[ Id ][ 0 ] * h_right +
                                  Gamma_vect_y[ Id ][ 1 ] * h_left;

          }



          for ( const UInt& Id : idStaggeredBoundaryVectSouth_excluded )
          {

            const Real h_left  = h_sd[ Id - N_cols ],
            h_right = 0;

            h_interface_y[ Id ] = Gamma_vect_y[ Id ][ 0 ] * h_right +
                                  Gamma_vect_y[ Id ][ 1 ] * h_left;
          }

          for ( const UInt& Id : idBasinVect_excluded )
          {
            Res_y[ Id ] = h_interface_y[ Id + N_cols ] - h_interface_y[ Id ];
          }



          for ( const auto & Id : idStaggeredInternalVectHorizontal )
          {

            const UInt i       = Id / ( N_cols + 1 ),

            IDleft  = Id - i - 1,
            IDright = Id - i;



            if ( std::get<0> ( excluded_ids[ IDleft ] ) )
            {
              const auto& k_pour = std::get<1> ( excluded_ids[ IDleft ] );
              if ( k_pour >= 0 )
              {
                additional_source_term[ k_pour ] += h_interface_x[ Id ] * std::abs (u[ Id ]) * c1_sed;
              }
            }

            if ( std::get<0> ( excluded_ids[ IDright ] ) )
            {
              const auto& k_pour = std::get<1> ( excluded_ids[ IDright ] );
              if ( k_pour >= 0 )
              {
                additional_source_term[ k_pour ] += h_interface_x[ Id ] * std::abs (u[ Id ]) * c1_sed;
              }
            }



          }


          for ( const auto & Id : idStaggeredInternalVectVertical )
          {

            const UInt IDleft  = Id - N_cols,
            IDright = Id;


            if ( std::get<0> ( excluded_ids[ IDleft ] ) )
            {
              const auto& k_pour = std::get<1> ( excluded_ids[ IDleft ] );
              if ( k_pour >= 0 )
              {
                additional_source_term[ k_pour ] += h_interface_y[ Id ] * std::abs (v[ Id ]) * c1_sed;
              }
            }

            if ( std::get<0> ( excluded_ids[ IDright ] ) )
            {
              const auto& k_pour = std::get<1> ( excluded_ids[ IDright ] );
              if ( k_pour >= 0 )
              {
                additional_source_term[ k_pour ] += h_interface_y[ Id ] * std::abs (v[ Id ]) * c1_sed;
              }
            }

          }



          for ( const UInt& Id : idBasinVect_excluded )
          {
            h_sd[ Id ] += - ( Res_x[ Id ] + Res_y[ Id ] ) + W_Gav[ Id ] + additional_source_term[ Id ];
          }


        }
        toc ("advance");
      }
      

      

      tic();

      // +-----------------------------------------------+
      // |               Update time                     |
      // +-----------------------------------------------+

      timedd = timed;
      timed = time;
      time += dt_DSV;

      if (check_last)
      {
        time = t_final;
        is_last_step = true;
      }
      
      
      // +-----------------------------------------------+
      // |             Save The Raster Solution          |
      // +-----------------------------------------------+

      
      if ( save_temporal_sequence )
      {

        for (Int number=1; number<=number_gauges; number++)
        {
          std::string filename_x = "discretization/X_gauges_", 
                      filename_y = "discretization/Y_gauges_";

          filename_x += std::to_string(number); 
          filename_y += std::to_string(number);

          const Real X_gauges = dataFile ( filename_x.c_str(), 0.);
          const Real Y_gauges = dataFile ( filename_y.c_str(), 0.);

          const Vector2D XX_gauges ( std::array<Real, 2> {{ X_gauges, Y_gauges }} );

          const Vector2D XX_O = std::array<Real, 2> {{ xllcorner, yllcorner + N_rows * pixel_size }};

          auto XX = ( XX_gauges - XX_O )/pixel_size; // coordinate in the matrix

          
          auto H_candidate = bilinearInterpolation (u, v, H, N_cols, N_rows, XX);


          double H_current = 0.;
          UInt kk_gauges_max = 0;
          for (const auto & candidate : kk_gauges[number-1])
          {
            const auto & cc = H[ candidate ];
            if (cc > H_current)
            {
              kk_gauges_max = candidate;
              H_current = cc;
            }
          }
          const UInt i = kk_gauges_max/N_cols;

          saveTemporalSequence ( XX_gauges, time, output_dir + "waterSurfaceHeight_" + std::to_string(number), H_candidate );
          //saveTemporalSequence ( XX_gauges, time, output_dir + "waterSurfaceHeight_" + std::to_string(number), H[ kk_gauges_max ] );
          saveTemporalSequence ( XX_gauges, time, output_dir + "waterSurfaceMassFlux_" + std::to_string(number),
            H[ kk_gauges_max ] * std::sqrt ( std::pow ( ( ( v[ kk_gauges_max ]     + v[ kk_gauges_max + N_cols ] ) / 2. ), 2. ) +
             std::pow ( ( ( u[ kk_gauges_max - i ] + u[ kk_gauges_max - i + 1 ]  ) / 2. ), 2. ) ) );
          saveTemporalSequence ( XX_gauges, time, output_dir + "SolidFlux_" + std::to_string(number),
            h_sd[ kk_gauges_max ] * std::sqrt ( std::pow ( ( ( v[ kk_gauges_max ]     + v[ kk_gauges_max + N_cols ] ) / 2. ), 2. ) +
              std::pow ( ( ( u[ kk_gauges_max - i ] + u[ kk_gauges_max - i + 1 ]  ) / 2. ), 2. ) ) );
        }
      }
      
      
      // sediment production zones,
      for ( const UInt& k : idBasinVect )
        {
          W_Gav_cum[ k ] += W_Gav[ k ] * ( dt_DSV/dt_sed );
        }

      
      if ( spit_out_solutions_each_time_step )
      {
        iter++;

        saveSolution ( output_dir + "u_",     "u", N_rows, N_cols, xllcorner_staggered_u, yllcorner_staggered_u, pixel_size, NODATA_value, iter, u, v, H                            );
        saveSolution ( output_dir + "v_",     "v", N_rows, N_cols, xllcorner_staggered_v, yllcorner_staggered_v, pixel_size, NODATA_value, iter, u, v, H                            );
        saveSolution ( output_dir + "H_",     " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, iter, u, v, H                            );
        saveSolution ( output_dir + "hsd_",   " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, iter, u, v, h_sd                         );
        saveSolution ( output_dir + "w_cum_", " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, iter, u, v, W_Gav_cum                    );
        saveSolution ( output_dir + "ET_",    " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, iter, u, v, ET.ET_vec                    );
        saveSolution ( output_dir + "q_",     " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, iter, u, v, precipitation.DP_cumulative  );
        saveSolution ( output_dir + "p_",     " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, iter, u, v, precipitation.DP_total       );
        saveSolution ( output_dir + "f_",     " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, iter, u, v, precipitation.DP_infiltrated );
        saveSolution ( output_dir + "hG_",    " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, iter, u, v, h_G                          );
        saveSolution ( output_dir + "hsn_",   " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, iter, u, v, h_sn                         );  

      }
      else
      {
        
        if ( std::floor ( time / (frequency_save * 3600) ) > std::floor ( (time - dt_DSV) / (frequency_save * 3600) ) )
        {
          iter++;

          const auto & currentDay = iter; // std::floor ( time / (frequency_save * 3600) )

          std::cout << "Saving solution ..., current day " << iter << std::endl;

          saveSolution ( output_dir + "u_",     "u", N_rows, N_cols, xllcorner_staggered_u, yllcorner_staggered_u, pixel_size, NODATA_value, currentDay, u, v, H                            );
          saveSolution ( output_dir + "v_",     "v", N_rows, N_cols, xllcorner_staggered_v, yllcorner_staggered_v, pixel_size, NODATA_value, currentDay, u, v, H                            );
          saveSolution ( output_dir + "H_",     " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, currentDay, u, v, H                            );
          saveSolution ( output_dir + "hsd_",   " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, currentDay, u, v, h_sd                         );
          saveSolution ( output_dir + "w_cum_", " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, currentDay, u, v, W_Gav_cum                    );
          saveSolution ( output_dir + "ET_",    " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, currentDay, u, v, ET.ET_vec                    );
          saveSolution ( output_dir + "q_",     " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, currentDay, u, v, precipitation.DP_cumulative  );
          saveSolution ( output_dir + "p_",     " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, currentDay, u, v, precipitation.DP_total       );
          saveSolution ( output_dir + "f_",     " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, currentDay, u, v, precipitation.DP_infiltrated );
          saveSolution ( output_dir + "hG_",    " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, currentDay, u, v, h_G                          );
          saveSolution ( output_dir + "hsn_",   " ", N_rows, N_cols, xllcorner,             yllcorner,             pixel_size, NODATA_value, currentDay, u, v, h_sn                         );

        }
      }
      
      // +-----------------------------------------------+
      // |               Update time step                |
      // +-----------------------------------------------+

      dt_DSV = maxdt(u, v, g, maxH, pixel_size);
      dt_DSV = dt_DSV < dt_DSV_given ? dt_DSV : dt_DSV_given;

      compute_dt_adaptive (H, H_old, H_oldold, idBasinVect_excluded, dt_DSV, 1.e-5,
        time, timed, timedd);
      
      c1_DSV_ = c1_DSV (dt_DSV, pixel_size);
      c2_DSV_ = c2_DSV (g, c1_DSV_); 
      c3_DSV_ = c3_DSV (g, c1_DSV_);
      
      
      if ((time+dt_DSV) > t_final)
      {
        dt_DSV = t_final - time;
        check_last = true;
      }
      
      
      toc ("file output");
      
      
      
      
    } // End Time Loop

  }


  MPI_Barrier (MPI_COMM_WORLD); 
  if (rank == 0) { print_timing_report (); }
  MPI_Finalize();

  return ( EXIT_SUCCESS );

}











