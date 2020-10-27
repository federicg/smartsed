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

    @author      Federico     Gatti               MOX Politecnico di Milano   <federico.gatti@polimi.it>
    @mantainer
    
    @supervisors       Luca           Bonaventura  MOX Politecnico di Milano   <luca.bonaventura@polimi.it>
                 Alessandra Menafoglio     MOX Politecnico di Milano   <alessandra.menafoglio@polimi.it>
                 Laura          Longoni          Applied Geology Politecnico di Milano   <laura.longoni@polimi.it>
    
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

#include "utils_H.h"

//! Parse library 
#include "GetPot.hpp" 

//! for parallelize with openmp
#if defined(_OPENMP)
#include <omp.h>
#endif

int main( int argc, char** argv )
{

    //
    #if defined(_OPENMP)
        std::clock_t start;
    #else
        Real start;
    #endif
    
    Real duration;
         
         
         

    // Reading parameters through GetPot
    GetPot command_line( argc, argv );
    const std::string dataFileName = command_line.follow( "SMARTSED_input", 2, "-f", "--file" );
    GetPot dataFile( dataFileName );


    const UInt        currentSimNumber       = command_line.follow( 2, "-sim" );
    const std::string friction_model         = dataFile( "physics/friction_model", "None" );
    const Real        n_manning              = dataFile( "physics/n_manning", 0.01 );
    
    const Real        height_thermometer     = dataFile( "files/meteo_data/height_thermometer", 200. );
    
    
    const UInt        steps_per_hour         = dataFile( "discretization/steps_per_hour", 10 ); // senza contare quello n = 0
    const Real        max_Days               = dataFile( "discretization/max_Days", 20 );
    const Real        H_min                  = dataFile( "discretization/H_min", 0.001 );  // threshold under which there is no mass flux  
    const Real        T_thr                  = dataFile( "discretization/T_thr", 0 );

    const bool        direct_method          = dataFile( "linear_solver/direct_method", true );
     
    const Real        X_gauges               = dataFile( "discretization/X_gauges", 530850.173 );
    const Real        Y_gauges               = dataFile( "discretization/Y_gauges", 5077721.741 );
    const Real        X_1                    = dataFile( "discretization/X_1", 530850.173 );
    const Real        Y_1                    = dataFile( "discretization/Y_1", 5077721.741 );
    const Real        X_2                    = dataFile( "discretization/X_2", 530850.173 );
    const Real        Y_2                    = dataFile( "discretization/Y_2", 5077721.741 );
    const bool        save_temporal_sequence = dataFile( "discretization/save_temporal_sequence", false );
    const bool        isNonReflectingBC      = dataFile( "discretization/isNonReflectingBC", false );

    const Vector2D    XX_gauges( std::array<Real,2>{{ X_gauges, Y_gauges }} );
    const Vector2D    XX_1     ( std::array<Real,2>{{ X_1,      Y_1 }} );
    const Vector2D    XX_2     ( std::array<Real,2>{{ X_2,      Y_2 }} );
    
    // max_Days non è detto che venga raggiunto --> fa da estremo superiore
     
//    const UInt ndata   = std::round( max_Days * 24 / time_spacing ); // if max_Days = 1 e time_spacing = 24 --> ndata = 1
//    const UInt nstep   = steps_per_hour * time_spacing * ndata;
//    const Real t_final = time_spacing * 3600 * ndata; // [seconds]    // initial condition is set at 0 h in the pluviometer time frame
//    const Real dt_DSV      = t_final / nstep; // [seconds]

    const UInt nstep   = steps_per_hour * max_Days * 24;
    const Real t_final = max_Days * 24 * 3600;
    const Real dt_DSV  = t_final / Real( nstep );
    
    const Real scale_coeff = dataFile( "discretization/scale_coeff", 1. );
    


    std::cout << "------------------------ "                            << std::endl;
    std::cout << "friction_model         = " << friction_model          << std::endl;
    std::cout << "n_manning              = " << n_manning               << std::endl;
    std::cout << "steps_per_hour         = " << steps_per_hour          << std::endl;
    std::cout << "max_Days               = " << max_Days                << std::endl;
    std::cout << "nstep                  = " << nstep                   << std::endl;
    std::cout << "t_final                = " << t_final    << " sec."   << std::endl;
    std::cout << "dt_DSV                 = " << dt_DSV     << " sec."   << std::endl;
    std::cout << "H_min                  = " << H_min                   << std::endl;
//    std::cout << "scale_coeff            = " << scale_coeff             << std::endl;
    std::cout << "------------------------ "                            << std::endl;
    

    // +-----------------------------------------------+
    // |                 Reading Input files           |
    // +-----------------------------------------------+


    const std::string file_dir       = "../Inputs/";
    const std::string output_dir     = "../Outputs/" + std::to_string( currentSimNumber ) + "/";

    const std::string orography_file = dataFile( "files/orography_file", "dt_DSVM1.txt" );
    const std::string mask_file      = dataFile( "files/mask_file", "Mask_bin.txt" );


    // precipitation files and temperature
    const std::string temperature_file   = dataFile( "files/meteo_data/temperature_file", "Temperature.txt" );



    const bool restart_H             = dataFile( "files/initial_conditions/restart_H",             false );      
    const bool restart_vel           = dataFile( "files/initial_conditions/restart_vel",           false ); 
    const bool restart_snow          = dataFile( "files/initial_conditions/restart_snow",          false );
    const bool restart_sediment      = dataFile( "files/initial_conditions/restart_sediment",      false );
    const bool restart_gravitational = dataFile( "files/initial_conditions/restart_gravitational", false );
    const bool restart_soilMoisture  = dataFile( "files/initial_conditions/restart_soilMoisture",  false );
   
//    std::cout << restart_H << " " << restart_vel << " " << restart_snow << " " << restart_sediment << " " << restart_gravitational << " " << restart_soilMoisture << std::endl;
    
    {
        const std::string cmd_str = "mkdir -p " + output_dir;  
        char * chararray_cmd = new char[ cmd_str.length() + 1 ];
        const char * cmd_bash = strcpy( chararray_cmd, cmd_str.c_str() ); 

        std::system( cmd_bash );
    }

    const std::string ET_model = dataFile( "files/evapotranspiration/ET_model", "None" ); 
    const Real phi_rad         = M_PI / 180 * dataFile( "files/evapotranspiration/latitude_deg", 45. ); 

    const std::string infiltrationModel = dataFile( "files/infiltration/infiltration_model", "None" );
    
       
    /* Static Variables */   
       
    UInt N_rows,
         N_cols,
         N,
         numberCellsInBasin;

    std::vector<UInt> idStaggeredBoundaryVectSouth,
                      idStaggeredBoundaryVectNorth,
                      idStaggeredBoundaryVectWest,
                      idStaggeredBoundaryVectEast,
                      idStaggeredInternalVectHorizontal,
                      idStaggeredInternalVectVertical,
                      idBasinVect,
                      idBasinVectReIndex;

    std::vector<std::array<Real, 2> > Gamma_vect_x,
                                      Gamma_vect_y;

    
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
                      soilMoistureRetention;
                               
    Eigen::VectorXd   eta,
                      H,
                      H_basin,
                      rhs;
                               
    Real pixel_size, // meter / pixel 
         xllcorner,
         yllcorner,
         xllcorner_staggered_u,
         yllcorner_staggered_u,
         xllcorner_staggered_v,
         yllcorner_staggered_v,
         NODATA_value;          
                               
    /*                   */     
    
    
    
    std::cout << " -- Reading Input Files ... " << std::endl;  
    
    
        
    //
    {
        Raster orographyMat( file_dir + orography_file );
        Raster basin_mask  ( file_dir + mask_file      );
        
        if ( basin_mask.cellsize != orographyMat.cellsize )
        {
            std::cout << mask_file << " cellsize and " << orography_file << " cellsize are not equal" << std::endl;
            exit ( -1 );
        }
        
        pixel_size   = Real( command_line.follow( 2, "-scale" ) ) * basin_mask.cellsize;
        
        std::cout << "cell resolution for the current simulation = " << pixel_size << " meters" << std::endl;
        std::cout << "-------------------- "                                                    << std::endl;
        
        xllcorner    = basin_mask.xllcorner;
        yllcorner    = basin_mask.yllcorner;
        NODATA_value = basin_mask.NODATA_value;
        if ( basin_mask.cellsize <= pixel_size )
        {

            const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
            "dem=raster('" + file_dir + orography_file + "');" +
            "basin=raster('" + file_dir + mask_file + "');" +
            "basin=aggregate(basin," + std::to_string( command_line.follow( 2, "-scale" ) ) + ");" +
            "dem=aggregate(dem," + std::to_string( command_line.follow( 2, "-scale" ) ) + ");" +
            "values(basin)[values(basin)>0]=1;" +
            "writeRaster( dem, file=paste0('" + output_dir + "DEM.asc'), overwrite=TRUE );" +
            "writeRaster( basin, file=paste0('" + output_dir + "basin_mask.asc'), overwrite=TRUE )\"";
            //std::cout << bashCommand << std::endl;
            std::system( bashCommand.c_str() );
        }
        else
        {
            std::cout << "Basin mask greater than simulation resolution i.e. " << pixel_size << std::endl;
            exit( -1. );
        }
        
        if ( dataFile( "discretization/FillSinks", false ) )
        {
            
            std::string bashCommand = std::string( "python3 -c " ) + "\"import os; import sys; import gdal; cwd = os.getcwd();" +
            "sys.path.append( cwd + '/../DeterministicProgram/include/richdem' ); import richdem as rd;" +
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

        Raster orographyMat( output_dir + "DEM.asc" );
        Raster basin_mask  ( output_dir + "basin_mask.asc" );
        
        
        N_rows = basin_mask.Coords.rows();
        N_cols = basin_mask.Coords.cols();
        
        //std::cout << orographyMat.Coords.rows() << " " << orographyMat.Coords.cols() << std::endl;
        
        N      = N_rows * N_cols;
        
        
        H                    .resize( N );
        orography            .resize( N );
        basin_mask_Vec       .resize( N );
        eta                  .resize( N );
        h_G                  .resize( N );
        h_sd                 .resize( N );
        h_sn                 .resize( N );
        S_coeff              .resize( N );
        W_Gav                .resize( N );
        W_Gav_cum            .resize( N );
        Res_x                .resize( N );
        Res_y                .resize( N );
        Z_Gav                .resize( N );
        d_90                 .resize( N );
        soilMoistureRetention.resize( N );
        hydraulic_conductivity.resize( N );
        
        u             .resize( ( N_cols + 1 ) * N_rows );
        v             .resize( ( N_rows + 1 ) * N_cols );
        n_x           .resize( u.size( ) );
        n_y           .resize( v.size( ) );
        u_star        .resize( u.size( ) );
        v_star        .resize( v.size( ) );
        
        Gamma_vect_x         .resize( u.size( ) );
        Gamma_vect_y         .resize( v.size( ) );
        
        h_interface_x .resize( u.size( ) );
        h_interface_y .resize( v.size( ) );
        
        slope_x   .resize( u.size( ) );
        slope_y   .resize( v.size( ) );
        
        

        
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto k = j + i * N_cols;
                basin_mask_Vec[ k ] = basin_mask.Coords.coeff( i,j ) > 0;
                orography[ k ]      = orographyMat.Coords.coeff( i,j );
            }
        }
        
    }
    
    if ( restart_H )
    {
        {
            const std::string H_file = dataFile( "files/initial_conditions/H_file", "H.txt" );
            Raster HMat( file_dir + H_file );
            
            if ( HMat.cellsize <= pixel_size )
            {
                const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
                "basin=raster('" + output_dir + "basin_mask.asc" + "');" +
                "H=raster('" + file_dir + H_file + "');" +
                "H=resample(H,basin,method='bilinear');" +
                "values(H)[is.na(values(H))]=0;" +
                "writeRaster( H, file=paste0('" + output_dir + "H_0.asc'), overwrite=TRUE )\"";
                std::system( bashCommand.c_str() );
                
            }
            else
            {
                std::cout << "Error! resolution of surface water is greater than simulation resolution i.e. " << pixel_size << std::endl;
                exit( -1. );
            }
        }
        Raster HMat( output_dir + "H_0.asc" );
        
        
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto k = j + i * N_cols;
                H( k )   = HMat.Coords.coeff( i,j ) * basin_mask_Vec[ k ];
                eta( k ) = H( k ) + orography[ k ];
            }
        }
        
    }
    else
    {
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto k = j + i * N_cols;
                H( k )   = 0.;
                eta( k ) = orography[ k ];
            }
        }
        
        saveSolution( output_dir + "H_0", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, H );
    }
    
    
    if ( restart_vel )
    {
        {
            const std::string restart_vel_u_file = dataFile( "files/initial_conditions/vel_file_u", "restart_vel_u.txt" );
            Raster vel_u_Mat( file_dir + restart_vel_u_file );
            
            if ( vel_u_Mat.cellsize <= pixel_size )
            {
                const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
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
                //std::cout << bashCommand << std::endl;
                std::system( bashCommand.c_str() );
                
            }
            else
            {
                std::cout << "Error! resolution of horizontal velocity is greater than simulation resolution i.e. " << pixel_size << std::endl;
                exit( -1. );
            }
        }
        Raster vel_u_Mat( output_dir + "u_0.asc" );
        //std::cout << vel_u_Mat.Coords.rows() << " " << vel_u_Mat.Coords.cols() << std::endl;
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j <= N_cols; j++ )
            {
                const auto Id = j + ( N_cols + 1 ) * i; // u
                u[ Id ] = vel_u_Mat.Coords.coeff( i,j );
            }
        }
                             
        xllcorner_staggered_u = vel_u_Mat.xllcorner;
        yllcorner_staggered_u = vel_u_Mat.yllcorner;
        
        
        {
            const std::string restart_vel_v_file = dataFile( "files/initial_conditions/vel_file_v", "restart_vel_v.txt" );
            Raster vel_v_Mat( file_dir + restart_vel_v_file );
            
            if ( vel_v_Mat.cellsize <= pixel_size )
            {
                const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
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
                //std::cout << bashCommand << std::endl;
                std::system( bashCommand.c_str() );
                
            }
            else
            {
                std::cout << "Error! resolution of vertical velocity is greater than simulation resolution i.e. " << pixel_size << std::endl;
                exit( -1. );
            }
        }
        Raster vel_v_Mat( output_dir + "v_0.asc" );
        
        xllcorner_staggered_v = vel_v_Mat.xllcorner;
        yllcorner_staggered_v = vel_v_Mat.yllcorner;
        
        for ( UInt i = 0; i <= N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto Id = j + i * N_cols;
                v[ Id ] = vel_v_Mat.Coords.coeff( i,j );
            }
        }
        
    }
    else
    {
        xllcorner_staggered_u = xllcorner - pixel_size / 2.;
        yllcorner_staggered_u = yllcorner;
        
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j <= N_cols; j++ )
            {
                const auto Id = j + ( N_cols + 1 ) * i; // u
                u[ Id ] = 0.;
            }
        }
        
        saveSolution( output_dir + "u_0", "u", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, u );
        
        
        xllcorner_staggered_v = xllcorner;
        yllcorner_staggered_v = yllcorner - pixel_size / 2.;
        
        for ( UInt i = 0; i <= N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto Id = j + i * N_cols;
                v[ Id ] = 0.;
            }
        }
        
        saveSolution( output_dir + "v_0", "v", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, v );
    }

    
    if ( restart_snow )
    {
        {
            const std::string restart_snow_file = dataFile( "files/initial_conditions/snow_file", "h_sn.asc" );
            Raster snow_Mat( file_dir + restart_snow_file );
            
            if ( snow_Mat.cellsize <= pixel_size )
            {
                const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
                "dem=raster('" + output_dir + "DEM.asc" + "');" +
                "hsn=raster('" + file_dir + restart_snow_file + "');" +
                "hsn=resample(hsn,dem,method='bilinear');" +
                "values(hsn)[is.na(values(hsn))]=0;" +
                "writeRaster( hsn, file=paste0('" + output_dir + "hsn_0.asc'), overwrite=TRUE )\"";
                std::system( bashCommand.c_str() );
                
            }
            else
            {
                std::cout << "Error! resolution of snow file is greater than simulation resolution i.e. " << pixel_size << std::endl;
                exit( -1. );
            }
        }
        Raster snow_Mat( output_dir + "hsn_0.asc" );
        
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto k = j + i * N_cols;
                h_sn[ k ] = snow_Mat.Coords.coeff( i,j ) * basin_mask_Vec[ k ];
            }
        }
        
        
    }
    else
    {
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto k = j + i * N_cols;
                h_sn[ k ] = 0.;
            }
        }
        
        saveSolution( output_dir + "hsn_0", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, h_sn );
        
    }

    if ( restart_sediment )
    {
        {
            const std::string restart_sediment_file = dataFile( "files/initial_conditions/sediment_file", "h_sd.asc" );
            Raster sediment_Mat( file_dir + restart_sediment_file );
            
            if ( sediment_Mat.cellsize <= pixel_size )
            {
                const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
                "dem=raster('" + output_dir + "DEM.asc" + "');" +
                "hsd=raster('" + file_dir + restart_sediment_file + "');" +
                "hsd=resample(hsd,dem,method='bilinear');" +
                "values(hsd)[is.na(values(hsd))]=0;" +
                "writeRaster( hsd, file=paste0('" + output_dir + "hsd_0.asc'), overwrite=TRUE )\"";
                std::system( bashCommand.c_str() );
                
            }
            else
            {
                std::cout << "Error! resolution of sediment file is greater than simulation resolution i.e. " << pixel_size << std::endl;
                exit( -1. );
            }
        }
        Raster sediment_Mat( output_dir + "hsd_0.asc" );
        
        
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto k = j + i * N_cols;
                h_sd[ k ] = sediment_Mat.Coords.coeff( i,j ) * basin_mask_Vec[ k ];
            }
        }
    }
    else
    {
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto k = j + i * N_cols;
                h_sd[ k ] = 0.;
            }
        }
        
        saveSolution( output_dir + "hsd_0", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, h_sd );
        
    }
    
    
    if ( restart_gravitational )
    {
        {
            const std::string restart_gravitational_file = dataFile( "files/initial_conditions/gravitational_file", "h_G.asc" );
            Raster gravitational_Mat( file_dir + restart_gravitational_file );
            
            if ( gravitational_Mat.cellsize <= pixel_size )
            {
                const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
                "dem=raster('" + output_dir + "DEM.asc" + "');" +
                "hG=raster('" + file_dir + restart_gravitational_file + "');" +
                "hG=resample(hG,dem,method='bilinear');" +
                "values(hG)[is.na(values(hG))]=0;" +
                "writeRaster( hG, file=paste0('" + output_dir + "hG_0.asc'), overwrite=TRUE )\"";
                std::system( bashCommand.c_str() );
                
            }
            else
            {
                std::cout << "Error! resolution of gravitational file is greater than simulation resolution i.e. " << pixel_size << std::endl;
                exit( -1. );
            }
        }
        Raster gravitational_Mat( output_dir + "hG_0.asc" );
        
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto k = j + i * N_cols;
                h_G[ k ] = gravitational_Mat.Coords.coeff( i,j ) * basin_mask_Vec[ k ];
            }
        }
        
    }
    else
    {
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto k = j + i * N_cols;
                h_G[ k ] = 0.;
            }
        }
        
        saveSolution( output_dir + "hG_0", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, h_G );
    }
    

    // +-----------------------------------------------+
    // |    Construct soilMoistureRetention vector     |
    // +-----------------------------------------------+
    
    
    if ( !restart_soilMoisture )
    {
       

        std::vector<Int> corineCode_Vec( N ),
                         HSG( N );

        std::vector<Real> clayPercentage_Vec( N ),
                          sandPercentage_Vec( N ),
                          X_Gav( N ),
                          Y_Gav( N );

    

        if ( infiltrationModel != "None" )
        {
            
            const std::string corineCode_file = dataFile( "files/infiltration/corineCode_file", "CLC_RASTER.txt" );
            
            
            // interpolate CLC to make sure to match correct dimensions
            {
                const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
                "dem=raster('" + output_dir + "DEM.asc" + "');" +
                "clc=raster('" + file_dir + corineCode_file + "');" +
                "clc=resample(clc,dem,method='ngb');" +
                "values(clc)[is.na(values(clc))]=0;" +
                "writeRaster( clc, file=paste0('" + output_dir + "CLC.asc'), overwrite=TRUE )\"";
                std::system( bashCommand.c_str() );
                
            }

            Raster corineCode    ( output_dir + "CLC.asc" ),
                   clayPercentage( output_dir + "/clay_sim_" + std::to_string( currentSimNumber ) + ".asc" ),
                   sandPercentage( output_dir + "/sand_sim_" + std::to_string( currentSimNumber ) + ".asc" );
            

            if ( corineCode.cellsize != pixel_size )
            {
                std::cout << "Please check that the " << corineCode_file << " cellsize is consistent with " 
                          << mask_file << " and "     << orography_file  << " ones" << std::endl;  
                exit( -1. );   
            }



            for ( UInt i = 0; i < N_rows; i++ )
            {
                for ( UInt j = 0; j < N_cols; j++ )
                {
                    const auto k = j + i * N_cols;
                    corineCode_Vec[ k ] = corineCode.Coords.coeff( i,j );
                }
            }
            
            
            
            

            if ( clayPercentage.cellsize != sandPercentage.cellsize )
            {
                std::cout << "Please check that the soil texture files have the same resolution" << std::endl;
                exit( -1. );
            }



            if ( clayPercentage.cellsize == pixel_size )
            {
                

                for ( UInt i = 0; i < N_rows; i++ )
                {
                    for ( UInt j = 0; j < N_cols; j++ )
                    {
      
                        
                        const auto k = j + i * N_cols;

                        clayPercentage_Vec[ k ] = clayPercentage.Coords.coeff( i,j );
                        sandPercentage_Vec[ k ] = sandPercentage.Coords.coeff( i,j );

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
                            
                            Vector2D point( std::array<Real,2>{{ clay,sand }} );

                            
                            Vector2D point_A( std::array<Real,2>{{ 0,1    }} );
                            Vector2D point_B( std::array<Real,2>{{ .1, 1  }} );
                            Vector2D point_C( std::array<Real,2>{{ .1, .9 }} );
                            Vector2D point_D( std::array<Real,2>{{ 0, .9  }} );

                            Vector2D point_E( std::array<Real,2>{{ .1, .5 }} );
                            Vector2D point_F( std::array<Real,2>{{ .2, .5 }} );
                            Vector2D point_G( std::array<Real,2>{{ .2, .9 }} );

                            Vector2D point_H( std::array<Real,2>{{ .2, 0  }} );
                            Vector2D point_I( std::array<Real,2>{{ .4, 0  }} );
                            Vector2D point_L( std::array<Real,2>{{ .4, .5 }} );

                            Vector2D point_M( std::array<Real,2>{{ 1, 0   }} );
                            Vector2D point_N( std::array<Real,2>{{ 1, .5  }} );


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




                            std::pair<Real, Int> min = std::make_pair( 1.e4, -1 );

                            for ( auto ii = 0; ii < vv.size(); ii += 4 )
                            {

                                const auto A = vv[ ii ],
                                           D = vv[ ii+1 ],
                                           C = vv[ ii+2 ],
                                           B = vv[ ii+3 ];

                                Real d1 = 1.e4, 
                                     d2 = 1.e4,
                                     d3 = 1.e4,
                                     d4 = 1.e4;


                                Vector2D e1( std::array<Real,2>{{ 1,0 }} ),
                                         e2( std::array<Real,2>{{ 0,1 }} );


                                if ( e1.dot( point - D ) >= 0 && e1.dot( point - D ) <= ( C( 0 ) - D( 0 ) ) ) d1 = std::abs( ( point - D ).dot( e2 ) );
                                if ( e2.dot( point - D ) >= 0 && e2.dot( point - D ) <= ( A( 1 ) - D( 1 ) ) ) d2 = std::abs( ( point - D ).dot( e1 ) );
                                if ( e1.dot( point - A ) >= 0 && e1.dot( point - A ) <= ( B( 0 ) - A( 0 ) ) ) d3 = std::abs( ( point - A ).dot( e2 ) );
                                if ( e2.dot( point - C ) >= 0 && e2.dot( point - C ) <= ( B( 1 ) - C( 1 ) ) ) d4 = std::abs( ( point - C ).dot( e1 ) );

                                if ( d1 < min.first ) min = std::pair<Real, Int>( d1, ii/4 );
                                if ( d2 < min.first ) min = std::pair<Real, Int>( d2, ii/4 );
                                if ( d3 < min.first ) min = std::pair<Real, Int>( d3, ii/4 );
                                if ( d4 < min.first ) min = std::pair<Real, Int>( d4, ii/4 );


                            }

                            const auto & id = min.second;
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
                                exit( -1. );
                            }    





                        }
                        else
                        {
                            //std::cout << clay << " " << sand << std::endl;
                            HSG[ k ] = -1;
                        }          


                    }

                }

            }
            else
            {
                std::cout << "Error! resolution of soil texture files are not equal to the simulation resolution i.e. " << pixel_size << std::endl;
                exit( -1. );
            }

        }
        else
        {
            saveSolution( output_dir + "CLC", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, corineCode_Vec );
        }


 
        const Real S_0 = .254;  // 254 mm

        const auto CN_map = createCN_map( );

        for ( UInt k = 0; k < N; k++ )
        {
            const auto key = std::array<Int,2>{{ corineCode_Vec[ k ], HSG[ k ] }};

            const auto it = CN_map.find( key );


            if ( it != CN_map.end( ) )
            {
                soilMoistureRetention[ k ] = S_0 * ( 100 / Real( it->second ) - 1 ) * basin_mask_Vec[ k ];
//                std::cout << soilMoistureRetention[ k ] << " " << Real( it->second ) << std::endl;
            }
            else
            {
                soilMoistureRetention[ k ] = 0;
            }
        }
//        exit(1);
        
        saveSolution( output_dir + "soilMoistureRetention", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, soilMoistureRetention );
        
        
        // build X_Gav and Y_Gav
        const auto CN_Gav_map = createCN_map_Gav( );
        
        for ( UInt k = 0; k < N; k++ )
        {
            const auto key = corineCode_Vec[ k ];

            const auto it = CN_Gav_map.find( key );
            

            if ( it != CN_Gav_map.end( ) ) // remains zero else
            {
                X_Gav[ k ] = it->second[ 0 ];
                Y_Gav[ k ] = it->second[ 1 ];
            }
            
            Z_Gav[ k ] = X_Gav[ k ] * Y_Gav[ k ];
            
        }
        
        
        
        // build d_10 (for k_c) and d_90 (frictionClass)
        const auto d_10 = compute_d_perc( clayPercentage_Vec, sandPercentage_Vec, 10 );
        
        //std::cout << "maximum and minimum d_10  " << *std::max_element( d_10.begin(), d_10.end() ) << " " << *std::min_element( d_10.begin(), d_10.end() ) << std::endl;
        
        for ( UInt i = 0; i < d_10.size(); i++ )
        {
            // Equations for hydraulic conductivity estimation from particle size distribution: A dimensional analysis
            // Ji-Peng Wang1, Bertrand François, and Pierre Lambert
            const Real C_H = 6.54e-4;
            const Real gravity = 9.81;
            const Real kin_visc = 0.89e-6;
            
            hydraulic_conductivity[ i ] = C_H * gravity / kin_visc * std::pow( d_10[ i ], 2. );
        }
        
        
        d_90 = compute_d_perc( clayPercentage_Vec, sandPercentage_Vec, 90 );
        
        saveSolution( output_dir + "d_10",                   " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, d_10 );
        saveSolution( output_dir + "d_90",                   " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, d_90 );
        saveSolution( output_dir + "k_c",                    " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, hydraulic_conductivity );

        
        //std::cout << "maximum and minimum d  " << *std::max_element( d_90.begin(), d_90.end() ) << " " << *std::min_element( d_90.begin(), d_90.end() ) << std::endl;
        
    }
    else
    {
        {
            const std::string restart_soilMoistureRetention_file = dataFile( "files/initial_conditions/soilMoistureFile", "soilMoisture.asc" ),
                              restart_CLC_file                   = dataFile( "files/initial_conditions/corineCode_file", "CLC.asc" ),
                              restart_clay_file                  = dataFile( "files/initial_conditions/clay_file", "clay.asc" ),
                              restart_sand_file                  = dataFile( "files/initial_conditions/sand_file", "sand.asc" );
            
            Raster soilMoistureRetention_Mat( file_dir + restart_soilMoistureRetention_file ),
                   CLC_Mat                  ( file_dir + restart_CLC_file ),
                   clayPercentage_Mat       ( file_dir + restart_clay_file ),
                   sandPercentage_Mat       ( file_dir + restart_sand_file );
            
            if ( soilMoistureRetention_Mat.cellsize <= pixel_size )
            {
                const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
                "dem=raster('" + output_dir + "DEM.asc" + "');" +
                "soilMoistureRetention=raster('" + file_dir + restart_soilMoistureRetention_file + "');" +
                "soilMoistureRetention=resample(soilMoistureRetention,dem,method='bilinear');" +
                "values(soilMoistureRetention)[is.na(values(soilMoistureRetention))]=0;" +
                "writeRaster( soilMoistureRetention, file=paste0('" + output_dir + "soilMoistureRetention.asc'), overwrite=TRUE )\"";
                std::system( bashCommand.c_str() );
                
            }
            else
            {
                std::cout << "Error! resolution of soilMoistureRetention file is greater than simulation resolution i.e. " << pixel_size << std::endl;
                exit( -1. );
            }
            
            if ( CLC_Mat.cellsize <= pixel_size )
            {
                const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
                "dem=raster('" + output_dir + "DEM.asc" + "');" +
                "CLC=raster('" + file_dir + restart_CLC_file + "');" +
                "CLC=resample(CLC,dem,method='bilinear');" +
                "values(CLC)[is.na(values(CLC))]=0;" +
                "writeRaster( CLC, file=paste0('" + output_dir + "CLC.asc'), overwrite=TRUE )\"";
                std::system( bashCommand.c_str() );
                
            }
            else
            {
                std::cout << "Error! resolution of corineCode_file file is greater than simulation resolution i.e. " << pixel_size << std::endl;
                exit( -1. );
            }
            
            if ( clayPercentage_Mat.cellsize <= pixel_size )
            {
                const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
                "dem=raster('" + output_dir + "DEM.asc" + "');" +
                "clay=raster('" + file_dir + restart_clay_file + "');" +
                "clay=resample(clay,dem,method='bilinear');" +
                "values(clay)[is.na(values(clay))]=0;" +
                "writeRaster( clay, file=paste0('" + output_dir + "clay.asc'), overwrite=TRUE )\"";
                std::system( bashCommand.c_str() );
                
            }
            else
            {
                std::cout << "Error! resolution of clay file is greater than simulation resolution i.e. " << pixel_size << std::endl;
                exit( -1. );
            }
            
            if ( sandPercentage_Mat.cellsize <= pixel_size )
            {
                const std::string bashCommand = std::string( "Rscript -e " ) + "\"library(raster);" +
                "dem=raster('" + output_dir + "DEM.asc" + "');" +
                "sand=raster('" + file_dir + restart_clay_file + "');" +
                "sand=resample(sand,dem,method='bilinear');" +
                "values(sand)[is.na(values(sand))]=0;" +
                "writeRaster( sand, file=paste0('" + output_dir + "sand.asc'), overwrite=TRUE )\"";
                std::system( bashCommand.c_str() );
                
            }
            else
            {
                std::cout << "Error! resolution of sand file is greater than simulation resolution i.e. " << pixel_size << std::endl;
                exit( -1. );
            }
            
        }
        
        Raster soilMoistureRetention_Mat( output_dir + "soilMoistureRetention.asc" ),
               CLC_Mat                  ( output_dir + "CLC.asc" ),
               clay_Mat                 ( output_dir + "clay.asc" ),
               sand_Mat                 ( output_dir + "sand.asc" );
        
        
        std::vector<Int> corineCode_Vec( N );

        std::vector<Real> clayPercentage_Vec( N ),
                          sandPercentage_Vec( N ),
                          X_Gav( N ),
                          Y_Gav( N );
        
        
        const auto CN_Gav_map = createCN_map_Gav( );
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const auto k = j + i * N_cols;
                soilMoistureRetention[ k ] = soilMoistureRetention_Mat.Coords.coeff( i,j ) * basin_mask_Vec[ k ];
                corineCode_Vec       [ k ] = CLC_Mat.Coords.coeff( i,j );
                clayPercentage_Vec[ k ]    = clay_Mat.Coords.coeff( i,j );
                sandPercentage_Vec[ k ]    = sand_Mat.Coords.coeff( i,j );
                
                const auto key = corineCode_Vec[ k ];

                const auto it = CN_Gav_map.find( key );
                
                if ( it != CN_Gav_map.end( ) )
                {
                    X_Gav[ k ] = it->second[ 0 ];
                    Y_Gav[ k ] = it->second[ 1 ];
                }
                
                Z_Gav[ k ] = X_Gav[ k ] * Y_Gav[ k ];
                
            }
        }
        
        

        // build d_10 (for k_c) and d_90 (frictionClass)
        const auto d_10 = compute_d_perc( clayPercentage_Vec, sandPercentage_Vec, 10 );
        
    
        
        for ( UInt i = 0; i < d_10.size(); i++ )
        {
            // Equations for hydraulic conductivity estimation from particle size distribution: A dimensional analysis
            // Ji-Peng Wang1, Bertrand François, and Pierre Lambert
            const Real C_H = 6.54e-4;
            const Real gravity = 9.81;
            const Real kin_visc = 0.89e-6;
            
            hydraulic_conductivity[ i ] = C_H * gravity / kin_visc * std::pow( d_10[ i ], 2. );
        }
        
        d_90 = compute_d_perc( clayPercentage_Vec, sandPercentage_Vec, 90 );
        
        
        saveSolution( output_dir + "d_10",                   " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, d_10 );
        saveSolution( output_dir + "d_90",                   " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, d_90 );
        saveSolution( output_dir + "k_c",                    " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, hydraulic_conductivity );
        
    }
           
        
    std::cout << "maximum and minimum hydraulic_conductivity  " << *std::max_element( hydraulic_conductivity.begin(), hydraulic_conductivity.end() ) << " " << *std::min_element( hydraulic_conductivity.begin(), hydraulic_conductivity.end() ) << std::endl;
        

    // +-----------------------------------------------+
    // |                Compute Slopes                 |
    // +-----------------------------------------------+
    
    

    for ( UInt i = 0; i < N_rows; i++ )
    {
        for ( UInt j = 1; j < N_cols; j++ )
        {
            const auto Id = j + i * ( N_cols + 1 );

            slope_x[ Id ] = ( orography[ Id - i ] - orography[ Id - 1 - i ] ) / pixel_size;
            n_x    [ Id ] = - ( orography[ Id - i ] - orography[ Id - 1 - i ] ) / std::abs( ( orography[ Id - i ] - orography[ Id - 1 - i ] ) );

            if ( std::isnan( n_x[ Id ] ) ) n_x[ Id ] = 0;
        }
    }

    for ( UInt i = 0, j = 0; i < N_rows; i++ )
    {
        const auto Id = j + i * ( N_cols + 1 );

        slope_x[ Id ] = slope_x[ Id + 1 ];
        n_x    [ Id ] = n_x    [ Id + 1 ]; //0;
    }

    for ( UInt i = 0, j = N_cols; i < N_rows; i++ )
    {
        const auto Id = j + i * ( N_cols + 1 );

        slope_x[ Id ] = slope_x[ Id - 1 ];
        n_x    [ Id ] = n_x    [ Id - 1 ]; //0;
    }




    for ( UInt i = 1; i < N_rows; i++ )
    {
        for ( UInt j = 0; j < N_cols; j++ )
        {
            const auto Id = j + i * N_cols;

            slope_y[ Id ] = ( orography[ Id ] - orography[ Id - N_cols ] ) / pixel_size;
            n_y    [ Id ] = - ( orography[ Id ] - orography[ Id - N_cols ] ) / std::abs( ( orography[ Id ] - orography[ Id - N_cols ] ) );
            
            if ( std::isnan( n_y[ Id ] ) ) n_y[ Id ] = 0;
        }
    }

    for ( UInt j = 0, i = 0; j < N_cols; j++ )
    {
        const auto Id = j + i * N_cols;

        slope_y[ Id ] = slope_y[ Id + N_cols ];
        n_y    [ Id ] = n_y    [ Id + N_cols ]; //0;
    }
    
    for ( UInt j = 0, i = N_rows; j < N_cols; j++ )
    {
        const auto Id = j + i * N_cols;

        slope_y[ Id ] = slope_y[ Id - N_cols ];
        n_y    [ Id ] = n_y    [ Id - N_cols ]; //0;
    }
    
    
    saveSolution( output_dir + "slope_x", "u", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, slope_x );
    saveSolution( output_dir + "slope_y", "v", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, slope_y );

    // +-----------------------------------------------+
    // |     Compute boundaries of basin domain        |
    // +-----------------------------------------------+
    
    

    computeAdjacencies( basin_mask_Vec,       
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

    numberCellsInBasin = idBasinVect.size();
    H_basin.resize( numberCellsInBasin );
    rhs    .resize( numberCellsInBasin );
    
    
    
    // maybe if H_basin is inizialized the restart
    if ( restart_H )
    {
        UInt h = 0;
        for ( UInt i = 0; i < N_rows; i++ )
        {
            for ( UInt j = 0; j < N_cols; j++ )
            {
                const UInt k = j + i * N_cols;
                
                if ( basin_mask_Vec[ k ] == 1 )
                {
                    H_basin( h ) = H( k ) * scale_coeff;
                    h++;
                }
            }
        }

    }
    else
    {
        for ( UInt k = 0; k < numberCellsInBasin; k++ )
        {
            H_basin( k ) = 0.;
        }
    }

    // +-----------------------------------------------+
    // |                Gavrilovic Coeff.              |
    // +-----------------------------------------------+


    for ( const auto & k : idBasinVect )
    {
        const UInt i = k / N_cols;
        const auto Sga = std::sqrt( std::pow( .5 * ( slope_x[ k + i ] + slope_x[ k + i + 1 ] ), 2. ) + std::pow( .5 * ( slope_y[ k ] + slope_y[ k + N_cols ] ), 2. ) );
        
        Z_Gav[ k ] *= ( .5 + std::sqrt( Sga ) );
        Z_Gav[ k ] = std::pow( Z_Gav[ k ], 1.5 );
    }


    // +-----------------------------------------------+
    // |                    Gauges i,j                 |
    // +-----------------------------------------------+


    UInt kk_gauges = 0;

    if ( save_temporal_sequence )
    {

        const Vector2D XX_O = std::array<Real,2>{{ xllcorner, yllcorner + N_rows * pixel_size }}; // 
    
        auto XX = ( XX_gauges - XX_O ) / pixel_size; // coordinate in the matrix

        XX( 0 ) =   std::round( XX( 0 ) );
        XX( 1 ) = - std::round( XX( 1 ) );

        if ( XX( 0 ) < 0 || XX( 1 ) < 0 ) 
        {
            std::cout << "The gauges in the input file are not good" << std::endl;
            exit( 1. );
        }

        kk_gauges = XX( 1 ) * N_cols + XX( 0 );  
    
    }
    


    // +-----------------------------------------------+
    // |                    Rain                       |
    // +-----------------------------------------------+
    
    const bool is_precipitation = dataFile( "files/meteo_data/precipitation", true );
    const bool constant_precipitation = dataFile( "files/meteo_data/constant_precipitation", true );

    
    Real dt_rain = 0;
    
    Rain precipitation( infiltrationModel, N, dataFile( "files/infiltration/isInitialLoss", false ), dataFile( "files/infiltration/perc_initialLoss", 0.05 ) );
    
    if ( constant_precipitation || !is_precipitation )
    {
        const std::string precipitation_file = dataFile( "files/meteo_data/rain_file", " " );
        
        const Real time_spacing_rain  = dataFile( "files/meteo_data/time_spacing_rain", 1. );
        
        dt_rain = time_spacing_rain * 3600;
        
        const auto ndata_rain = std::round( max_Days * 24 / time_spacing_rain );
        
        precipitation.constant_precipitation( file_dir + precipitation_file, ndata_rain, is_precipitation, time_spacing_rain );
    }
    else // IDW
    {
        
        const std::string precipitation_file_1 = dataFile( "files/meteo_data/rain_file_1", " " );
        const std::string precipitation_file_2 = dataFile( "files/meteo_data/rain_file_2", " " );
        const std::string precipitation_file_3 = dataFile( "files/meteo_data/rain_file_3", " " );
        const std::string precipitation_file_4 = dataFile( "files/meteo_data/rain_file_4", " " );
        const std::string precipitation_file_5 = dataFile( "files/meteo_data/rain_file_5", " " );
        const std::string precipitation_file_6 = dataFile( "files/meteo_data/rain_file_6", " " );
        const std::string precipitation_file_7 = dataFile( "files/meteo_data/rain_file_7", " " );
        const std::string precipitation_file_8 = dataFile( "files/meteo_data/rain_file_8", " " );
        const std::string precipitation_file_9 = dataFile( "files/meteo_data/rain_file_9", " " );
        const std::vector<std::string> precipitation_file = { file_dir + precipitation_file_1,
            file_dir + precipitation_file_2,
            file_dir + precipitation_file_3,
            file_dir + precipitation_file_4,
            file_dir + precipitation_file_5,
            file_dir + precipitation_file_6,
            file_dir + precipitation_file_7,
            file_dir + precipitation_file_8,
            file_dir + precipitation_file_9
        };
        
        const Real time_spacing_rain_1 = dataFile( "files/meteo_data/time_spacing_rain_1", 1. );
        const Real time_spacing_rain_2 = dataFile( "files/meteo_data/time_spacing_rain_2", 1. );
        const Real time_spacing_rain_3 = dataFile( "files/meteo_data/time_spacing_rain_3", 1. );
        const Real time_spacing_rain_4 = dataFile( "files/meteo_data/time_spacing_rain_4", 1. );
        const Real time_spacing_rain_5 = dataFile( "files/meteo_data/time_spacing_rain_5", 1. );
        const Real time_spacing_rain_6 = dataFile( "files/meteo_data/time_spacing_rain_6", 1. );
        const Real time_spacing_rain_7 = dataFile( "files/meteo_data/time_spacing_rain_7", 1. );
        const Real time_spacing_rain_8 = dataFile( "files/meteo_data/time_spacing_rain_8", 1. );
        const Real time_spacing_rain_9 = dataFile( "files/meteo_data/time_spacing_rain_9", 1. );
        const std::vector<Real> time_spacing_rain = { time_spacing_rain_1,
            time_spacing_rain_2,
            time_spacing_rain_3,
            time_spacing_rain_4,
            time_spacing_rain_5,
            time_spacing_rain_6,
            time_spacing_rain_7,
            time_spacing_rain_8,
            time_spacing_rain_9
        };
        
        dt_rain = *std::min_element( time_spacing_rain.begin(), time_spacing_rain.end() ) * 3600;
        
        const Real X_1 = dataFile( "files/meteo_data/X_1", 1. );
        const Real Y_1 = dataFile( "files/meteo_data/Y_1", 1. );
        const Real X_2 = dataFile( "files/meteo_data/X_2", 1. );
        const Real Y_2 = dataFile( "files/meteo_data/Y_2", 1. );
        const Real X_3 = dataFile( "files/meteo_data/X_3", 1. );
        const Real Y_3 = dataFile( "files/meteo_data/Y_3", 1. );
        const Real X_4 = dataFile( "files/meteo_data/X_4", 1. );
        const Real Y_4 = dataFile( "files/meteo_data/Y_4", 1. );
        const Real X_5 = dataFile( "files/meteo_data/X_5", 1. );
        const Real Y_5 = dataFile( "files/meteo_data/Y_5", 1. );
        const Real X_6 = dataFile( "files/meteo_data/X_6", 1. );
        const Real Y_6 = dataFile( "files/meteo_data/Y_6", 1. );
        const Real X_7 = dataFile( "files/meteo_data/X_7", 1. );
        const Real Y_7 = dataFile( "files/meteo_data/Y_7", 1. );
        const Real X_8 = dataFile( "files/meteo_data/X_8", 1. );
        const Real Y_8 = dataFile( "files/meteo_data/Y_8", 1. );
        const Real X_9 = dataFile( "files/meteo_data/X_9", 1. );
        const Real Y_9 = dataFile( "files/meteo_data/Y_9", 1. );
        const std::vector<Real> X = { X_1,
            X_2,
            X_3,
            X_4,
            X_5,
            X_6,
            X_7,
            X_8,
            X_9
        };
        const std::vector<Real> Y = { Y_1,
            Y_2,
            Y_3,
            Y_4,
            Y_5,
            Y_6,
            Y_7,
            Y_8,
            Y_9
        };
        
        const UInt ndata_rain_1 = std::round( max_Days * 24 / time_spacing_rain_1 );
        const UInt ndata_rain_2 = std::round( max_Days * 24 / time_spacing_rain_2 );
        const UInt ndata_rain_3 = std::round( max_Days * 24 / time_spacing_rain_3 );
        const UInt ndata_rain_4 = std::round( max_Days * 24 / time_spacing_rain_4 );
        const UInt ndata_rain_5 = std::round( max_Days * 24 / time_spacing_rain_5 );
        const UInt ndata_rain_6 = std::round( max_Days * 24 / time_spacing_rain_6 );
        const UInt ndata_rain_7 = std::round( max_Days * 24 / time_spacing_rain_7 );
        const UInt ndata_rain_8 = std::round( max_Days * 24 / time_spacing_rain_8 );
        const UInt ndata_rain_9 = std::round( max_Days * 24 / time_spacing_rain_9 );
        const std::vector<UInt> ndata_rain = { ndata_rain_1,
            ndata_rain_2,
            ndata_rain_3,
            ndata_rain_4,
            ndata_rain_5,
            ndata_rain_6,
            ndata_rain_7,
            ndata_rain_8,
            ndata_rain_9
        };
        
        
        
        precipitation.IDW_precipitation( precipitation_file, ndata_rain, time_spacing_rain, X, Y, xllcorner, yllcorner, pixel_size, N_rows, N_cols, idBasinVect );
    }


    // +-----------------------------------------------+
    // |                 Temperature                   |
    // +-----------------------------------------------+

    const Real        time_spacing_temp = dataFile( "files/meteo_data/time_spacing_temp", 1. );
    const std::string format_temp       = dataFile( "files/meteo_data/format_temp", "arpa" );
    const Real        dt_temp           = time_spacing_temp * 3600;
    
    Temperature temp( file_dir + temperature_file,
                      N,
                      max_Days,
                      T_thr,
                      orography,
                      std::round( max_Days * 24 / time_spacing_temp ),
                      steps_per_hour,
                      time_spacing_temp,
                      height_thermometer,
                      format_temp );


    // +-----------------------------------------------+
    // |              Evapotranspiration               |
    // +-----------------------------------------------+

    evapoTranspiration ET( ET_model, N, orography, temp.J, max_Days, phi_rad, height_thermometer );
    

    // +-----------------------------------------------+
    // |                   Core Part                   |
    // +-----------------------------------------------+


            
    upwind H_interface( u, v, N_rows, N_cols );

    frictionClass alfa( friction_model, n_manning, dt_DSV, d_90, dataFile( "files/infiltration/roughness_scale_factor", 100. ), H_min, N_rows, N_cols, slope_x, slope_y );
    
    
    
    const Real dt_min = std::min( dt_rain, dt_temp );
    for ( const auto & kk : hydraulic_conductivity )
    {
        if ( kk * ( dt_min / pixel_size ) > 1. )
        {
            std::cout << "Error! Courant number for gravitational layer is greater than 1" << std::endl;
            exit( -1. );
        }
    }
    
    const Real c1_min = dt_min / pixel_size;
    
    
    
    
    const Real g = 9.81;
    const Real c1_DSV = dt_DSV / pixel_size;
    const Real c2_DSV = g * c1_DSV;
    const Real c3_DSV = g * std::pow( c1_DSV, 2 );

    const Real area = std::pow( pixel_size, 2 ) * 1.e-6; // km^2
    
    std::cout << "c1_DSV " << c1_DSV << "   c2_DSV " << c2_DSV << "   c3_DSV " << c3_DSV << std::endl;
    
    
    const Real slope_y_max = std::max( *std::max_element( slope_y.begin(), slope_y.end() ), std::abs( *std::min_element( slope_y.begin(), slope_y.end() ) ) );
    const Real slope_x_max = std::max( *std::max_element( slope_x.begin(), slope_x.end() ), std::abs( *std::min_element( slope_x.begin(), slope_x.end() ) ) );

    Real dt_sed,
         c1_sed;

    UInt numberOfSteps = 1;
    
    bool isHNegative =  false;

    #if defined(_OPENMP)
       start = omp_get_wtime();
    #else
       start = std::clock();
    #endif
    
    
    
    for ( Int n = 1; n <= nstep; n++ )
    {

        std::cout << "Time step number: " << n << "/" << nstep << " Simulation progress: " << n/Real(nstep)*100 << " %" << " max surface runoff vel. based Courant: " << maxCourant( u, v, c1_DSV ) << " max surface runoff cel. based Courant: " << maxCourant( H, g, c1_DSV ) << " H has been negative: " << isHNegative << std::endl;

        
        
        // Compute interface fluxes via upwind method
        H_interface.computeHorizontal( H,
                                       u,
                                       idStaggeredInternalVectHorizontal,
                                       idStaggeredBoundaryVectWest,
                                       idStaggeredBoundaryVectEast );

        H_interface.computeVertical  ( H,
                                       v,
                                       idStaggeredInternalVectVertical,
                                       idStaggeredBoundaryVectNorth,
                                       idStaggeredBoundaryVectSouth );
        
        

        // Compute alfa coefficients
        alfa.f_x( H_interface.horizontal,
                  u,
                  idStaggeredInternalVectHorizontal,
                  idStaggeredBoundaryVectWest,
                  idStaggeredBoundaryVectEast );
        
        alfa.f_y( H_interface.vertical,
                  v,
                  idStaggeredInternalVectVertical,
                  idStaggeredBoundaryVectNorth,
                  idStaggeredBoundaryVectSouth );
        
    
        
        // update only if necessary  --> governed by temperature dynamics, i.e. time_spacing_temp
        if ( std::floor( ( n - 1 ) * ( dt_DSV / dt_temp ) ) > std::floor( ( n - 2 ) * ( dt_DSV / dt_temp ) ) )
        {
            // Compute temperature map
            temp.computeTemperature( n, steps_per_hour, time_spacing_temp, orography, idBasinVect );
        }
        
        // ET varies daily
        if ( std::floor( ( n - 1 ) * ( dt_DSV / (24*3600) ) ) > std::floor( ( n - 2 ) * ( dt_DSV / (24*3600) ) ) )
        {
            // Get ET rate at the current time
            ET.ET( temp.T_dailyMean, temp.T_dailyMin, temp.T_dailyMax, n, idBasinVect, orography, steps_per_hour );
        }

        
        // update only if necessary
        if ( std::floor( ( n - 1 ) * ( dt_DSV / dt_min ) ) > std::floor( ( n - 2 ) * ( dt_DSV / dt_min ) ) )
        {
            
            // +-----------------------------------------------+
            // |            Gravitational Layer                |
            // +-----------------------------------------------+
            
            // h_G
            // vertical and horizontal residuals for Gravitational Layer 
            computeResiduals( n_x,
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
        }
        
        
        //
        precipitation.computePrecipitation( n,
                                            steps_per_hour,
                                            soilMoistureRetention,
                                            temp.melt_mask,
                                            h_G,
                                            H,
                                            N_rows,
                                            N_cols,
                                            idBasinVect );
        
        // update only if necessary
        if ( std::floor( ( n - 1 ) * ( dt_DSV / dt_min ) ) > std::floor( ( n - 2 ) * ( dt_DSV / dt_min ) ) )
        {
            
            for ( const UInt & k : idBasinVect )
            {
                S_coeff[ k ] = 4.62e-10 * h_sn[ k ] * ( temp.T_raster[ k ] - T_thr ) * temp.melt_mask[ k ];
            }
            
            for ( const auto & k : idBasinVect )
            {
                h_G[ k ] += ( S_coeff[ k ] - ET.ET_vec[ k ] ) * dt_min + precipitation.DP_infiltrated[ k ] * dt_min - c1_min * ( Res_x[ k ] + Res_y[ k ] );
                h_G[ k ] *= ( h_G[ k ] >= 0 ); // to account for evapotranspiration
            }
            
            // +-----------------------------------------------+
            // |              Snow Accumulation                |
            // +-----------------------------------------------+

            
            for ( const auto & k : idBasinVect )
            {
                const auto snow_acc = precipitation.DP_total[ k ] * ( 1 - temp.melt_mask[ k ] ) * dt_min - S_coeff[ k ] * dt_min;
                h_sn[ k ] += snow_acc;
            }
            
        }
        


        // +-----------------------------------------------+
        // |               De Saint Venant                 |
        // +-----------------------------------------------+
        
        
        
        // fill u_star and v_star with a Bilinear Interpolation
        bilinearInterpolation( u,
                               v,
                               u_star,
                               v_star,
                               N_rows,
                               N_cols,
                               dt_DSV,
                               pixel_size,
                               idStaggeredInternalVectHorizontal,
                               idStaggeredInternalVectVertical,
                               idStaggeredBoundaryVectWest,
                               idStaggeredBoundaryVectEast,
                               idStaggeredBoundaryVectNorth,
                               idStaggeredBoundaryVectSouth );
         
         /*
         bilinearInterpolation( u,
                                v,
                                u_star,
                                v_star,
                                N_rows,
                                N_cols,
                                dt_DSV,
                                pixel_size );*/
        
        
         
       
        // row, column and value in the Triplet
        std::vector<Eigen::Triplet<Real> > coefficients;           
        
        
        
        
        buildMatrix( H_interface.horizontal,
                     H_interface.vertical,
                     orography,
                     u_star,
                     v_star,
                     H,
                     N_cols,
                     N_rows,
                     c1_DSV,
                     c3_DSV,
                     H_min,  // 0
                     precipitation.DP_cumulative,
                     dt_DSV,
                     alfa.alfa_x,
                     alfa.alfa_y,
                     idStaggeredInternalVectHorizontal,
                     idStaggeredInternalVectVertical,
                     idStaggeredBoundaryVectWest,
                     idStaggeredBoundaryVectEast,
                     idStaggeredBoundaryVectNorth,
                     idStaggeredBoundaryVectSouth,
                     idBasinVect,
                     idBasinVectReIndex,
                     isNonReflectingBC,
                     scale_coeff,
                     coefficients,
                     rhs );
        
         
        
        SpMat A( numberCellsInBasin, numberCellsInBasin );
        A.setFromTriplets( coefficients.begin(), coefficients.end() );         
        A.makeCompressed();
        
    
        
        
        
         
        if ( direct_method )  // Direct Sparse method: Cholesky being A spd 
        {
            Eigen::SimplicialLDLT<SpMat, Eigen::Upper> solver;   // se metti Eigen::Upper|Eigen::Lower non funziona per niente
            
   
            solver.compute( A );
   
            if ( solver.info() != Eigen::Success ) 
            {
                std::cout << "Decomposition Failed" << std::endl;
                exit( -1. );
            }      
        
            H_basin = solver.solve( rhs );
            
        }
        else
        {
                        
            //std::cout << Eigen::nbThreads( ) << std::endl;

            Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<Real>> cg;
            cg.compute( A );
            

            H_basin = cg.solveWithGuess( rhs, H_basin );

            std::cout << "# iterations:    " << cg.iterations() << std::endl;
            std::cout << "estimated error: " << cg.error()      << std::endl;

        }
        
//        std::cout << "save matrix" << std::endl;
//        saveMatrix(A, output_dir+"A_"+std::to_string(n)+".txt");
//        saveVector(rhs, output_dir+"b_"+std::to_string(n)+".txt");
//        saveVector(H_basin, output_dir+"H_"+std::to_string(n)+".txt");
        
        for ( const UInt & Id : idBasinVect )
        {
            const UInt IDreIndex = idBasinVectReIndex[ Id ];
            H( Id ) = H_basin( IDreIndex ) / scale_coeff;
            eta( Id ) = H( Id ) + orography[ Id ];
        }
        
        
        // +-----------------------------------------------+
        // |                  Update u, v                  |
        // +-----------------------------------------------+
        

        
        updateVel( u,
                   v,
                   u_star,
                   v_star,
                   H_interface.horizontal,
                   H_interface.vertical,
                   alfa.alfa_x,
                   alfa.alfa_y,
                   N_rows,
                   N_cols,
                   c2_DSV,
                   H_min,
                   eta,
                   orography,
                   idStaggeredInternalVectHorizontal,
                   idStaggeredInternalVectVertical,
                   idStaggeredBoundaryVectWest,
                   idStaggeredBoundaryVectEast,
                   idStaggeredBoundaryVectNorth,
                   idStaggeredBoundaryVectSouth,
                   isNonReflectingBC );
        
        
        std::cout << "min H: " << H_basin.minCoeff() /scale_coeff << " max H: " << H_basin.maxCoeff() /scale_coeff << std::endl;
        
        for ( const UInt & k : idBasinVect )
        {
            if ( H( k ) < 0 )
            {
                H( k ) *= 0.;
                isHNegative = true;
//                std::cout << "negative H!" << std::endl;
            }
            eta( k ) = H( k ) + orography[ k ];
        }

        

        // +-----------------------------------------------+
        // |             Sediment Transport                |
        // +-----------------------------------------------+


        dt_sed = compute_dt_sediment( 2.5, 1.6, slope_x_max, slope_y_max, u, v, pixel_size, dt_DSV, numberOfSteps );
        c1_sed = dt_sed / pixel_size;

        
        
        // vertical and horizontal residuals truncated for Sediment Transport
        computeResidualsTruncated( u,
                                   v,
                                   N_cols,
                                   N_rows,
                                   N,
                                   c1_sed,
                                   slope_x,
                                   slope_y,
                                   2.5,   // alfa
                                   1.6,   // beta
                                   1, // gamma
                                   idStaggeredInternalVectHorizontal,
                                   idStaggeredInternalVectVertical,
                                   idStaggeredBoundaryVectWest,
                                   idStaggeredBoundaryVectEast,
                                   idStaggeredBoundaryVectNorth,
                                   idStaggeredBoundaryVectSouth,
                                   Gamma_vect_x,
                                   Gamma_vect_y
                                 );
        
        
        
        

        for ( const UInt & k : idBasinVect )
        {
            W_Gav[ k ] = M_PI * Z_Gav[ k ] * std::sqrt( std::abs( ( .1 + .1 * temp.T_raster[ k ] ) * temp.melt_mask[ k ] ) ) * precipitation.DP_total[ k ] * dt_sed;
        }
        
        
        
        for ( UInt kk = 0; kk < numberOfSteps; kk++ )
        {
        
            
            for ( const UInt & Id : idStaggeredInternalVectHorizontal )
            {
                const UInt i       = Id / ( N_cols + 1 ),
                           IDeast  = Id - i,                  // H
                           IDwest  = Id - i - 1;              // H
                
                const Real & h_left  = h_sd[ IDwest ],
                           & h_right = h_sd[ IDeast ];
                       
                
                h_interface_x[ Id ] = Gamma_vect_x[ Id ][ 0 ] * h_right +
                                      Gamma_vect_x[ Id ][ 1 ] * h_left;
                
            }


            for ( const UInt & Id : idStaggeredBoundaryVectWest )
            {
            
                const UInt i = Id / ( N_cols + 1 );

                const Real h_left  = 0, //0,
                           h_right = h_sd[ Id - i ];
                                               

                h_interface_x[ Id ] = Gamma_vect_x[ Id ][ 0 ] * h_right +
                                      Gamma_vect_x[ Id ][ 1 ] * h_left;

               
            }

           
            for ( const UInt & Id : idStaggeredBoundaryVectEast )
            {
            
                const UInt i = Id / ( N_cols + 1 );


                const Real h_left  = h_sd[ Id - i - 1 ],
                           h_right = 0; //0;
                                               

                h_interface_x[ Id ] = Gamma_vect_x[ Id ][ 0 ] * h_right +
                                      Gamma_vect_x[ Id ][ 1 ] * h_left;
        
            }

            for ( const UInt & Id : idBasinVect )
            {
                const UInt i = Id / N_cols;
                
                Res_x[ Id ] = h_interface_x[ Id + 1 + i ] - h_interface_x[ Id + i ];
            }
            
             
            
            
            

            
            for ( const UInt & Id : idStaggeredInternalVectVertical )
            {
                const UInt IDsouth = Id,              // H
                           IDnorth = Id - N_cols;     // H
                
                const Real & h_left  = h_sd[ IDnorth ],
                           & h_right = h_sd[ IDsouth ];
                       
                
                h_interface_y[ Id ] = Gamma_vect_y[ Id ][ 0 ] * h_right +
                                      Gamma_vect_y[ Id ][ 1 ] * h_left;
                
            }


            for ( const UInt & Id : idStaggeredBoundaryVectNorth )
            {

                const Real h_left  = 0, // 0
                           h_right = h_sd[ Id ];
                                                                                                            

                h_interface_y[ Id ] = Gamma_vect_y[ Id ][ 0 ] * h_right +
                                      Gamma_vect_y[ Id ][ 1 ] * h_left;
                
            }


            
            for ( const UInt & Id : idStaggeredBoundaryVectSouth )
            {

                const Real h_left  = h_sd[ Id - N_cols ],
                           h_right = 0; //0;
                                                   

                h_interface_y[ Id ] = Gamma_vect_y[ Id ][ 0 ] * h_right +
                                      Gamma_vect_y[ Id ][ 1 ] * h_left;
                
            }
            
            for ( const UInt & Id : idBasinVect )
            {
                Res_y[ Id ] = h_interface_y[ Id + N_cols ] - h_interface_y[ Id ]; // controllare questo anche in h_G
            }
            

            
            for ( const UInt & Id : idBasinVect )
            {
                h_sd[ Id ] += - ( Res_x[ Id ] + Res_y[ Id ] ) + W_Gav[ Id ];
            }
            
            
            

        }
        
        
        
    
        
        

    
        if ( save_temporal_sequence )
        {
            const UInt i = kk_gauges / N_cols;
            
            saveTemporalSequence( XX_gauges, dt_DSV, n, output_dir + "waterSurfaceHeight", H[ kk_gauges ] );
            saveTemporalSequence( XX_gauges, dt_DSV, n, output_dir + "waterSurfaceMassFlux",
                                 H[ kk_gauges ] * std::sqrt( std::pow( ( ( v[ kk_gauges ]     + v[ kk_gauges + N_cols ] ) / 2. ), 2. ) +
                                                             std::pow( ( ( u[ kk_gauges - i ] + u[ kk_gauges - i + 1 ]  ) / 2. ), 2. ) ) );
            
            
        }
        
        
        
        // +-----------------------------------------------+
        // |             Save The Raster Solution          |
        // +-----------------------------------------------+            
        
        // zone sorgenti
        for ( const UInt & k : idBasinVect )
        {
            W_Gav_cum[ k ] += W_Gav[ k ] * ( dt_DSV / dt_sed );
        }
        
        
        
//        if ( std::floor( ( n - 1 ) * ( dt_DSV / dt_min ) ) > std::floor( ( n - 2 ) * ( dt_DSV / dt_min ) ) )
//        {
//            const auto currentTime = std::floor( n * ( dt_DSV / dt_min ) );
//
//            saveSolution( output_dir + "q_",  " ",  N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentTime, u, v, precipitation.DP_cumulative );
//            saveSolution( output_dir + "p_",  " ",  N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentTime, u, v, precipitation.DP_total );
//            saveSolution( output_dir + "f_",  " ",  N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentTime, u, v, precipitation.DP_infiltrated );
//        }
        
//        std::cout << int(36000 * ( dt_DSV / (24*3600) )) << " " << std::floor( n * ( dt_DSV / (24*3600) ) ) << " " << std::floor( ( n - 1 ) * ( dt_DSV / (24*3600) ) ) << std::endl;
        
        if ( std::floor( ( n + 1 ) * ( dt_DSV / (24*3600) ) ) > std::floor( n * ( dt_DSV / (24*3600) ) ) )
        {
            const auto currentDay = std::floor( ( n + 1 ) * ( dt_DSV / (24*3600) ) );
            
            saveSolution( output_dir + "u_",   "u", N_rows, N_cols, xllcorner_staggered_u, yllcorner_staggered_u, pixel_size, NODATA_value, currentDay, u, v, H );
            saveSolution( output_dir + "v_",   "v", N_rows, N_cols, xllcorner_staggered_v, yllcorner_staggered_v, pixel_size, NODATA_value, currentDay, u, v, H );
            saveSolution( output_dir + "H_",   " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentDay, u, v, H );
            saveSolution( output_dir + "hsd_", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentDay, u, v, h_sd );
            saveSolution( output_dir + "w_cum_",   " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentDay, u, v, W_Gav_cum );
            saveSolution( output_dir + "ET_",  " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentDay, u, v, ET.ET_vec );
            saveSolution( output_dir + "q_",  " ",  N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentDay, u, v, precipitation.DP_cumulative );
            saveSolution( output_dir + "p_",  " ",  N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentDay, u, v, precipitation.DP_total );
            saveSolution( output_dir + "f_",  " ",  N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentDay, u, v, precipitation.DP_infiltrated );
            saveSolution( output_dir + "hG_",  " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentDay, u, v, h_G );
            saveSolution( output_dir + "hsn_", " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, currentDay, u, v, h_sn );
        }
        
//        saveSolution( output_dir + "H_",   " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, n, u, v, H );
//        saveSolution( output_dir + "u_",   "u", N_rows, N_cols, xllcorner_staggered_u, yllcorner_staggered_u, pixel_size, NODATA_value, n, u, v, H );
//        saveSolution( output_dir + "v_",   "v", N_rows, N_cols, xllcorner_staggered_v, yllcorner_staggered_v, pixel_size, NODATA_value, n, u, v, H );
        
//        saveSolution( output_dir + "hG_",  " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, n, u, v, h_G );
//        saveSolution( output_dir + "ET_",  " ", N_rows, N_cols, xllcorner, yllcorner, pixel_size, NODATA_value, n, u, v, ET.ET_vec );
        
//        if ( n == 2 )
//            exit(1);
        
        
         
         
    } // End Time Loop


     
    #if defined(_OPENMP)
        duration = ( omp_get_wtime() - start );
    #else
        duration = ( std::clock() - start ) / ( Real ) CLOCKS_PER_SEC;
    #endif
    
    std::cout << "Operation took " << duration << " seconds" << std::endl;

    return ( EXIT_SUCCESS );
    
}











