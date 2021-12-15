#include "typedefs_H.h"
#include "utils_H.h"

//! std library
#include <cstdlib>

#include <cmath>
#include <numeric>
#include <iosfwd>
#include <string>
#include <limits>

#include <chrono>
#include <ctime>
#include <functional>

//! Eigen library
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"

//! Python interface --> correctly the QGIS path in Makefile
//#include "Python.h"


Vector2D
operator* (Vector2D const& vector, Real const& factor)
{
  Vector2D tmp( vector );
  return tmp *= factor;
}

Vector2D
operator* (Real const& factor, Vector2D const& vector)
{
  Vector2D tmp( vector );
  return tmp *= factor;
}


std::map<Int, std::array<Real,2>>
createCN_map_Gav ()
{
  std::map<Int, std::array<Real,2> > CN;

  // Marrazzo evento Ottobre,
  CN[ 111 ] = std::array<Real,2>{{ 0.05,0.2 }};
  CN[ 112 ] = std::array<Real,2>{{ 0.2,0.2 }};
  CN[ 121 ] = std::array<Real,2>{{ 0.05,0.2 }};
  CN[ 122 ] = std::array<Real,2>{{ 0.2,0.2 }};
  CN[ 123 ] = std::array<Real,2>{{ 0.2,0.2 }};
  CN[ 124 ] = std::array<Real,2>{{ 0.2,0.2 }};
  CN[ 131 ] = std::array<Real,2>{{ 1,1.8 }};
  CN[ 132 ] = std::array<Real,2>{{ 0.8,1.8 }};
  CN[ 133 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 141 ] = std::array<Real,2>{{ 0.5,1.5 }};
  CN[ 142 ] = std::array<Real,2>{{ 0.5,1.5 }};
  CN[ 211 ] = std::array<Real,2>{{ 0.6,1.8 }};
  CN[ 212 ] = std::array<Real,2>{{ 0.6,1.8 }};
  CN[ 213 ] = std::array<Real,2>{{ 0.6,1.8 }};
  CN[ 221 ] = std::array<Real,2>{{ 0.6,2 }};
  CN[ 222 ] = std::array<Real,2>{{ 0.6,1.8 }};
  CN[ 223 ] = std::array<Real,2>{{ 0.6,1.8 }};
  CN[ 231 ] = std::array<Real,2>{{ 0.4,1.6 }};
  CN[ 241 ] = std::array<Real,2>{{ 0.6,1.8 }};
  CN[ 242 ] = std::array<Real,2>{{ 0.6,1.8 }};
  CN[ 243 ] = std::array<Real,2>{{ 0.4,1.8 }};
  CN[ 244 ] = std::array<Real,2>{{ 0.4,1.6 }};
  CN[ 311 ] = std::array<Real,2>{{ 0.05,1.6 }};
  CN[ 312 ] = std::array<Real,2>{{ 0.05,1.6 }};
  CN[ 313 ] = std::array<Real,2>{{ 0.05,1.6 }};
  CN[ 321 ] = std::array<Real,2>{{ 0.6,1.6 }};
  CN[ 322 ] = std::array<Real,2>{{ 0.2,1.6 }};
  CN[ 323 ] = std::array<Real,2>{{ 0.2,1.6 }};
  CN[ 324 ] = std::array<Real,2>{{ 0.4,1.6 }};
  CN[ 331 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 332 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 333 ] = std::array<Real,2>{{ 0.8,1.5 }};
  CN[ 334 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 335 ] = std::array<Real,2>{{ 1,0.2 }};
  CN[ 411 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 412 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 421 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 422 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 423 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 511 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 512 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 521 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 522 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 523 ] = std::array<Real,2>{{ 1,0 }};


  /*
  // Graziotto Luglio-Ottobre,
  CN[ 111 ] = std::array<Real,2>{{ 0.2,0.6 }};
  CN[ 112 ] = std::array<Real,2>{{ 0.6,0.6 }};
  CN[ 121 ] = std::array<Real,2>{{ 0.4,0.6 }};
  CN[ 122 ] = std::array<Real,2>{{ 1,0.6 }};
  CN[ 123 ] = std::array<Real,2>{{ 1,0.6 }};
  CN[ 124 ] = std::array<Real,2>{{ 1,0.6 }};
  CN[ 131 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 132 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 133 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 141 ] = std::array<Real,2>{{ 0.6,2 }};
  CN[ 142 ] = std::array<Real,2>{{ 0.6,2 }};
  CN[ 211 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 212 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 213 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 221 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 222 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 223 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 231 ] = std::array<Real,2>{{ 1,1.8 }};
  CN[ 241 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 242 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 243 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 244 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 311 ] = std::array<Real,2>{{ 1,1.8 }};
  CN[ 312 ] = std::array<Real,2>{{ 0.1,1.8 }};
  CN[ 313 ] = std::array<Real,2>{{ 0.1,1.8 }};
  CN[ 321 ] = std::array<Real,2>{{ 0.6,1.8 }};
  CN[ 322 ] = std::array<Real,2>{{ 0.6,1.8 }};
  CN[ 323 ] = std::array<Real,2>{{ 0.6,1.8 }};
  CN[ 324 ] = std::array<Real,2>{{ 0.6,1.8 }};
  CN[ 331 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 332 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 333 ] = std::array<Real,2>{{ 1,1.8 }};
  CN[ 334 ] = std::array<Real,2>{{ 1,2 }};
  CN[ 335 ] = std::array<Real,2>{{ 1,0.2 }};
  CN[ 411 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 412 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 421 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 422 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 423 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 511 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 512 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 521 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 522 ] = std::array<Real,2>{{ 1,0 }};
  CN[ 523 ] = std::array<Real,2>{{ 1,0 }};*/

  return CN;

}


std::map<std::array<Int,2>, Int>
createCN_map ()
{

  std::map<std::array<Int,2>, Int > CN;

  CN[ std::array<Int,2>{{ 111,0 }} ] = 89;
  CN[ std::array<Int,2>{{ 111,1 }} ] = 92;
  CN[ std::array<Int,2>{{ 111,2 }} ] = 94;
  CN[ std::array<Int,2>{{ 111,3 }} ] = 95;

  CN[ std::array<Int,2>{{ 112,0 }} ] = 77;
  CN[ std::array<Int,2>{{ 112,1 }} ] = 85;
  CN[ std::array<Int,2>{{ 112,2 }} ] = 90;
  CN[ std::array<Int,2>{{ 112,3 }} ] = 92;

  CN[ std::array<Int,2>{{ 121,0 }} ] = 81;
  CN[ std::array<Int,2>{{ 121,1 }} ] = 88;
  CN[ std::array<Int,2>{{ 121,2 }} ] = 91;
  CN[ std::array<Int,2>{{ 121,3 }} ] = 93;

  CN[ std::array<Int,2>{{ 122,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 122,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 122,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 122,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 123,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 123,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 123,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 123,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 124,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 124,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 124,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 124,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 131,0 }} ] = 76;
  CN[ std::array<Int,2>{{ 131,1 }} ] = 85;
  CN[ std::array<Int,2>{{ 131,2 }} ] = 89;
  CN[ std::array<Int,2>{{ 131,3 }} ] = 91;

  CN[ std::array<Int,2>{{ 132,0 }} ] = 81;
  CN[ std::array<Int,2>{{ 132,1 }} ] = 88;
  CN[ std::array<Int,2>{{ 132,2 }} ] = 91;
  CN[ std::array<Int,2>{{ 132,3 }} ] = 93;

  CN[ std::array<Int,2>{{ 133,0 }} ] = 77;
  CN[ std::array<Int,2>{{ 133,1 }} ] = 86;
  CN[ std::array<Int,2>{{ 133,2 }} ] = 91;
  CN[ std::array<Int,2>{{ 133,3 }} ] = 94;

  CN[ std::array<Int,2>{{ 141,0 }} ] = 49;
  CN[ std::array<Int,2>{{ 141,1 }} ] = 69;
  CN[ std::array<Int,2>{{ 141,2 }} ] = 79;
  CN[ std::array<Int,2>{{ 141,3 }} ] = 84;

  CN[ std::array<Int,2>{{ 142,0 }} ] = 68;
  CN[ std::array<Int,2>{{ 142,1 }} ] = 79;
  CN[ std::array<Int,2>{{ 142,2 }} ] = 86;
  CN[ std::array<Int,2>{{ 142,3 }} ] = 89;

  CN[ std::array<Int,2>{{ 211,0 }} ] = 61;
  CN[ std::array<Int,2>{{ 211,1 }} ] = 73;
  CN[ std::array<Int,2>{{ 211,2 }} ] = 81;
  CN[ std::array<Int,2>{{ 211,3 }} ] = 84;

  CN[ std::array<Int,2>{{ 212,0 }} ] = 67;
  CN[ std::array<Int,2>{{ 212,1 }} ] = 78;
  CN[ std::array<Int,2>{{ 212,2 }} ] = 85;
  CN[ std::array<Int,2>{{ 212,3 }} ] = 89;

  CN[ std::array<Int,2>{{ 213,0 }} ] = 62;
  CN[ std::array<Int,2>{{ 213,1 }} ] = 71;
  CN[ std::array<Int,2>{{ 213,2 }} ] = 78;
  CN[ std::array<Int,2>{{ 213,3 }} ] = 81;

  CN[ std::array<Int,2>{{ 221,0 }} ] = 76;
  CN[ std::array<Int,2>{{ 221,1 }} ] = 85;
  CN[ std::array<Int,2>{{ 221,2 }} ] = 90;
  CN[ std::array<Int,2>{{ 221,3 }} ] = 93;

  CN[ std::array<Int,2>{{ 222,0 }} ] = 43;
  CN[ std::array<Int,2>{{ 222,1 }} ] = 65;
  CN[ std::array<Int,2>{{ 222,2 }} ] = 76;
  CN[ std::array<Int,2>{{ 222,3 }} ] = 82;

  CN[ std::array<Int,2>{{ 223,0 }} ] = 43;
  CN[ std::array<Int,2>{{ 223,1 }} ] = 65;
  CN[ std::array<Int,2>{{ 223,2 }} ] = 76;
  CN[ std::array<Int,2>{{ 223,3 }} ] = 82;

  CN[ std::array<Int,2>{{ 231,0 }} ] = 49;
  CN[ std::array<Int,2>{{ 231,1 }} ] = 69;
  CN[ std::array<Int,2>{{ 231,2 }} ] = 79;
  CN[ std::array<Int,2>{{ 231,3 }} ] = 84;

  CN[ std::array<Int,2>{{ 241,0 }} ] = 61;
  CN[ std::array<Int,2>{{ 241,1 }} ] = 73;
  CN[ std::array<Int,2>{{ 241,2 }} ] = 81;
  CN[ std::array<Int,2>{{ 241,3 }} ] = 84;

  CN[ std::array<Int,2>{{ 242,0 }} ] = 61;
  CN[ std::array<Int,2>{{ 242,1 }} ] = 73;
  CN[ std::array<Int,2>{{ 242,2 }} ] = 81;
  CN[ std::array<Int,2>{{ 242,3 }} ] = 84;

  CN[ std::array<Int,2>{{ 243,0 }} ] = 61;
  CN[ std::array<Int,2>{{ 243,1 }} ] = 73;
  CN[ std::array<Int,2>{{ 243,2 }} ] = 81;
  CN[ std::array<Int,2>{{ 243,3 }} ] = 84;

  CN[ std::array<Int,2>{{ 244,0 }} ] = 43;
  CN[ std::array<Int,2>{{ 244,1 }} ] = 65;
  CN[ std::array<Int,2>{{ 244,2 }} ] = 76;
  CN[ std::array<Int,2>{{ 244,3 }} ] = 82;

  CN[ std::array<Int,2>{{ 311,0 }} ] = 36;
  CN[ std::array<Int,2>{{ 311,1 }} ] = 60;
  CN[ std::array<Int,2>{{ 311,2 }} ] = 73;
  CN[ std::array<Int,2>{{ 311,3 }} ] = 79;

  CN[ std::array<Int,2>{{ 312,0 }} ] = 36;
  CN[ std::array<Int,2>{{ 312,1 }} ] = 60;
  CN[ std::array<Int,2>{{ 312,2 }} ] = 73;
  CN[ std::array<Int,2>{{ 312,3 }} ] = 79;

  CN[ std::array<Int,2>{{ 313,0 }} ] = 36;
  CN[ std::array<Int,2>{{ 313,1 }} ] = 60;
  CN[ std::array<Int,2>{{ 313,2 }} ] = 73;
  CN[ std::array<Int,2>{{ 313,3 }} ] = 79;

  CN[ std::array<Int,2>{{ 321,0 }} ] = 49;
  CN[ std::array<Int,2>{{ 321,1 }} ] = 69;
  CN[ std::array<Int,2>{{ 321,2 }} ] = 79;
  CN[ std::array<Int,2>{{ 321,3 }} ] = 84;

  CN[ std::array<Int,2>{{ 322,0 }} ] = 49;
  CN[ std::array<Int,2>{{ 322,1 }} ] = 69;
  CN[ std::array<Int,2>{{ 322,2 }} ] = 79;
  CN[ std::array<Int,2>{{ 322,3 }} ] = 84;

  CN[ std::array<Int,2>{{ 323,0 }} ] = 35;
  CN[ std::array<Int,2>{{ 323,1 }} ] = 56;
  CN[ std::array<Int,2>{{ 323,2 }} ] = 70;
  CN[ std::array<Int,2>{{ 323,3 }} ] = 77;

  CN[ std::array<Int,2>{{ 324,0 }} ] = 35;
  CN[ std::array<Int,2>{{ 324,1 }} ] = 56;
  CN[ std::array<Int,2>{{ 324,2 }} ] = 70;
  CN[ std::array<Int,2>{{ 324,3 }} ] = 77;

  CN[ std::array<Int,2>{{ 331,0 }} ] = 46;
  CN[ std::array<Int,2>{{ 331,1 }} ] = 65;
  CN[ std::array<Int,2>{{ 331,2 }} ] = 77;
  CN[ std::array<Int,2>{{ 331,3 }} ] = 82;

  CN[ std::array<Int,2>{{ 332,0 }} ] = 96;
  CN[ std::array<Int,2>{{ 332,1 }} ] = 96;
  CN[ std::array<Int,2>{{ 332,2 }} ] = 96;
  CN[ std::array<Int,2>{{ 332,3 }} ] = 96;

  CN[ std::array<Int,2>{{ 333,0 }} ] = 63;
  CN[ std::array<Int,2>{{ 333,1 }} ] = 77;
  CN[ std::array<Int,2>{{ 333,2 }} ] = 85;
  CN[ std::array<Int,2>{{ 333,3 }} ] = 88;

  CN[ std::array<Int,2>{{ 334,0 }} ] = 63;
  CN[ std::array<Int,2>{{ 334,1 }} ] = 77;
  CN[ std::array<Int,2>{{ 334,2 }} ] = 85;
  CN[ std::array<Int,2>{{ 334,3 }} ] = 88;

  CN[ std::array<Int,2>{{ 335,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 335,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 335,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 335,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 411,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 411,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 411,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 411,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 412,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 412,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 412,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 412,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 421,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 421,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 421,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 421,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 422,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 422,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 422,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 422,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 423,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 423,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 423,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 423,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 511,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 511,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 511,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 511,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 512,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 512,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 512,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 512,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 521,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 521,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 521,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 521,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 522,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 522,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 522,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 522,3 }} ] = 98;

  CN[ std::array<Int,2>{{ 523,0 }} ] = 98;
  CN[ std::array<Int,2>{{ 523,1 }} ] = 98;
  CN[ std::array<Int,2>{{ 523,2 }} ] = 98;
  CN[ std::array<Int,2>{{ 523,3 }} ] = 98;

  return CN;
}



Raster::Raster (const std::string& file)
{
  std::cout << "Reading file, " << file << std::endl;
  
  std::vector<Eigen::Triplet<Real> > cc;

  std::ifstream ff ( file );
  if ( ff.is_open() )
    {
      std::string str;

      ff >> str;
      if ( str != "ncols" && str != "NCOLS" )
        {
          std::cout << "Wrong file format in the Raster .txt files" << std::endl;
          exit( -1. );
        }
      ff >> ncols;

      //std::cout << ncols << std::endl;

      ff >> str;
      if ( str != "nrows" && str != "NROWS" )
        {
          std::cout << "Wrong file format in the Raster .txt files" << std::endl;
          exit( -1. );
        }
      ff >> nrows;



      ff >> str;
      if ( str != "xllcorner" && str != "XLLCORNER" )
        {
          std::cout << "Wrong file format in the Raster .txt files" << std::endl;
          exit( -1. );
        }
      ff >> xllcorner;


      ff >> str;
      if ( str != "yllcorner" && str != "YLLCORNER" )
        {
          std::cout << "Wrong file format in the Raster .txt files" << std::endl;
          exit( -1. );
        }
      ff >> yllcorner;


      ff >> str;
      if ( str != "cellsize" && str != "CELLSIZE" )
        {
          std::cout << "Wrong file format in the Raster .txt files" << std::endl;
          exit( -1. );
        }
      ff >> cellsize;


      ff >> str;
      if ( str != "NODATA_value" )
        {
          std::cout << "Wrong file format in the Raster .txt files" << std::endl;
          exit( -1. );
        }
      ff >> NODATA_value;


      cc.reserve( nrows * ncols );


      Real value;
      for ( UInt i = 0; i < nrows; i++ )
        {
          for ( UInt j = 0; j < ncols; j++ )
            {

              ff >> value;

              cc.push_back( Eigen::Triplet<Real>( i, j, value ) );


            }
        }
    }
  else
    {
      std::cout << "Unable to open the file, check " << file << " in the SMARTSED_input file" << std::endl;
      exit( -1. );
    }

  Coords.resize (nrows, ncols);
  Coords.setFromTriplets (cc.begin(), cc.end());

};


Real
signum (const Real& x)
{ return ((x > 0) ? 1.0 : (x < 0) ? -1.0 : 0.0); }


Rain::Rain (const std::string& infiltrationModel,
            const UInt&        N,
            const bool&        isInitialLoss,
            const Real&        perc_initialLoss)
{
  M_isInitialLoss = isInitialLoss;
  c = perc_initialLoss;
    
    
  if ( infiltrationModel != "None" && infiltrationModel  != "SCS-CN" )
  {
    std::cout << "Insert a valid infiltration model, STOP!" << std::endl;
    exit( 1. );
  }
  

  DP_total.resize( N );
  DP_cumulative.resize( N );
  DP_infiltrated.resize( N );
  IDW_weights.resize( N );

}


void
Rain::constant_precipitation (const std::string&       file,
                              const UInt&              ndata,
                              const bool&              is_precipitation,
                              const Real&              time_spacing)
{
  M_time_spacing_vect.resize( 1 );
  M_time_spacing_vect[ 0 ] = time_spacing;
  
  Hyetograph.resize( 1 );
  
  
  
  if ( is_precipitation )
{
    std::ifstream ff( file );

      if ( ff.is_open() )
        {

          std::string str;
          std::getline( ff, str );


          for ( UInt i = 0; i < ndata; i++ )
            {
              std::string str1;
              ff >> str1;

              std::string hour;
              ff >> hour;

              std::string hour_string  ( hour.begin(),   hour.begin()+2 ),
                          minute_string( hour.begin()+3, hour.begin()+5 ),
                          second_string( hour.begin()+6, hour.begin()+8 );


              const Int hour_number   = std::stoi( hour_string   ),
                        minute_number = std::stoi( minute_string ),
                        second_number = std::stoi( second_string );

              const Real hour_full    = Real(hour_number) + Real(minute_number)/60 + Real(second_number)/3600;


              Real value;
              ff >> value;

              const int rr = int(std::round(24./time_spacing));
              if ( std::abs( (i%rr)*time_spacing - hour_full ) > .1 * (i%rr)*time_spacing ) //( std::abs( ((i%rr)+1)*time_spacing - hour_full ) > .1 * ((i%rr)+1)*time_spacing )
                {
                  std::cout << str1 << " " << hour << " " << value << std::endl;
                  std::cout << "Invalid rain file" << std::endl;
                  exit( 1. );
                }
              

              // mm/h  --> m/sec.
              Hyetograph[ 0 ].push_back( value * 1.e-3 / ( time_spacing * 3600 ) );
              //ff >> value; // second pluviometer data --> not used yet

              //std::cout << str1 << " " << hour << " " << value << std::endl;

            }

          ff.close();





        }
      else
        {
          std::cout << "Unable to open the file, check rain_file in the SMARTSED_input file" << std::endl;
          exit( -1. );
        }
    }
  else
    {
      for ( UInt i = 0; i < ndata; i++ )
        {
          Hyetograph[ 0 ].push_back( 0. );
        }
    }

  for ( UInt i = 0; i < IDW_weights.size(); i++ )
    {
      IDW_weights[ i ].push_back( 1. );
    }

}


void
Rain::IDW_precipitation (const std::vector<std::string>& file_vect,
                         const std::vector<UInt>&        ndata_vec,
                         const std::vector<Real>&        time_spacing_vect,
                         const std::vector<Real>&        X,
                         const std::vector<Real>&        Y,
                         const Real&                     xllcorner,
                         const Real&                     yllcorner,
                         const Real&                     pixel_size,
                         const UInt&                     N_rows,
                         const UInt&                     N_cols,
                         const std::vector<UInt>&        idBasinVect )
{
  
  M_time_spacing_vect = time_spacing_vect;
  
  Hyetograph.resize( ndata_vec.size() );


  for ( UInt k = 0; k < file_vect.size(); k++ )
    {
      const auto file = file_vect[ k ];
      const auto time_spacing = time_spacing_vect[ k ];
      std::ifstream ff( file );

      if ( ff.is_open() )
        {

          std::string str;
          std::getline( ff, str );


          for ( UInt i = 0; i < ndata_vec[ k ]; i++ )
            {
              std::string Id_sensor;
              ff >> Id_sensor;

              std::string day;
              ff >> day;

              std::string hour;
              ff >> hour;

              std::string hour_string  ( hour.begin(),   hour.begin()+2 ),
                minute_string( hour.begin()+3, hour.begin()+5 );

              //std::cout << Id_sensor << " " << day << " " << hour << " " << hour_string << " " << minute_string << std::endl;

              const Int hour_number   = std::stoi( hour_string   ),
                        minute_number = std::stoi( minute_string );
              
              const Real hour_full    = Real(hour_number) + Real(minute_number)/60;


              Real value;
              ff >> value;

              const int rr = int(std::round(24./time_spacing));
              if ( std::abs( (i%rr)*time_spacing - hour_full ) > .1 * (i%rr)*time_spacing )
                {
                  std::cout << hour << std::endl;
                  std::cout << "Invalid rain file" << std::endl;
                  exit( 1. );
                }

              

              // mm/h  --> m/sec.
              Hyetograph[ k ].push_back( value * 1.e-3 / ( time_spacing * 3600 ) );
              //ff >> value; // second pluviometer data --> not used yet

            }

          ff.close();


        }
      else
        {
          std::cout << "Unable to open the " << k+1 << "-th rain file, check the SMARTSED_input file" << std::endl;
          exit( -1. );
        }
    }


  // fill NO_DATA values
  UInt k = 0;
  for ( auto & it : Hyetograph )
    {

      for ( UInt i = 0; i < ndata_vec[ k ]; i++ )
        {
          auto & value = it[ i ];

          if ( value < 0 ) // -999.0 in ARPA files
            {
              std::vector<Real> otherStationsRain;
              UInt kk = 0;
              for ( const auto & itt : Hyetograph )
                {

                  const Int ii = std::floor( i * ( time_spacing_vect[ k ] / time_spacing_vect[ kk ] ) );
                  if ( itt[ ii ] >= 0. )
                    {
                      otherStationsRain.push_back( itt[ ii ] );
                    }
                  kk++;
                }

              if ( otherStationsRain.size() != 0 )
                {
                  Real sum = 0;
                  for ( const auto & iter : otherStationsRain )
                    {
                      sum += iter;
                    }
                  value = sum / otherStationsRain.size();
                }
              else
                {
                  value = 0.;
                }
            }

        }
      k++;
    }


  for ( const auto & it : Hyetograph )
    {
      for ( const auto & value : it )
        {
          if ( value < 0. )
            {
              std::cout << value << " One negative value in class rain!" << std::endl;
              exit( 1. );
            }
        }
    }


  // compute distances for IDW method
  for ( const auto & IDcenter : idBasinVect )
    {
      const Int i = IDcenter / N_cols,
                j = IDcenter % N_cols;

      const Real X_cell =   j * pixel_size + xllcorner,
                 Y_cell = - i * pixel_size + yllcorner + N_rows * pixel_size;


      for ( UInt ii = 0; ii < X.size(); ii++ )
        {
          const Real p = 3.;

          // divide by 1000 to obtain better number, the result does not change
          const auto delta_x = ( X_cell - X[ ii ] )/1000,
                     delta_y = ( Y_cell - Y[ ii ] )/1000;

          const auto dist = std::pow( std::sqrt( std::pow( delta_x, 2 ) + std::pow( delta_y, 2 ) ), p );


          if ( dist == 0 )
            {
              IDW_weights[ IDcenter ].push_back( 1. );
            }
          else
            {
              IDW_weights[ IDcenter ].push_back( 1./dist );
            }

        }

      Real sum = 0;
      for ( const auto & iter : IDW_weights[ IDcenter ] )
        {
          sum += iter;
        }

      for ( UInt ii = 0; ii < IDW_weights[ IDcenter ].size(); ii++ )
        {
          IDW_weights[ IDcenter ][ ii ] /= sum;
        }


    }



}


void
Rain::computePrecipitation (const Real&              time,
                            const std::vector<Real>& S,
                            const std::vector<Real>& melt_mask,
                            const std::vector<Real>& h_G,
                            const Eigen::VectorXd&   H,
                            const UInt&              N_rows,
                            const UInt&              N_cols,
                            const std::vector<UInt>& idBasinVect )
{

  // SCS-CN method and Initial and Constant Loss Model
  for ( UInt Id = 0; Id < Hyetograph.size(); Id++ )
    {
            
            
      const UInt i_index = std::floor( time / (M_time_spacing_vect[ Id ]*3600) );
            

      for ( const auto & IDcenter : idBasinVect )
        {
                
          if ( Id == 0 )
            {
              DP_total      [ IDcenter ] = 0.;
              DP_cumulative [ IDcenter ] = 0.;
              DP_infiltrated[ IDcenter ] = 0.;
            }

          rainfall_intensity = Hyetograph[ Id ][ i_index ] * IDW_weights[ IDcenter ][ Id ];
                
                        
          const Real deltaSoilMoisture = h_G[ IDcenter ] - S[ IDcenter ];
                
          Real weight = 0.;
          if ( S[ IDcenter ] > 0 && deltaSoilMoisture < 0. )
            {
              weight = std::pow( deltaSoilMoisture / S[ IDcenter ], 2. );
            }
                
                
                
          if ( weight > 1. || h_G[ IDcenter ] < 0 )
            {
              std::cout << "Error in weight infiltration model\n" <<
                "h_G = " << h_G[ IDcenter ] << "\nweight = " << weight <<
                "\nmax soil moisture ret. = " << S[ IDcenter ] << std::endl;
              exit( -1 );
            }
                

                
                
          Real infiltrationRate = weight * rainfall_intensity * melt_mask[ IDcenter ],
            potential_runoff = std::max( rainfall_intensity * melt_mask[ IDcenter ] - infiltrationRate, 0. );
                
                
          if ( ( ( H( IDcenter ) + h_G[ IDcenter ] ) < ( c * S[ IDcenter ] ) ) && M_isInitialLoss ) // initial loss
            {
              potential_runoff = 0.;
              infiltrationRate = rainfall_intensity * melt_mask[ IDcenter ];
            }
                
          //std::cout << potential_runoff << std::endl;
                
          DP_total      [ IDcenter ] += rainfall_intensity;
          DP_cumulative [ IDcenter ] += potential_runoff;
          DP_infiltrated[ IDcenter ] += infiltrationRate;
                
                
        }
        
    }
                    
}



Temperature::Temperature (const std::string&       file,
                          const UInt&              N,
                          const UInt&              max_Days,
                          const Real&              T_crit,
                          const std::vector<Real>& orography,
                          const UInt&              ndata,
                          const UInt&              steps_per_hour,
                          const Real&              time_spacing,
                          const Real&              height_thermometer,
                          const std::string        format_temp)
: T_crit( T_crit ),
  height_th( height_thermometer )
{

  T_raster          .resize( N );
  melt_mask         .resize( N );
  
  T_dailyMean       .resize( max_Days );
  T_dailyMin        .resize( max_Days );
  T_dailyMax        .resize( max_Days );
  J                 .resize( max_Days );
  
  Temperature_Graph .reserve( ndata );
  
  std::vector<Real> J_ndata;
  J_ndata.reserve( ndata );


  std::ifstream ff( file );
  
  

  if ( ff.is_open() )
    {

      if ( format_temp == "comune" )
        {
          std::string str;
          std::getline( ff, str );


          for ( UInt i = 0; i < ndata; i++ )
            {

              std::string str1;
              ff >> str1;

              std::vector<UInt> nn;
              
              if ( str1.length() != 10 )
                {
                  std::cout << str1 << std::endl;
                  std::cout << "Wrong Temperature file format" << std::endl;
                  exit( 1. );
                }


              std::string day_string  ( str1.begin(),   str1.begin()+2 ),
                          month_string( str1.begin()+3, str1.begin()+5 ),
                          year_string ( str1.begin()+6, str1.end()     );


              //std::cout << day_string << " " << month_string << " " << year_string << std::endl;



              const Int day   = std::stoi( day_string   ),
                        month = std::stoi( month_string ),
                        year  = std::stoi( year_string  );


              //std::cout << day << " " << month << " " << year << std::endl;
              //std::cout << day_string << " " << month_string << " " << year_string << std::endl;




              for ( UInt uu = 1; uu < month; uu++ )
                {
                  if ( uu == 4 || uu == 6 || uu == 9 || uu == 11 )
                    {
                      nn.push_back( 30 );
                    }
                  else if ( uu != 2 )
                    {
                      nn.push_back( 31 );
                    }
                  else if ( ( year % 4 == 0 && year % 100 != 0 ) || year % 400 == 0 ) // febbraio bisestile
                    {
                      nn.push_back( 29 );
                    }
                  else // febbraio non bisestile
                    {
                      nn.push_back( 28 );
                    }
                }

              Real scalar_result = 0;
              for ( UInt uu = 0; uu < nn.size(); uu++ )
                {
                  scalar_result += nn[ uu ];
                }

              J_ndata.push_back( day + scalar_result );

              //std::cout << day + scalar_result << std::endl;

              std::string hour;
              ff >> hour;

              std::string hour_string( hour.begin(),   hour.begin()+2 ),
                minute_string( hour.begin()+3, hour.begin()+5 ),
                second_string( hour.begin()+6, hour.begin()+8 );


              const Int hour_number   = std::stoi( hour_string   ),
                        minute_number = std::stoi( minute_string ),
                        second_number = std::stoi( second_string );

              const Real hour_full    = Real(hour_number) + Real(minute_number)/60 + Real(second_number)/3600;


              Real value;  // temperature data
              ff >> value;
              
              
              
              const int rr = int(std::round(24./time_spacing));
              if ( std::abs( (i%rr)*time_spacing - hour_full ) > .1 * (i%rr)*time_spacing ) //( std::abs( ((i%rr)+1)*time_spacing - hour_full ) > .1 * ((i%rr)+1)*time_spacing )
                {
                  std::cout << str1 << " " << hour << " " << value << " " << i+1 << std::endl;
                  std::cout << "Invalid temperature file, maybe check time spacing in SMARTSED_input file" << std::endl;
                  exit( 1. );
                }

              if ( value == -999 && i == 0 )
                {
                  std::cout << "Correct temperature file to eliminate first element as NODATA value" << std::endl;
                  exit( -1. );
                }

              if ( value == -999 )
                {
                  value = Temperature_Graph[ i - 1 ];
                }



              Temperature_Graph.push_back( value );


            }

          ff.close();
        }
      else  if ( format_temp == "arpa" )
        {

          std::string str;
          std::getline( ff, str );


          for ( UInt i = 0; i < ndata; i++ )
            {

              UInt Id_sensor;
              ff >> Id_sensor;

              std::string str1;
              ff >> str1;

              std::vector<UInt> nn;


              if ( str1.length() != 10 )
                {
                  std::cout << str1 << std::endl;
                  std::cout << "Wrong Temperature file format" << std::endl;
                  exit( 1. );
                }


              std::string year_string ( str1.begin(),   str1.begin()+4 ),
                month_string( str1.begin()+5, str1.begin()+7 ),
                day_string  ( str1.begin()+8, str1.end()     );




              const Int day   = std::stoi( day_string   ),
                        month = std::stoi( month_string ),
                        year  = std::stoi( year_string  );


              //std::cout << day << " " << month << " " << year << std::endl;
              //std::cout << day_string << " " << month_string << " " << year_string << std::endl;




              for ( UInt uu = 1; uu < month; uu++ )
                {
                  if ( uu == 4 || uu == 6 || uu == 9 || uu == 11 )
                    {
                      nn.push_back( 30 );
                    }
                  else if ( uu != 2 )
                    {
                      nn.push_back( 31 );
                    }
                  else if ( ( year % 4 == 0 && year % 100 != 0 ) || year % 400 == 0 ) // febbraio bisestile
                    {
                      nn.push_back( 29 );
                    }
                  else // febbraio non bisestile
                    {
                      nn.push_back( 28 );
                    }
                }

              Real scalar_result = 0;
              for ( UInt uu = 0; uu < nn.size(); uu++ )
                {
                  scalar_result += nn[ uu ];
                }

              J_ndata.push_back( day + scalar_result );

              //std::cout << day + scalar_result << std::endl;

              std::string hour;
              ff >> hour;



              std::string hour_string  ( hour.begin(),   hour.begin()+2 ),
                minute_string( hour.begin()+3, hour.begin()+5 );


              const Int hour_number   = std::stoi( hour_string   ),
                minute_number = std::stoi( minute_string );

              const Real hour_full    = Real(hour_number) + Real(minute_number)/60;


              Real value;  // temperature data
              ff >> value;

              //std::cout << str1 << " " << hour << " " << value << std::endl;

              const int rr = int(std::round(24./time_spacing));
              if ( std::abs( (i%rr)*time_spacing - hour_full ) > .1 * (i%rr)*time_spacing )
                {
                  std::cout << str1 << " " << hour << " " << value << std::endl;
                  std::cout << "Invalid temperature file" << std::endl;
                  exit( 1. );
                }

              if ( value == -999 && i == 0 )
                {
                  std::cout << "Correct temperature file to eliminate first element as NODATA value" << std::endl;
                  exit( -1. );
                }

              if ( value == -999 )
                {
                  value = Temperature_Graph[ i - 1 ];
                }


              Temperature_Graph.push_back( value );


            }

          ff.close();
        }
      else
        {
          std::cout << "Temperature format non recognized" << std::endl;
          exit( -1. );
        }

    }
  else
    {
      std::cout << "Unable to open the file, check temperature_file in the SMARTSED_input file" << std::endl;
      exit( -1. );
    }


  // --------------------------------------------- //
  for ( UInt n = 1; n <= max_Days; n++ )
    {

      const auto i = std::floor( ( n - 1 ) * ( 24./time_spacing ) );

      UInt k = 0,
        h = i; // i


      while ( k != UInt( std::round( (24./time_spacing) ) ) )
        {

          T_dailyMean[ n - 1 ] += Temperature_Graph[ h ];


          if ( k != 0 )
            {

              if ( Temperature_Graph[ h ] < T_dailyMin[ n - 1 ] )
                {
                  T_dailyMin[ n - 1 ] = Temperature_Graph[ h ];
                }

              if ( Temperature_Graph[ h ] > T_dailyMax[ n - 1 ] )
                {
                  T_dailyMax[ n - 1 ] = Temperature_Graph[ h ];
                }
            }
          else
            {
              T_dailyMin [ n - 1 ] = Temperature_Graph[ h ];
              T_dailyMax [ n - 1 ] = Temperature_Graph[ h ];
            }

          h++;
          k++;
        }

      if ( k == 0 )
        {
          std::cout << "Something wrong in Temperature class constructor" << std::endl;
          exit ( -1. );
        }

      T_dailyMean[ n - 1 ] /= k;

    }
  // --------------------------------------------- //


  // Now fill J starting from J_ndata
  for ( UInt n = 1; n <= max_Days; n++ )
    {
      const UInt i = std::floor( ( n - 1 ) * ( 24./time_spacing ) );
      J[ n - 1 ] = J_ndata[ i ];
    }


}

void
Temperature::computeTemperature (const UInt&              i,
                                 const std::vector<Real>& orography,
                                 const std::vector<UInt>& idBasinVect )
{
  const Real T = Temperature_Graph[ i ];

  if ( std::isnan(T) )
    {
      std::cout << "nan in computeTemperature" << std::endl;
      exit( -1 );
    }

  for ( const auto & j : idBasinVect )
    {
      T_raster [ j ] = T + Temp_diff * ( orography[ j ] - height_th );
      melt_mask[ j ] = ( T_raster[ j ] > T_crit );    // melt_mask = 1 -\mu in the paper
    }
}



evapoTranspiration::evapoTranspiration (const std::string&       ET_model,
                                        const UInt&              N,
                                        const std::vector<Real>& orography,
                                        const std::vector<Real>& J,
                                        const UInt&              max_Days,
                                        const Real&              phi_rad,
                                        const Real&              height_thermometer)
  : height_th (height_thermometer)
{

  ET_vec.resize( N );

  if ( ET_model == "None" )
    {
      M_evapoTranspiration_model = 0;
    }
  else if ( ET_model == "Hargreaves" )
    {
      M_evapoTranspiration_model = 1;
    }
  else
    {
      std::cout << "No evapo-transpiration model inserted!!" << std::endl;
      exit( -1. );
    }


  Ra.resize( max_Days );
  for ( UInt n = 1; n <= max_Days; n++ )
    {
      const auto dr    = 1 +.033 * std::cos( 2 * M_PI * J[ n - 1 ] / 365 ),
        delta =    .409 * std::sin( 2 * M_PI * J[ n - 1 ] / 365 - 1.39 ),
        ws    = std::acos( - std::tan( phi_rad ) * std::tan( delta ) );

      Ra[ n - 1 ] = (24*60/M_PI) * M_Gsc * dr * ( ws * std::sin( phi_rad ) * std::sin( delta ) + std::cos( phi_rad ) * std::cos( delta ) * std::sin( ws ) );
    }


}

void
evapoTranspiration::ET (const std::vector<Real>& T_mean, // lungo nstep: vettore delle temperature in Â°C
                        const std::vector<Real>& T_min,  // lungo nstep
                        const std::vector<Real>& T_max,  // lungo nstep
                        const Int&               i,
                        const std::vector<UInt>& idBasinVect,
                        const std::vector<Real>& orography)
{

  switch ( M_evapoTranspiration_model )
    {
    case 0:

      break;

    case 1:

      for ( const auto & k : idBasinVect )
        {

          const auto t_mean = T_mean[ i ] + Temp_diff * ( orography[ k ] - height_th ),
            t_max  = T_max [ i ] + Temp_diff * ( orography[ k ] - height_th ),
            t_min  = T_min [ i ] + Temp_diff * ( orography[ k ] - height_th );

          // unity: mm/day --> m/sec.
          ET_vec[ k ] = .0023 * Ra[ i ] * ( t_mean + 17.8 ) * std::pow( ( t_max - t_min ), .5 ) * ( 1.e-3/(24*3600) ); // più che altro mettere la T del giorno media per poter calcolare T_min e T_max

        }

      break;
    }
}


frictionClass::frictionClass (const std::vector<Real>& H_interface_horizontal, 
                              const std::vector<Real>& H_interface_vertical, 
                              const std::vector<Real>& u, const std::vector<Real>& v, 
                              const std::vector<UInt>& idStaggeredInternalVectHorizontal,
                              const std::vector<UInt>& idStaggeredBoundaryVectWest,
                              const std::vector<UInt>& idStaggeredBoundaryVectEast,
                              const std::vector<UInt>& idStaggeredInternalVectVertical,
                              const std::vector<UInt>& idStaggeredBoundaryVectNorth,
                              const std::vector<UInt>& idStaggeredBoundaryVectSouth,
                              const std::string& friction_model,
                              const Real& n_manning,
                              const Real& dt_DSV,
                              const std::vector<Real>& d_90,
                              const std::vector<Real>& rough,
                              const Real& H_min,
                              const UInt& N_rows,
                              const UInt& N_cols,
                              const std::vector<Real>& S_x,
                              const std::vector<Real>& S_y)
: H_interface_horizontal(H_interface_horizontal), H_interface_vertical(H_interface_vertical), u(u), v(v),
  idStaggeredInternalVectHorizontal(idStaggeredInternalVectHorizontal), idStaggeredBoundaryVectWest(idStaggeredBoundaryVectWest), idStaggeredBoundaryVectEast(idStaggeredBoundaryVectEast),
  idStaggeredInternalVectVertical(idStaggeredInternalVectVertical), idStaggeredBoundaryVectNorth(idStaggeredBoundaryVectNorth), idStaggeredBoundaryVectSouth(idStaggeredBoundaryVectSouth),
  M_n_manning( n_manning ),
  M_dt_DSV( dt_DSV ),
  N_rows( N_rows ),
  N_cols( N_cols )
{

  M_H_min =  std::pow(H_min, M_expo);

  std::vector<Real> M_fc0_lower_x( S_x.size() ),
    M_fc0_greater_x( S_x.size() ),
    M_fc0_lower_y( S_y.size() ),
    M_fc0_greater_y( S_y.size() );

  for ( UInt i = 0; i < N_rows; i++ )
    {
      for ( UInt j = 0; j < N_cols; j++ )
        {
          const UInt IDcell = j + i * N_cols,
            IDleft = IDcell + i,
            IDright = IDleft + 1,
            IDup = IDcell,
            IDdown = IDcell + N_cols;

          const auto d_90_cell = rough[ IDcell ] * d_90[ IDcell ];

          M_fc0_greater_x[ IDleft ] = std::pow( d_90_cell, .45 ) / ( .56  * std::pow( M_g, .44 ) );
          M_fc0_lower_x  [ IDleft ] = std::pow( d_90_cell, .24 ) / ( 2.73 * std::pow( M_g, .49 ) );
          M_fc0_greater_y[ IDup ]   = std::pow( d_90_cell, .45 ) / ( .56  * std::pow( M_g, .44 ) );
          M_fc0_lower_y  [ IDup ]   = std::pow( d_90_cell, .24 ) / ( 2.73 * std::pow( M_g, .49 ) );

          M_fc0_greater_x[ IDright ] = std::pow( d_90_cell, .45 ) / ( .56  * std::pow( M_g, .44 ) );
          M_fc0_lower_x  [ IDright ] = std::pow( d_90_cell, .24 ) / ( 2.73 * std::pow( M_g, .49 ) );
          M_fc0_greater_y[ IDdown ]  = std::pow( d_90_cell, .45 ) / ( .56  * std::pow( M_g, .44 ) );
          M_fc0_lower_y  [ IDdown ]  = std::pow( d_90_cell, .24 ) / ( 2.73 * std::pow( M_g, .49 ) );
        }
    }




  alfa_x.resize( S_x.size() );
  alfa_y.resize( S_y.size() );

  M_coeff = M_g * std::pow( M_n_manning, 2. );


  M_expo_r_x_vect.resize( S_x.size() );
  M_gamma_dt_DSV_x_.resize( S_x.size() );
  for ( UInt k = 0; k < S_x.size(); k++ )
    {
        M_expo_r_x_vect[ k ] = M_expo_r1 * ( std::abs( S_x[ k ] ) >   .006 ) +
        M_expo_r2 * ( std::abs( S_x[ k ] ) <=  .006 );

        const auto M_Rick_x = ( M_fc0_greater_x[ k ] * std::pow( std::abs( S_x[ k ] ), .33 ) ) * ( std::abs( S_x[ k ] ) >  .006 ) +
          ( M_fc0_lower_x[ k ]   * std::pow( std::abs( S_x[ k ] ), .08 ) ) * ( std::abs( S_x[ k ] ) <= .006 );
        M_gamma_dt_DSV_x_[ k ] = M_g * std::pow( M_Rick_x, 2 );

    }

  M_expo_r_y_vect.resize( S_y.size() );
  M_gamma_dt_DSV_y_.resize( S_y.size() );
  for ( UInt k = 0; k < S_y.size(); k++ )
    {
        M_expo_r_y_vect[ k ] = M_expo_r1 * ( std::abs( S_y[ k ] ) >   .006 ) +
        M_expo_r2 * ( std::abs( S_y[ k ] ) <=  .006 );

        const auto M_Rick_y = ( M_fc0_greater_y[ k ] * std::pow( std::abs( S_y[ k ] ), .33 ) ) * ( std::abs( S_y[ k ] ) >  .006 ) +
          ( M_fc0_lower_y[ k ]   * std::pow( std::abs( S_y[ k ] ), .08 ) ) * ( std::abs( S_y[ k ] ) <= .006 );
        M_gamma_dt_DSV_y_ [ k ] = M_g * std::pow( M_Rick_y, 2 );

    }

  if ( friction_model == "None" )
    {
      M_frictionModel = 0;
    }
  else if ( friction_model == "Manning" )
    {
      M_frictionModel = 1;
    }
  else if ( friction_model == "Rickenmann" )
    {
      M_frictionModel = 2;
    }
  else
    {
      std::cout << "No friction model inserted!!" << std::endl;
      exit( -1. );
    }
}

void
frictionClass::f_x ()
{

  
  for ( const auto & Id : idStaggeredInternalVectHorizontal )
    {
        Real alfa = 1.;

        const auto & H_int = H_interface_horizontal[ Id ];
        const auto & exponent = M_expo_r_x_vect[ Id ];
        const auto den = std::pow( H_int, M_expo + exponent * (M_frictionModel == 2) );

        if ( den > M_H_min )
          {

              const auto u_abs = std::abs( u[ Id ] );

              Real coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * u_abs / den * (M_frictionModel > 0);
              coeff = std::max( coeff, M_dt_DSV * M_gamma_dt_DSV_x_[ Id ] * std::pow( u_abs, 1. - exponent * (M_frictionModel == 2) ) / den );
              alfa = 1. / ( 1. + coeff );
          }

        alfa_x[ Id ] = alfa;

    }

    for ( const auto & Id : idStaggeredBoundaryVectWest )
      {
          Real alfa = 1.;

          const auto & H_int = H_interface_horizontal[ Id ];
          const auto & exponent = M_expo_r_x_vect[ Id ];
          const auto den = std::pow( H_int, M_expo + exponent * (M_frictionModel == 2) );

          if ( den > M_H_min )
            {

                const auto u_abs = std::abs( u[ Id ] );

                Real coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * u_abs / den * (M_frictionModel > 0);
                coeff = std::max( coeff, M_dt_DSV * M_gamma_dt_DSV_x_[ Id ] * std::pow( u_abs, 1. - exponent * (M_frictionModel == 2) ) / den );
                alfa = 1. / ( 1. + coeff );
            }

          alfa_x[ Id ] = alfa;

      }

    for ( const auto & Id : idStaggeredBoundaryVectEast )
      {
          Real alfa = 1.;

          const auto & H_int = H_interface_horizontal[ Id ];
          const auto & exponent = M_expo_r_x_vect[ Id ];
          const auto den = std::pow( H_int, M_expo + exponent * (M_frictionModel == 2) );

          if ( den > M_H_min )
            {

                const auto u_abs = std::abs( u[ Id ] );

                Real coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * u_abs / den * (M_frictionModel > 0);
                coeff = std::max( coeff, M_dt_DSV * M_gamma_dt_DSV_x_[ Id ] * std::pow( u_abs, 1. - exponent * (M_frictionModel == 2) ) / den );
                alfa = 1. / ( 1. + coeff );
            }

          alfa_x[ Id ] = alfa;

      }

}


void
frictionClass::f_y ()
{

    for ( const auto & Id : idStaggeredInternalVectVertical )
      {
          Real alfa = 1.;

          const auto & H_int = H_interface_vertical[ Id ];
          const auto & exponent = M_expo_r_y_vect[ Id ];
          const auto den = std::pow( H_int, M_expo + exponent * (M_frictionModel == 2) );

          if ( den > M_H_min )
            {

                const auto v_abs = std::abs( v[ Id ] );

                Real coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * v_abs / den * (M_frictionModel > 0);
                coeff = std::max( coeff, M_dt_DSV * M_gamma_dt_DSV_y_[ Id ] * std::pow( v_abs, 1. - exponent * (M_frictionModel == 2) ) / den );
                alfa = 1. / ( 1. + coeff );
            }


          alfa_y[ Id ] = alfa;

      }

      for ( const auto & Id : idStaggeredBoundaryVectNorth )
        {
            Real alfa = 1.;

            const auto & H_int = H_interface_vertical[ Id ];
            const auto & exponent = M_expo_r_y_vect[ Id ];
            const auto den = std::pow( H_int, M_expo + exponent * (M_frictionModel == 2) );

            if ( den > M_H_min )
              {

                  const auto v_abs = std::abs( v[ Id ] );

                  Real coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * v_abs / den * (M_frictionModel > 0);
                  coeff = std::max( coeff, M_dt_DSV * M_gamma_dt_DSV_y_[ Id ] * std::pow( v_abs, 1. - exponent * (M_frictionModel == 2) ) / den );
                  alfa = 1. / ( 1. + coeff );
              }


            alfa_y[ Id ] = alfa;

        }

      for ( const auto & Id : idStaggeredBoundaryVectSouth )
        {
            Real alfa = 1.;

            const auto & H_int = H_interface_vertical[ Id ];
            const auto & exponent = M_expo_r_y_vect[ Id ];
            const auto den = std::pow( H_int, M_expo + exponent * (M_frictionModel == 2) );

            if ( den > M_H_min )
              {

                  const auto v_abs = std::abs( v[ Id ] );

                  Real coeff = M_gamma_dt_DSV(M_dt_DSV, M_coeff) * v_abs / den * (M_frictionModel > 0);
                  coeff = std::max( coeff, M_dt_DSV * M_gamma_dt_DSV_y_[ Id ] * std::pow( v_abs, 1. - exponent * (M_frictionModel == 2) ) / den );
                  alfa = 1. / ( 1. + coeff );
              }


            alfa_y[ Id ] = alfa;

        }

}

void
upwind::computeHorizontal ()
{

  
  for ( const auto & Id : idStaggeredInternalVectHorizontal )
    {
      const UInt i      = Id / ( N_cols + 1 ), // u
        IDeast = Id - i,              // H
        IDwest = Id - i - 1;          // H


      const Real & H_left  = H( IDwest ),
        & H_right = H( IDeast );
      
       
      horizontal[ Id ] = ( H_left + H_right ) * .5 + signum( u[ Id ] ) * ( H_left - H_right ) * .5;
    }



  for ( const auto & Id : idStaggeredBoundaryVectWest )
    {

      const UInt i = Id / ( N_cols + 1 );



      const Real   H_left  = 0,
        & H_right = H( Id - i );


      horizontal[ Id ] = ( H_left + H_right ) * .5 + signum( u[ Id + 1 ] ) * ( H_left - H_right ) * .5; // messo il segno della velocitÃ  interna di fianco

    }


  for ( const auto & Id : idStaggeredBoundaryVectEast )
    {

      const UInt i = Id / ( N_cols + 1 );



      const Real & H_left  = H( Id - i - 1 ),
        H_right = 0;


      horizontal[ Id ] = ( H_left + H_right ) * .5 + signum( u[ Id - 1 ] ) * ( H_left - H_right ) * .5; // messo il segno della velocitÃ  interna di fianco


    }




}

void
upwind::computeVertical ()
{


  for ( const auto & Id : idStaggeredInternalVectVertical )
    {

      const auto IDsouth = Id,              // H
        IDnorth = Id - N_cols;     // H



      const Real & H_left  = H( IDnorth ),
        & H_right = H( IDsouth );


      vertical[ Id ] = ( H_left + H_right ) * .5 + signum( v[ Id ] ) * ( H_left - H_right ) * .5;


    }


  for ( const auto & Id : idStaggeredBoundaryVectNorth )
    {

      const Real   H_left  = 0,
        & H_right = H( Id );


      vertical[ Id ] = ( H_left + H_right ) * .5 + signum( v[ Id + N_cols ] ) * ( H_left - H_right ) * .5;

    }


  for ( const auto & Id : idStaggeredBoundaryVectSouth )
    {

      const Real & H_left  = H( Id - N_cols ),
        H_right = 0;


      vertical[ Id ] = ( H_left + H_right ) * .5 + signum( v[ Id - N_cols ] ) * ( H_left - H_right ) * .5;


    }

}

void
bilinearInterpolation (const std::vector<Real>& u,
                       const std::vector<Real>& v,
                       std::vector<Real>& u_star,
                       std::vector<Real>& v_star,
                       const UInt&              nrows,
                       const UInt&              ncols,
                       const Real&              dt,
                       const Real&              pixel_size)
{


  // +-----------------------------------------------+
  // |              Horizontal Velocity              |
  // +-----------------------------------------------+


  for ( UInt i = 0; i < nrows; i++ )
    {

      for ( UInt j = 1; j < ncols; j++ )
        {

          const auto Id    = j + i * ( ncols + 1 ), // u

            // ID della velocitÃ
            ID_NE = Id - i,        // v
            ID_NW = ID_NE - 1,     // v
            ID_SE = ID_NE + ncols, // v
            ID_SW = ID_NW + ncols; // v


          //std::cout << i << std::endl;


          Vector2D vel( std::array<Real,2>{{ u[ Id ], ( v[ ID_NE ] + v[ ID_NW ] + v[ ID_SE ] + v[ ID_SW ] ) / 4. }} );

          const auto Dx = vel / pixel_size * dt;

          Vector2D xx( std::array<Real,2>{{ Real( j ), Real( i ) }} );  // Top-Left reference frame

          xx = xx - Dx;

          auto x = xx( 0 ),
            y = xx( 1 );

          auto x_1 = std::floor( x ),           //  11    21     ----> x
            y_1 = std::floor( y ),              //
            x_2 = x_1 + 1,                      //
            y_2 = y_1 + 1;                      //  12    22
          //
          // |
          // | y

          if ( x_1 == ncols )
            {
              x_2 -= 1;
              x_1 -= 1;
              x   -= 1;
            }
          if ( x_2 == 0 )
            {
              x_2 += 1;
              x_1 += 1;
              x   += 1;
            }


          if ( y_1 == ( nrows - 1 ) )
            {
              y_2 -= 1;
              y_1 -= 1;
              y   -= 1;
            }
          if ( y_2 == 0 )
            {
              y_2 += 1;
              y_1 += 1;
              y   += 1;
            }

          // return 0 for target values that are out of bounds
          if ( x_2 < 0 || x_1 > ncols || y_2 < 0 || y_1 > ( nrows - 1 ) )
            {
              u_star[ Id ] = 0.;
            }
          else
            {

              const auto Id_11 = x_1 + y_1 * ( ncols + 1 ), // u
                Id_12 = x_1 + y_2 * ( ncols + 1 ), // u
                Id_21 = x_2 + y_1 * ( ncols + 1 ), // u
                Id_22 = x_2 + y_2 * ( ncols + 1 ); // u

              // compute weights
              const Real w_x2 = x_2 - x,
                w_x1 = x   - x_1,
                w_y2 = y_2 - y,
                w_y1 = y   - y_1;

              const auto a = u[ Id_11 ] * w_x2 + u[ Id_21 ] * w_x1,
                b = u[ Id_12 ] * w_x2 + u[ Id_22 ] * w_x1;

              u_star[ Id ] = a * w_y2 + b * w_y1;

            }

        }

    }

  // 1st column
  {
    UInt j = 0;
    for ( UInt i = 0; i < nrows; i++ )
      {

        const auto Id  = j + i * ( ncols + 1 ),
          Idd = Id + 1;

        u_star[ Id ] = u_star[ Idd ] * ( u[ Idd ] < 0. );
        //u_star[ Id ] *= ( u[ Idd ] < 0. );

      }
  }


  // last column
  {
    UInt j = ncols;
    for ( UInt i = 0; i < nrows; i++ )
      {

        const auto Id  = j + i * ( ncols + 1 ),
          Idd = Id - 1;

        u_star[ Id ] = u_star[ Idd ] * ( u[ Idd ] > 0. );
        //u_star[ Id ] *= ( u[ Idd ] > 0. );

      }
  }






  // +-----------------------------------------------+
  // |              Vertical Velocity                |
  // +-----------------------------------------------+


  for ( UInt i = 1; i < nrows; i++ )
    {

      for ( UInt j = 0; j < ncols; j++ )
        {

          const auto Id    = j + i * ncols, // v

            // ID della velocitÃ
            ID_SW = Id + i,                // u
            ID_SE = ID_SW + 1,             // u
            ID_NW = ID_SW - ( ncols + 1 ), // u
            ID_NE = ID_NW + 1;             // u




          Vector2D vel( std::array<Real,2>{{ ( u[ ID_SW ] + u[ ID_SE ] + u[ ID_NW ] + u[ ID_NE ] ) / 4., v[ Id ] }} );



          const auto Dx = vel / pixel_size * dt;





          Vector2D xx( std::array<Real,2>{{ Real( j ), Real( i ) }} );  // Top-Left reference frame



          xx = xx - Dx;




          auto x = xx( 0 ),
            y = xx( 1 );


          auto x_1 = std::floor( x ),              //  11    21     ----> x
            y_1 = std::floor( y ),              //
            x_2 = x_1 + 1,                      //
            y_2 = y_1 + 1;                      //  12    22
          //
          // |
          // | y




          if ( x_1 == ( ncols - 1 ) )
            {
              x_2 -= 1;
              x_1 -= 1;
              x   -= 1;
            }
          if ( x_2 == 0 )
            {
              x_2 += 1;
              x_1 += 1;
              x   += 1;
            }


          if ( y_1 == nrows )
            {
              y_2 -= 1;
              y_1 -= 1;
              y   -= 1;
            }
          if ( y_2 == 0 )
            {
              y_2 += 1;
              x_2 += 1;
              y   += 1;
            }



          // return 0 for target values that are out of bounds
          if ( x_2 < 0 || x_1 > ( ncols - 1 ) || y_2 < 0 || y_1 > nrows )
            {
              v_star[ Id ] = 0.;
            }
          else
            {


              const auto Id_11 = x_1 + y_1 * ncols,
                Id_12 = x_1 + y_2 * ncols,
                Id_21 = x_2 + y_1 * ncols,
                Id_22 = x_2 + y_2 * ncols;



              // compute weights
              const Real w_x2 = x_2 - x,
                w_x1 = x   - x_1,
                w_y2 = y_2 - y,
                w_y1 = y   - y_1;


              const auto a = v[ Id_11 ] * w_x2 + v[ Id_21 ] * w_x1,
                b = v[ Id_12 ] * w_x2 + v[ Id_22 ] * w_x1;


              v_star[ Id ] = a * w_y2 + b * w_y1;

            }


        }

    }



  // 1st row
  {
    UInt i = 0;
    for ( UInt j = 0; j < ncols; j++ )
      {

        const auto Id  = j + i * ncols,
          Idd = Id + ncols;

        v_star[ Id ] = v_star[ Idd ] * ( v[ Idd ] < 0. );
        //v_star[ Id ] *= ( v[ Idd ] < 0. );

      }
  }



  // last row
  {
    UInt i = nrows;
    for ( UInt j = 0; j < ncols; j++ )
      {

        const auto Id  = j + i * ncols,
          Idd = Id - ncols;

        v_star[ Id ] = v_star[ Idd ] * ( v[ Idd ] > 0. );
        //v_star[ Id ] *= ( v[ Idd ] > 0. );

      }
  }

}

bool 
is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

void
bilinearInterpolation (const std::vector<Real>& u,
                       const std::vector<Real>& v,
                       std::vector<Real>& u_star,
                       std::vector<Real>& v_star,
                       const UInt&              nrows,
                       const UInt&              ncols,
                       const Real&              dt_DSV,
                       const Real&              pixel_size,
                       const std::vector<UInt>& idStaggeredInternalVectHorizontal,
                       const std::vector<UInt>& idStaggeredInternalVectVertical,
                       const std::vector<UInt>& idStaggeredBoundaryVectWest,
                       const std::vector<UInt>& idStaggeredBoundaryVectEast,
                       const std::vector<UInt>& idStaggeredBoundaryVectNorth,
                       const std::vector<UInt>& idStaggeredBoundaryVectSouth)
{


  // +-----------------------------------------------+
  // |              Horizontal Velocity              |
  // +-----------------------------------------------+


  for ( const auto & Id : idStaggeredInternalVectHorizontal )
    {


      const UInt i     = Id / ( ncols + 1 ),
        j     = Id % ( ncols + 1 ),
        ID_NE = Id - i,        // v
        ID_NW = ID_NE - 1,     // v
        ID_SE = ID_NE + ncols, // v
        ID_SW = ID_NW + ncols; // v



      Vector2D vel( std::array<Real,2>{{ u[ Id ], ( v[ ID_NE ] + v[ ID_NW ] + v[ ID_SE ] + v[ ID_SW ] ) / 4. }} );





      const auto Dx = vel / pixel_size * dt_DSV;




      Vector2D xx( std::array<Real,2>{{ Real( j ), Real( i ) }} );  // Top-Left reference frame



      xx = xx - Dx;




      auto x = xx( 0 ),
        y = xx( 1 );


      auto x_1 = std::floor( x ),              //  11    21     ----> x
        y_1 = std::floor( y ),              //
        x_2 = x_1 + 1,                      //
        y_2 = y_1 + 1;                      //  12    22
      //
      // |
      // | y




      if ( x_1 == ncols )
        {
          x_2 -= 1;
          x_1 -= 1;
          x   -= 1;
        }
      if ( x_2 == 0 )
        {
          x_2 += 1;
          x_1 += 1;
          x   += 1;
        }


      if ( y_1 == ( nrows - 1 ) )
        {
          y_2 -= 1;
          y_1 -= 1;
          y   -= 1;
        }
      if ( y_2 == 0 )
        {
          y_2 += 1;
          y_1 += 1;
          y   += 1;
        }



      // return 0 for target values that are out of bounds
      if ( x_2 < 0 || x_1 > ncols || y_2 < 0 || y_1 > ( nrows - 1 ) )
        {
          u_star[ Id ] = 0.;
        }
      else
        {

          const auto Id_11 = x_1 + y_1 * ( ncols + 1 ), // u
            Id_12 = x_1 + y_2 * ( ncols + 1 ), // u
            Id_21 = x_2 + y_1 * ( ncols + 1 ), // u
            Id_22 = x_2 + y_2 * ( ncols + 1 ); // u




          // compute weights
          const Real w_x2 = x_2 - x,
            w_x1 = x   - x_1,
            w_y2 = y_2 - y,
            w_y1 = y   - y_1;



          const auto a = u[ Id_11 ] * w_x2 + u[ Id_21 ] * w_x1,
            b = u[ Id_12 ] * w_x2 + u[ Id_22 ] * w_x1;


          u_star[ Id ] = a * w_y2 + b * w_y1;


        }

    }




  for ( const auto & Id : idStaggeredBoundaryVectWest )
    {

      const auto Idd = Id + 1;

      u_star[ Id ] = u_star[ Idd ] * ( u[ Idd ] < 0. );

    }


  for ( const auto & Id : idStaggeredBoundaryVectEast )
    {

      const auto Idd = Id - 1;

      u_star[ Id ] = u_star[ Idd ] * ( u[ Idd ] > 0. );

    }



  // +-----------------------------------------------+
  // |              Vertical Velocity                |
  // +-----------------------------------------------+


  for ( const auto & Id : idStaggeredInternalVectVertical )
    {


      const auto i     = Id / ncols,
        j     = Id % ncols,
        ID_SW = Id + i,                // u
        ID_SE = ID_SW + 1,             // u
        ID_NW = ID_SW - ( ncols + 1 ), // u
        ID_NE = ID_NW + 1;             // u




      Vector2D vel( std::array<Real,2>{{ ( u[ ID_SW ] + u[ ID_SE ] + u[ ID_NW ] + u[ ID_NE ] ) / 4., v[ Id ] }} );



      const auto Dx = vel / pixel_size * dt_DSV;





      Vector2D xx( std::array<Real,2>{{ Real( j ), Real( i ) }} );  // Top-Left reference frame



      xx = xx - Dx;




      auto x = xx( 0 ),
        y = xx( 1 );


      auto x_1 = std::floor( x ),              //  11    21     ----> x
        y_1 = std::floor( y ),              //
        x_2 = x_1 + 1,                      //
        y_2 = y_1 + 1;                      //  12    22
      //
      // |
      // | y




      if ( x_1 == ( ncols - 1 ) )
        {
          x_2 -= 1;
          x_1 -= 1;
          x   -= 1;
        }
      if ( x_2 == 0 )
        {
          x_2 += 1;
          x_1 += 1;
          x   += 1;
        }


      if ( y_1 == nrows )
        {
          y_2 -= 1;
          y_1 -= 1;
          y   -= 1;
        }
      if ( y_2 == 0 )
        {
          y_2 += 1;
          x_2 += 1;
          y   += 1;
        }



      // return 0 for target values that are out of bounds
      if ( x_2 < 0 || x_1 > ( ncols - 1 ) || y_2 < 0 || y_1 > nrows )
        {
          v_star[ Id ] = 0.;
        }
      else
        {


          const auto Id_11 = x_1 + y_1 * ncols,
            Id_12 = x_1 + y_2 * ncols,
            Id_21 = x_2 + y_1 * ncols,
            Id_22 = x_2 + y_2 * ncols;



          // compute weights
          const Real w_x2 = x_2 - x,
            w_x1 = x   - x_1,
            w_y2 = y_2 - y,
            w_y1 = y   - y_1;


          const auto a = v[ Id_11 ] * w_x2 + v[ Id_21 ] * w_x1,
            b = v[ Id_12 ] * w_x2 + v[ Id_22 ] * w_x1;


          v_star[ Id ] = a * w_y2 + b * w_y1;

        }




    }




  for ( const auto & Id : idStaggeredBoundaryVectNorth )
    {

      const auto Idd = Id + ncols;

      v_star[ Id ] = v_star[ Idd ] * ( v[ Idd ] < 0. );

    }

  for ( const auto & Id : idStaggeredBoundaryVectSouth )
    {

      const auto Idd = Id - ncols;

      v_star[ Id ] = v_star[ Idd ] * ( v[ Idd ] > 0. );

    }


}


void
bilinearInterpolation (const std::vector<Real>& u,
                       const std::vector<Real>& v,
                       const std::vector<Real>& slope_x,
                       const std::vector<Real>& slope_y,
                       const Real&              slope_thr,
                       std::vector<Real>& u_star,
                       std::vector<Real>& v_star,
                       const UInt&              nrows,
                       const UInt&              ncols,
                       const Real&              dt_DSV,
                       const Real&              pixel_size,
                       const std::vector<UInt>& idStaggeredInternalVectHorizontal,
                       const std::vector<UInt>& idStaggeredInternalVectVertical,
                       const std::vector<UInt>& idStaggeredBoundaryVectWest,
                       const std::vector<UInt>& idStaggeredBoundaryVectEast,
                       const std::vector<UInt>& idStaggeredBoundaryVectNorth,
                       const std::vector<UInt>& idStaggeredBoundaryVectSouth)
{
  
  
    // +-----------------------------------------------+
    // |              Horizontal Velocity              |
    // +-----------------------------------------------+
  
  
  for ( const auto & Id : idStaggeredInternalVectHorizontal )
  {
    
    
    const UInt i     = Id / ( ncols + 1 ),
    j     = Id % ( ncols + 1 ),
    ID_NE = Id - i,        // v
    ID_NW = ID_NE - 1,     // v
    ID_SE = ID_NE + ncols, // v
    ID_SW = ID_NW + ncols; // v
    
    
    
    Vector2D vel( std::array<Real,2>{{ u[ Id ], ( v[ ID_NE ] + v[ ID_NW ] + v[ ID_SE ] + v[ ID_SW ] ) / 4. }} );
    
    
    
    
    
    const auto Dx = vel / pixel_size * dt_DSV;
    
    
    
    
    Vector2D xx( std::array<Real,2>{{ Real( j ), Real( i ) }} );  // Top-Left reference frame
    
    
    
    xx = xx - Dx;
    
    
    
    
    auto x = xx( 0 ),
    y = xx( 1 );
    
    
    auto x_1 = std::floor( x ),              //  11    21     ----> x
    y_1 = std::floor( y ),              //
    x_2 = x_1 + 1,                      //
    y_2 = y_1 + 1;                      //  12    22
      //
      // |
      // | y
    
    
    
    
    if ( x_1 == ncols )
    {
      x_2 -= 1;
      x_1 -= 1;
      x   -= 1;
    }
    if ( x_2 == 0 )
    {
      x_2 += 1;
      x_1 += 1;
      x   += 1;
    }
    
    
    if ( y_1 == ( nrows - 1 ) )
    {
      y_2 -= 1;
      y_1 -= 1;
      y   -= 1;
    }
    if ( y_2 == 0 )
    {
      y_2 += 1;
      y_1 += 1;
      y   += 1;
    }
    
    
    
      // return 0 for target values that are out of bounds
    if ( x_2 < 0 || x_1 > ncols || y_2 < 0 || y_1 > ( nrows - 1 ) )
    {
      u_star[ Id ] = 0.;
    }
    else
    {
      
      const auto Id_11 = x_1 + y_1 * ( ncols + 1 ), // u
      Id_12 = x_1 + y_2 * ( ncols + 1 ), // u
      Id_21 = x_2 + y_1 * ( ncols + 1 ), // u
      Id_22 = x_2 + y_2 * ( ncols + 1 ); // u
      
      
      
      
        // compute weights
      const Real w_x2 = x_2 - x,
      w_x1 = x   - x_1,
      w_y2 = y_2 - y,
      w_y1 = y   - y_1;
      
      
      
      const auto a = u[ Id_11 ] * w_x2 + u[ Id_21 ] * w_x1,
      b = u[ Id_12 ] * w_x2 + u[ Id_22 ] * w_x1;
      
      
      u_star[ Id ] = ((std::abs(slope_x[ Id ]) < slope_thr) ? (a * w_y2 + b * w_y1) : 0);
      
      
    }
    
  }
  
  
  
  
  for ( const auto & Id : idStaggeredBoundaryVectWest )
  {
    
    const auto Idd = Id + 1;
    
    u_star[ Id ] = ((std::abs(slope_x[ Id ]) < slope_thr) ? (u_star[ Idd ] * ( u[ Idd ] < 0. )) : 0);
    
  }
  
  
  for ( const auto & Id : idStaggeredBoundaryVectEast )
  {
    
    const auto Idd = Id - 1;
    
    u_star[ Id ] = ((std::abs(slope_x[ Id ]) < slope_thr) ? (u_star[ Idd ] * ( u[ Idd ] > 0. )) : 0);
    
  }
  
  
  
    // +-----------------------------------------------+
    // |              Vertical Velocity                |
    // +-----------------------------------------------+
  
  
  for ( const auto & Id : idStaggeredInternalVectVertical )
  {
    
    
    const auto i     = Id / ncols,
    j     = Id % ncols,
    ID_SW = Id + i,                // u
    ID_SE = ID_SW + 1,             // u
    ID_NW = ID_SW - ( ncols + 1 ), // u
    ID_NE = ID_NW + 1;             // u
    
    
    
    
    Vector2D vel( std::array<Real,2>{{ ( u[ ID_SW ] + u[ ID_SE ] + u[ ID_NW ] + u[ ID_NE ] ) / 4., v[ Id ] }} );
    
    
    
    const auto Dx = vel / pixel_size * dt_DSV;
    
    
    
    
    
    Vector2D xx( std::array<Real,2>{{ Real( j ), Real( i ) }} );  // Top-Left reference frame
    
    
    
    xx = xx - Dx;
    
    
    
    
    auto x = xx( 0 ),
    y = xx( 1 );
    
    
    auto x_1 = std::floor( x ),              //  11    21     ----> x
    y_1 = std::floor( y ),              //
    x_2 = x_1 + 1,                      //
    y_2 = y_1 + 1;                      //  12    22
      //
      // |
      // | y
    
    
    
    
    if ( x_1 == ( ncols - 1 ) )
    {
      x_2 -= 1;
      x_1 -= 1;
      x   -= 1;
    }
    if ( x_2 == 0 )
    {
      x_2 += 1;
      x_1 += 1;
      x   += 1;
    }
    
    
    if ( y_1 == nrows )
    {
      y_2 -= 1;
      y_1 -= 1;
      y   -= 1;
    }
    if ( y_2 == 0 )
    {
      y_2 += 1;
      x_2 += 1;
      y   += 1;
    }
    
    
    
      // return 0 for target values that are out of bounds
    if ( x_2 < 0 || x_1 > ( ncols - 1 ) || y_2 < 0 || y_1 > nrows )
    {
      v_star[ Id ] = 0.;
    }
    else
    {
      
      
      const auto Id_11 = x_1 + y_1 * ncols,
      Id_12 = x_1 + y_2 * ncols,
      Id_21 = x_2 + y_1 * ncols,
      Id_22 = x_2 + y_2 * ncols;
      
      
      
        // compute weights
      const Real w_x2 = x_2 - x,
      w_x1 = x   - x_1,
      w_y2 = y_2 - y,
      w_y1 = y   - y_1;
      
      
      const auto a = v[ Id_11 ] * w_x2 + v[ Id_21 ] * w_x1,
      b = v[ Id_12 ] * w_x2 + v[ Id_22 ] * w_x1;
      
      
      v_star[ Id ] = ((std::abs(slope_y[ Id ]) < slope_thr) ? (a * w_y2 + b * w_y1) : 0);
      
    }
    
    
    
    
  }
  
  
  
  
  for ( const auto & Id : idStaggeredBoundaryVectNorth )
  {
    
    const auto Idd = Id + ncols;
    
    v_star[ Id ] = ((std::abs(slope_y[ Id ]) < slope_thr) ? (v_star[ Idd ] * ( v[ Idd ] < 0. )) : 0);
    
  }
  
  for ( const auto & Id : idStaggeredBoundaryVectSouth )
  {
    
    const auto Idd = Id - ncols;
    
    v_star[ Id ] = ((std::abs(slope_y[ Id ]) < slope_thr) ? (v_star[ Idd ] * ( v[ Idd ] > 0. )) : 0);
    
  }
  
  
}




int
computePourCell(const int& IDcell,
                const UInt& N_cols,
                const std::vector<Real>& oro,
                const std::set<UInt>& idBasinVect,
                const std::set<UInt>& idStaggeredBoundaryVectSouth,
                const std::set<UInt>& idStaggeredBoundaryVectNorth,
                const std::set<UInt>& idStaggeredBoundaryVectWest,
                const std::set<UInt>& idStaggeredBoundaryVectEast)
{
  
  int candidate_id = -1;
    
  // +-----------------------------------------------+
  // |               "Cartesian" cells               |
  // +-----------------------------------------------+
  
  const int i = IDcell / N_cols;

  const int IDsouth = IDcell + N_cols,
             IDnorth = IDcell,
             IDwest  = IDcell + i,
             IDeast  = IDcell + i + 1;

  const auto iterator_south = idStaggeredBoundaryVectSouth.find(IDsouth),
             iterator_north = idStaggeredBoundaryVectNorth.find(IDnorth),
             iterator_west  = idStaggeredBoundaryVectWest .find(IDwest),
             iterator_east  = idStaggeredBoundaryVectEast .find(IDeast);

  std::vector<UInt> candidates;
  candidates.reserve(8);

  bool is_north = false,
       is_south = false,
       is_west  = false,
       is_east  = false;

  if (iterator_south == idStaggeredBoundaryVectSouth.end())
  {
    const int IDcell_south = IDcell + N_cols;
    candidates.push_back( IDcell_south );

    is_south = true;
  }

  if (iterator_north == idStaggeredBoundaryVectNorth.end())
  {
    const int IDcell_north = IDcell - N_cols;
    candidates.push_back( IDcell_north );

    is_north = true;
  }

  if (iterator_west == idStaggeredBoundaryVectWest.end())
  {
    const int IDcell_west = IDcell - 1;
    candidates.push_back( IDcell_west );

    is_west = true;
  }

  if (iterator_east == idStaggeredBoundaryVectEast.end())
  {
    const int IDcell_east = IDcell + 1;
    candidates.push_back( IDcell_east );

    is_east = true;
  }

  // +-----------------------------------------------+
  // |               Diagonal cells                  |
  // +-----------------------------------------------+
  
  
  if (is_north && is_west)
  {
    // exist north-west cell for sure
    const int IDcel_nw = IDcell - 1 - N_cols;

    if (idBasinVect.find(IDcel_nw) != idBasinVect.end())
    {
      candidates.push_back( IDcel_nw );
    }
  }

  if (is_north && is_east)
  {
    // exist north-east cell for sure
    const int IDcel_ne = IDcell + 1 - N_cols;
    
    if (idBasinVect.find(IDcel_ne) != idBasinVect.end())
    {
      candidates.push_back( IDcel_ne );
    }
  }
  
  if (is_south && is_west)
  {
    // exist south-west cell for sure
    const int IDcel_sw = IDcell - 1 + N_cols;
    
    if (idBasinVect.find(IDcel_sw) != idBasinVect.end())
    {
      candidates.push_back( IDcel_sw );
    }
  }
  
  if (is_south && is_east)
  {
    // exist south-east cell for sure
    const int IDcel_se = IDcell + 1 + N_cols;
    
    if (idBasinVect.find(IDcel_se) != idBasinVect.end())
    {
      candidates.push_back( IDcel_se );
    }
  }

  // +-----------------------------------------------+
  // |       Compute min oro cell (pour point)       |
  // +-----------------------------------------------+

  Real current_minimum_oro = oro[IDcell];
  for (const auto & it : candidates)
  {
    const auto & candidate_minimum_oro = oro[it];
    if ( candidate_minimum_oro < current_minimum_oro )
    {
      current_minimum_oro = candidate_minimum_oro;
      candidate_id = it;
    }
  }


  return candidate_id;
}





void
computeAdjacencies (const std::vector<Real>& basin_mask_Vec,

                    std::vector<UInt>& idStaggeredBoundaryVectSouth,
                    std::vector<UInt>& idStaggeredBoundaryVectNorth,
                    std::vector<UInt>& idStaggeredBoundaryVectWest,
                    std::vector<UInt>& idStaggeredBoundaryVectEast,

                    std::vector<UInt>& idStaggeredInternalVectHorizontal,
                    std::vector<UInt>& idStaggeredInternalVectVertical,

                    std::vector<UInt>& idBasinVect,
                    std::vector<UInt>& idBasinVectReIndex,

                    const UInt&              N_rows,
                    const UInt&              N_cols)
{

  // +-----------------------------------------------+
  // |                 Basin H IDs                   |
  // +-----------------------------------------------+

  UInt h = 0;
  for ( UInt i = 0; i < N_rows; i++ )
    {
      for ( UInt j = 0; j < N_cols; j++ )
        {
          const UInt k = j + i * N_cols;

          idBasinVectReIndex.push_back( h );

          if ( basin_mask_Vec[ k ] == 1 )
            {
              idBasinVect.push_back( k );
              //                std::cout << k << std::endl;
              h++;
            }
        }
    }

  // +-----------------------------------------------+
  // |         Vertical Vel. Staggered IDs           |
  // +-----------------------------------------------+

  // cycle on centered cells
  for ( UInt i = 0; i < N_rows; i++ )
    {
      for ( UInt j = 0; j < N_cols; j++ )
        {
          const UInt IDcell       = j + i * N_cols,
            IDcell_south = IDcell + N_cols,

            IDvel        = IDcell + N_cols; // interface between IDcell and IDcell_south


          if ( i != ( N_rows - 1 ) )
            {
              if ( ( basin_mask_Vec[ IDcell ] + basin_mask_Vec[ IDcell_south ] ) == 1 ) // interface cell
                {
                  if ( basin_mask_Vec[ IDcell ] == 0 )
                    {
                      idStaggeredBoundaryVectNorth.push_back( IDvel );
                    }
                  else
                    {
                      idStaggeredBoundaryVectSouth.push_back( IDvel );
                    }
                }

              if ( i == 0 && basin_mask_Vec[ IDcell ] == 1 )
                {
                  const auto IDvel_north = IDcell;

                  idStaggeredBoundaryVectNorth.push_back( IDvel_north ); // it is ok
                }

              if ( ( basin_mask_Vec[ IDcell ] + basin_mask_Vec[ IDcell_south ] ) == 2 )
                {
                  idStaggeredInternalVectVertical.push_back( IDvel ); // it is ok no repetition
                }

            }
          else
            {

              if ( basin_mask_Vec[ IDcell ] == 1 )
                {
                  const auto IDvel_south = IDvel;

                  idStaggeredBoundaryVectSouth.push_back( IDvel_south );
                }



            }


        }
    }

  // +-----------------------------------------------+
  // |         Horizontal Vel. Staggered IDs         |
  // +-----------------------------------------------+


  // cicle on centered cells
  for ( UInt i = 0; i < N_rows; i++ )
    {
      for ( UInt j = 0; j < N_cols; j++ )
        {
          const UInt IDcell      = j + i * N_cols,
            IDcell_east = IDcell + 1,

            IDvel       = IDcell + i + 1; // interface between IDcell and IDcell_east


          if ( j != ( N_cols - 1 ) )
            {
              if ( ( basin_mask_Vec[ IDcell ] + basin_mask_Vec[ IDcell_east ] ) == 1 )
                {
                  if ( basin_mask_Vec[ IDcell ] == 0 )
                    {
                      idStaggeredBoundaryVectWest.push_back( IDvel );
                    }
                  else
                    {
                      idStaggeredBoundaryVectEast.push_back( IDvel );
                    }
                }

              if ( j == 0 && basin_mask_Vec[ IDcell ] == 1 )
                {
                  const auto IDvel_west = IDcell + i;

                  idStaggeredBoundaryVectWest.push_back( IDvel_west );
                }

              if ( ( basin_mask_Vec[ IDcell ] + basin_mask_Vec[ IDcell_east ] ) == 2 )
                {
                  idStaggeredInternalVectHorizontal.push_back( IDvel );
                }

            }
          else
            {

              if ( basin_mask_Vec[ IDcell ] == 1 )
                {
                  const auto IDvel_east = IDvel;

                  idStaggeredBoundaryVectEast.push_back( IDvel_east );
                }



            }


        }
    }

  /*
    for (UInt i = 0; i < idStaggeredBoundaryVectEast.size(); i++)
    {
    //        if ( idBasinVect[ i ] != idBasinVectReIndex[ i ])
    std::cout << idStaggeredBoundaryVectEast[ i ] << std::endl;
    }
    exit(1);*/


}


void
computeAdjacencies (const std::vector<Real>& basin_mask_Vec_input,
                    const std::vector<std::tuple<bool, int> >& excluded_ids,
                    
                    std::vector<UInt>& idStaggeredBoundaryVectSouth,
                    std::vector<UInt>& idStaggeredBoundaryVectNorth,
                    std::vector<UInt>& idStaggeredBoundaryVectWest,
                    std::vector<UInt>& idStaggeredBoundaryVectEast,
                    
                    std::vector<UInt>& idStaggeredInternalVectHorizontal,
                    std::vector<UInt>& idStaggeredInternalVectVertical,
                    
                    std::vector<UInt>& idBasinVect,
                    std::vector<UInt>& idBasinVectReIndex,
                    
                    const UInt&              N_rows,
                    const UInt&              N_cols)
{
  
  
  std::vector<Real> basin_mask_Vec = basin_mask_Vec_input;
  
  for ( UInt i = 0; i < basin_mask_Vec.size(); i++ )
  {
    if ( std::get<0>( excluded_ids[ i ] ) )
    {
      basin_mask_Vec[ i ] = 0.;
    }
  }
  
  
    // +-----------------------------------------------+
    // |                 Basin H IDs                   |
    // +-----------------------------------------------+
  
  UInt h = 0;
  for ( UInt i = 0; i < N_rows; i++ )
  {
    for ( UInt j = 0; j < N_cols; j++ )
    {
      const UInt k = j + i * N_cols;
      
      idBasinVectReIndex.push_back( h );
      
      if ( basin_mask_Vec[ k ] == 1 )
      {
        idBasinVect.push_back( k );
          //                std::cout << k << std::endl;
        h++;
      }
    }
  }
  
    // +-----------------------------------------------+
    // |         Vertical Vel. Staggered IDs           |
    // +-----------------------------------------------+
  
    // cycle on centered cells
  for ( UInt i = 0; i < N_rows; i++ )
  {
    for ( UInt j = 0; j < N_cols; j++ )
    {
      const UInt IDcell       = j + i * N_cols,
      IDcell_south = IDcell + N_cols,
      
      IDvel        = IDcell + N_cols; // interface between IDcell and IDcell_south
      
      
      if ( i != ( N_rows - 1 ) )
      {
        if ( ( basin_mask_Vec[ IDcell ] + basin_mask_Vec[ IDcell_south ] ) == 1 ) // interface cell
        {
          if ( basin_mask_Vec[ IDcell ] == 0 )
          {
            idStaggeredBoundaryVectNorth.push_back( IDvel );
          }
          else
          {
            idStaggeredBoundaryVectSouth.push_back( IDvel );
          }
        }
        
        if ( i == 0 && basin_mask_Vec[ IDcell ] == 1 )
        {
          const auto IDvel_north = IDcell;
          
          idStaggeredBoundaryVectNorth.push_back( IDvel_north ); // it is ok
        }
        
        if ( ( basin_mask_Vec[ IDcell ] + basin_mask_Vec[ IDcell_south ] ) == 2 )
        {
          idStaggeredInternalVectVertical.push_back( IDvel ); // it is ok no repetition
        }
        
      }
      else
      {
        
        if ( basin_mask_Vec[ IDcell ] == 1 )
        {
          const auto IDvel_south = IDvel;
          
          idStaggeredBoundaryVectSouth.push_back( IDvel_south );
        }
        
        
        
      }
      
      
    }
  }
  
    // +-----------------------------------------------+
    // |         Horizontal Vel. Staggered IDs         |
    // +-----------------------------------------------+
  
  
    // cicle on centered cells
  for ( UInt i = 0; i < N_rows; i++ )
  {
    for ( UInt j = 0; j < N_cols; j++ )
    {
      const UInt IDcell      = j + i * N_cols,
      IDcell_east = IDcell + 1,
      
      IDvel       = IDcell + i + 1; // interface between IDcell and IDcell_east
      
      
      if ( j != ( N_cols - 1 ) )
      {
        if ( ( basin_mask_Vec[ IDcell ] + basin_mask_Vec[ IDcell_east ] ) == 1 )
        {
          if ( basin_mask_Vec[ IDcell ] == 0 )
          {
            idStaggeredBoundaryVectWest.push_back( IDvel );
          }
          else
          {
            idStaggeredBoundaryVectEast.push_back( IDvel );
          }
        }
        
        if ( j == 0 && basin_mask_Vec[ IDcell ] == 1 )
        {
          const auto IDvel_west = IDcell + i;
          
          idStaggeredBoundaryVectWest.push_back( IDvel_west );
        }
        
        if ( ( basin_mask_Vec[ IDcell ] + basin_mask_Vec[ IDcell_east ] ) == 2 )
        {
          idStaggeredInternalVectHorizontal.push_back( IDvel );
        }
        
      }
      else
      {
        
        if ( basin_mask_Vec[ IDcell ] == 1 )
        {
          const auto IDvel_east = IDvel;
          
          idStaggeredBoundaryVectEast.push_back( IDvel_east );
        }
        
        
        
      }
      
      
    }
  }
  
  /*
   for (UInt i = 0; i < idStaggeredBoundaryVectEast.size(); i++)
   {
   //        if ( idBasinVect[ i ] != idBasinVectReIndex[ i ])
   std::cout << idStaggeredBoundaryVectEast[ i ] << std::endl;
   }
   exit(1);*/
  
  
}






void
buildMatrix (const std::vector<Real>& H_int_x,
             const std::vector<Real>& H_int_y,
             const std::vector<Real>& orography,
             const std::vector<Real>& u_star,
             const std::vector<Real>& v_star,
             const std::vector<Real>& u,
             const std::vector<Real>& v,
             const Eigen::VectorXd&   H,
             const UInt&              N_cols,
             const UInt&              N_rows,
             const UInt&              N,
             const Real&              c1,
             const Real&              c3,
             const Real&              H_min,
             const std::vector<Real>& precipitation,
             const Real&              dt_DSV,
             const std::vector<Real>& alfa_x,
             const std::vector<Real>& alfa_y,
             const std::vector<UInt>& idStaggeredInternalVectHorizontal,
             const std::vector<UInt>& idStaggeredInternalVectVertical,
             const std::vector<UInt>& idStaggeredBoundaryVectWest,
             const std::vector<UInt>& idStaggeredBoundaryVectEast,
             const std::vector<UInt>& idStaggeredBoundaryVectNorth,
             const std::vector<UInt>& idStaggeredBoundaryVectSouth,
             const std::vector<UInt>& idBasinVect,
             const std::vector<UInt>& idBasinVect_not_excluded,
             const std::vector<UInt>& idStaggeredInternalVectHorizontal_not_excluded,
             const std::vector<UInt>& idStaggeredInternalVectVertical_not_excluded,
             const std::vector<UInt>& idBasinVectReIndex,
             const bool&              isNonReflectingBC,
             const bool&              isH,
             
             const std::vector<std::tuple<bool, int> >& excluded_ids, // bool, value 1 in excluded ids  UInt, Id del pour point
             std::vector<Real>&                   additional_source_term,
             
             std::vector<Eigen::Triplet<Real> >& coefficients,
             Eigen::VectorXd&                    rhs)
{
  
  
  // Be careful to maintain the conservativity!
  // cycle over boundary interfaces and check if one of the left, right cells are excluded ones
  
  for ( const auto & k_ex : idBasinVect_not_excluded )
  {
    const auto & current_tuple = excluded_ids[ k_ex ];
    if ( std::get<0>( current_tuple ) )
    {
      const auto & k_pour = std::get<1>( current_tuple );
      if ( k_pour >= 0 )
      {
        additional_source_term[ k_pour ] += precipitation[ k_ex ] * dt_DSV;
      }
    }
  }
  
  
  for ( const auto & Id : idBasinVect )
  {
    const auto IDreIndex = idBasinVectReIndex[ Id ];
    
    coefficients.push_back( Eigen::Triplet<Real>( IDreIndex,  IDreIndex,  1. ) );
    rhs( IDreIndex ) = H( Id ) + precipitation[ Id ] * dt_DSV;
  }
  
  
  for ( const auto & Id : idStaggeredInternalVectHorizontal )
  {
     
    const UInt i       = Id / ( N_cols + 1 ),
    
    IDleft  = Id - i - 1, // H
    IDright = Id - i, // H
    IDleftReIndex = idBasinVectReIndex[ IDleft ], // H
    IDrightReIndex = idBasinVectReIndex[ IDright ]; // H
    
    
    
      // define H at interfaces
    const auto H_interface = H_int_x[ Id ];
    
    const Real coeff_m = H_interface * alfa_x[ Id ];
    
    
    
    if ( H_interface > H_min )
    {
      
      coefficients.push_back( Eigen::Triplet<Real>( IDleftReIndex, IDrightReIndex, - c3 * coeff_m ) );
      coefficients.push_back( Eigen::Triplet<Real>( IDleftReIndex,  IDleftReIndex,   c3 * coeff_m ) );
      
      rhs( IDleftReIndex )  += - c1 * ( + coeff_m  * u_star[ Id ] ) - ( orography[ IDleft ] - orography[ IDright ] ) * c3 * coeff_m * isH;
      
      coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDleftReIndex, - c3 * coeff_m ) );
      coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDrightReIndex,  c3 * coeff_m ) );
      
      rhs( IDrightReIndex ) += - c1 * ( - coeff_m  * u_star[ Id ] ) - ( orography[ IDright ] - orography[ IDleft ] ) * c3 * coeff_m * isH;
      
      
    }
    
  }
  
  
  
  for ( const auto & Id : idStaggeredBoundaryVectWest )
  {
    
    const UInt i            = Id / ( N_cols + 1 ),
    
    IDright      = Id - i,
    IDrightright = IDright + 1, // H
    IDrightReIndex = idBasinVectReIndex[ IDright ]; // H
    
    
    
      // define H at interfaces
    const auto H_interface = H_int_x[ Id ];
    
    
    
    const Real coeff_m = H_interface * alfa_x[ Id ];
    
    
    
    if ( H_interface > H_min )
    {
//      rhs ( IDrightReIndex ) += isNonReflectingBC * ( - c1 * ( - coeff_m  * u_star[ Id ] ) - ( orography[ IDrightright ] - orography[ IDright ] ) * c3 * coeff_m );
      rhs ( IDrightReIndex ) += isNonReflectingBC * ( - c1 * ( - coeff_m  * u_star[ Id ] ) );
    }
    
    
    
  }
  
  
  
  for ( const auto & Id : idStaggeredBoundaryVectEast )
  {
    
    const UInt i          = Id / ( N_cols + 1 ),
    
    IDleft     = Id - i - 1,
    IDleftleft = IDleft - 1, // H
    IDleftReIndex = idBasinVectReIndex[ IDleft ]; // H
    
    
    
      // define H at interfaces
    const auto H_interface = H_int_x[ Id ];
    
    
    
    const Real coeff_m = H_interface * alfa_x[ Id ];
    
    
    
    if ( H_interface > H_min )
    {
//      rhs ( IDleftReIndex ) += isNonReflectingBC * ( - c1 * ( + coeff_m  * u_star[ Id ] ) - ( orography[ IDleftleft ] - orography[ IDleft ] ) * c3 * coeff_m );
      rhs ( IDleftReIndex ) += isNonReflectingBC * ( - c1 * ( + coeff_m  * u_star[ Id ] ) );
    }
    
  }
  
  for ( const auto & Id : idStaggeredInternalVectVertical )
  {
    
    const UInt IDleft  = Id - N_cols, // H
    IDright = Id, // H
    IDleftReIndex = idBasinVectReIndex[ IDleft ], // H
    IDrightReIndex = idBasinVectReIndex[ IDright ]; // H
    
    
    
      // define H at interfaces
    const auto H_interface = H_int_y[ Id ];
    
    
    const Real coeff_m = H_interface * alfa_y[ Id ];
    
    
    
    if ( H_interface > H_min )
    {
      
      coefficients.push_back( Eigen::Triplet<Real>( IDleftReIndex, IDrightReIndex, - c3 * coeff_m ) );
      coefficients.push_back( Eigen::Triplet<Real>( IDleftReIndex,  IDleftReIndex,   c3 * coeff_m ) );
      
      rhs( IDleftReIndex )  += - c1 * ( + coeff_m  * v_star[ Id ] ) - ( orography[ IDleft ] - orography[ IDright ] ) * c3 * coeff_m * isH;
      
      coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDleftReIndex, - c3 * coeff_m ) );
      coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDrightReIndex,  c3 * coeff_m ) );
      
      rhs( IDrightReIndex ) += - c1 * ( - coeff_m  * v_star[ Id ] ) - ( orography[ IDright ] - orography[ IDleft ] ) * c3 * coeff_m * isH;
      
      
    }
    
  }
  
  for ( const auto & Id : idStaggeredBoundaryVectNorth )
  {
    
    const UInt IDright      = Id,
    IDrightright = Id + N_cols, //
    IDrightReIndex = idBasinVectReIndex[ IDright ]; // H
    
    
    
      // define H at interfaces
    const auto H_interface = H_int_y[ Id ];
    
    
    const Real coeff_m = H_interface * alfa_y[ Id ];
    

    
    if ( H_interface > H_min )
    {
//      rhs ( IDrightReIndex ) += isNonReflectingBC * ( - c1 * ( - coeff_m  * v_star[ Id ] ) - ( orography[ IDrightright ] - orography[ IDright ] ) * c3 * coeff_m );
      rhs ( IDrightReIndex ) += isNonReflectingBC * ( - c1 * ( - coeff_m  * v_star[ Id ] ) );

    }
    
  }
  
  
  for ( const auto & Id : idStaggeredBoundaryVectSouth )
  {
    
    const UInt IDleft     = Id - N_cols, // H
    IDleftleft = IDleft - N_cols, // H
    IDleftReIndex = idBasinVectReIndex[ IDleft ]; // H
    
    
      // define H at interfaces
    const auto H_interface = H_int_y[ Id ];
    
    
    const Real coeff_m = H_interface * alfa_y[ Id ];
    
    
    if ( H_interface > H_min )
    {
//      rhs ( IDleftReIndex )  += isNonReflectingBC * ( - c1 * ( + coeff_m  * v_star[ Id ] ) - ( orography[ IDleftleft ] - orography[ IDleft ] ) * c3 * coeff_m );
      rhs ( IDleftReIndex )  += isNonReflectingBC * ( - c1 * ( + coeff_m  * v_star[ Id ] ) );
    }
    
  }
  
  
  
  
  
  for ( const auto & Id : idStaggeredInternalVectHorizontal_not_excluded )
  {
    
    const UInt i       = Id / ( N_cols + 1 ),
    
    IDleft  = Id - i - 1, // H
    IDright = Id - i; // H
    
    
    
      // define H at interfaces
    const auto H_interface = H_int_x[ Id ];
    
  
    
    if ( H_interface > H_min )
    {
      if ( std::get<0>( excluded_ids[ IDleft ] ) )
      {
        const auto & k_pour = std::get<1>( excluded_ids[ IDleft ] );
        if ( k_pour >= 0 )
        {
          additional_source_term[ k_pour ] += H_interface * std::abs(u[ Id ]) * c1;
        }
      }
      
      if ( std::get<0>( excluded_ids[ IDright ] ) )
      {
        const auto & k_pour = std::get<1>( excluded_ids[ IDright ] );
        if ( k_pour >= 0 )
        {
          additional_source_term[ k_pour ] += H_interface * std::abs(u[ Id ]) * c1;
        }
      }
      
    }
    
  }
  
  
  for ( const auto & Id : idStaggeredInternalVectVertical_not_excluded )
  {
    
    const UInt IDleft  = Id - N_cols, // H
    IDright = Id; // H
    
    
    
      // define H at interfaces
    const auto H_interface = H_int_y[ Id ];
    
    
    if ( H_interface > H_min )
    {
      if ( std::get<0>( excluded_ids[ IDleft ] ) )
      {
        const auto & k_pour = std::get<1>( excluded_ids[ IDleft ] );
        if ( k_pour >= 0 )
        {
          additional_source_term[ k_pour ] += H_interface * std::abs(v[ Id ]) * c1;
        }
      }
      
      if ( std::get<0>( excluded_ids[ IDright ] ) )
      {
        const auto & k_pour = std::get<1>( excluded_ids[ IDright ] );
        if ( k_pour >= 0 )
        {
          additional_source_term[ k_pour ] += H_interface * std::abs(v[ Id ]) * c1;
        }
      }
      
      
    }
    
  }
  
  
  
  for ( const auto & Id : idBasinVect )
  {
    const auto IDreIndex = idBasinVectReIndex[ Id ];
    
    if (!std::get<0>(excluded_ids[Id]))
    {
      rhs( IDreIndex ) += additional_source_term[ Id ];
    }
  }
  
  
}






void
putDry_excludedNodes( const std::vector<UInt>& idStaggeredInternalVectHorizontal,
                      const std::vector<UInt>& idStaggeredInternalVectVertical,
                      const std::vector<UInt>& idStaggeredBoundaryVectWest,
                      const std::vector<UInt>& idStaggeredBoundaryVectEast,
                      const std::vector<UInt>& idStaggeredBoundaryVectNorth,
                      const std::vector<UInt>& idStaggeredBoundaryVectSouth,
                      const std::vector<UInt>& idBasinVect,
                      const UInt& N_cols,
                      const std::vector<std::tuple<bool, int> >& excluded_ids,

                            Eigen::VectorXd&   H,
                            Eigen::VectorXd&   eta,
                      const std::vector<Real>& orography,
                            std::vector<Real>& u,
                            std::vector<Real>& v )
{
  /*
   The static mass balance is needed to make the first set of cells at the boundary of the exlcuded region as an accumulation zone that stores the water coming from the dynamic part near cells
   */

  for ( const auto & Id : idBasinVect )
  {
    const auto & current_tuple = excluded_ids[ Id ];
    if ( std::get<0>( current_tuple ) )
    {
      H( Id ) = 0.;
      eta ( Id ) = orography[ Id ];
      
      
//      const auto& k_pour = std::get<1> ( current_tuple );
//      if (k_pour >= 0)
//      {
//        H_basin ( k_pour ) += H ( IDreIndex ); H( k_pour ) += H ( Id );
//        H_basin ( IDreIndex ) = 0.;
//      }
    }
  }
  

  for ( const auto & Id : idStaggeredInternalVectHorizontal )
  {

    const UInt i       = Id / ( N_cols + 1 ),

    IDleft  = Id - i - 1, // H
    IDright = Id - i; // H


    if ( std::get<0> (excluded_ids[IDleft]) && std::get<0> (excluded_ids[IDright]) ) // se entrambe le celle sono escluse
    {
      H( IDleft )  = 0.;
      H( IDright ) = 0.;
      u[ Id ]      = 0.;
    }

  }


  for ( const auto & Id : idStaggeredBoundaryVectWest )
  {

    const UInt i            = Id / ( N_cols + 1 ),
               IDright      = Id - i;

    if ( std::get<0> (excluded_ids[IDright]) )
    {
      H( IDright ) = 0.;
      u[ Id ]      = 0.;
    }

  }



  for ( const auto & Id : idStaggeredBoundaryVectEast )
  {

    const UInt i          = Id / ( N_cols + 1 ),
               IDleft     = Id - i - 1;


    if ( std::get<0> (excluded_ids[IDleft]) )
    {
      H( IDleft ) = 0.;
      u[ Id ]     = 0.;
    }

  }



  for ( const auto & Id : idStaggeredInternalVectVertical )
  {

    const UInt IDleft  = Id - N_cols, // H
               IDright = Id; // H


    if ( std::get<0> (excluded_ids[IDleft]) && std::get<0> (excluded_ids[IDright]) ) // se entrambe le celle sono escluse
    {
      H( IDleft )  = 0.;
      H( IDright ) = 0.;
      v[ Id ]      = 0.;
    }

  }



  for ( const auto & Id : idStaggeredBoundaryVectNorth )
  {

    const UInt IDright = Id;

    if ( std::get<0> (excluded_ids[IDright]) ) // se entrambe le celle sono escluse
    {
      H( IDright ) = 0.;
      v[ Id ]      = 0.;
    }

  }


  for ( const auto & Id : idStaggeredBoundaryVectSouth )
  {

    const UInt IDleft = Id - N_cols;

    if ( std::get<0> (excluded_ids[IDleft]) ) // se entrambe le celle sono escluse
    {
      H( IDleft ) = 0.;
      v[ Id ]     = 0.;
    }

  }


}



void
updateVel (std::vector<Real>& u,
           std::vector<Real>& v,
           const std::vector<Real>& u_star,
           const std::vector<Real>& v_star,
           const std::vector<Real>& alfa_x,
           const std::vector<Real>& alfa_y,
           const Real&              N_rows,
           const Real&              N_cols,
           const Real&              c2,
           const Real&              H_min,
           const Eigen::VectorXd&   eta,
           const Eigen::VectorXd&   H,
           const std::vector<Real>& orography,
           const std::vector<UInt>& idStaggeredInternalVectHorizontal,
           const std::vector<UInt>& idStaggeredInternalVectVertical,
           const std::vector<UInt>& idStaggeredBoundaryVectWest,
           const std::vector<UInt>& idStaggeredBoundaryVectEast,
           const std::vector<UInt>& idStaggeredBoundaryVectNorth,
           const std::vector<UInt>& idStaggeredBoundaryVectSouth,
           const bool&              isNonReflectingBC)
{


  // +-----------------------------------------------+
  // |              Update Vertical Velocity         |
  // +-----------------------------------------------+

  for ( const auto & Id : idStaggeredInternalVectVertical )
    {

      const UInt IDsouth = Id,
        IDnorth = Id - N_cols;



      const auto & H_interface = (H( IDsouth ) + H( IDnorth ))*.5 + signum(-eta(IDsouth)+eta(IDnorth)) * (-H( IDsouth ) + H( IDnorth ))*.5;
      if ( H_interface > H_min )
        {
          v[ Id ] = alfa_y[ Id ] * ( v_star[ Id ] - c2 * ( eta( IDsouth ) - eta( IDnorth ) ) );
        }
      else
        {
          v[ Id ] = 0.;
        }


    }



  // first row
  for ( const auto & Id : idStaggeredBoundaryVectNorth )
    {

      const UInt IDsouth = Id + N_cols,
      IDnorth = Id;
      
      const auto & H_interface = (H( IDnorth ))*.5 + signum(-eta(IDsouth)+eta(IDnorth)) * (-H( IDnorth ))*.5;
      if ( H_interface > H_min )
        {
//          v[ Id ] = isNonReflectingBC * alfa_y[ Id ] * ( v_star[ Id ] - c2 * ( orography[ IDsouth ] - orography[ IDnorth ] ) );
          v[ Id ] = isNonReflectingBC * alfa_y[ Id ] * ( v_star[ Id ] );
        }
      else
        {
          v[ Id ] = 0.;
        }

    }


  // last row
  for ( const auto & Id : idStaggeredBoundaryVectSouth )
    {

      const UInt IDsouth = Id - N_cols,
      IDnorth = Id - 2*N_cols;
      
      const auto & H_interface = (H( IDsouth ))*.5 + signum(-eta(IDsouth)+eta(IDnorth)) * (H( IDsouth ))*.5;
      if ( H_interface > H_min )
        {
//          v[ Id ] = isNonReflectingBC * alfa_y[ Id ] * ( v_star[ Id ] - c2 * ( orography[ IDsouth ] - orography[ IDnorth ] ) );
          v[ Id ] = isNonReflectingBC * alfa_y[ Id ] * ( v_star[ Id ] );
        }
      else
        {
          v[ Id ] = 0.;
        }

    }


  // +-----------------------------------------------+
  // |              Update Horizontal Velocity       |
  // +-----------------------------------------------+


  for ( const auto & Id : idStaggeredInternalVectHorizontal )
    {

      const UInt i       = Id / ( N_cols + 1 ),
        IDeast  = Id - i,
        IDwest  = Id - i - 1;

      
      const auto & H_interface = (H( IDeast ) + H( IDwest ))*.5 + signum(-eta(IDeast)+eta(IDwest)) * (-H( IDeast ) + H( IDwest ))*.5;
      if ( H_interface > H_min )
        {
          u[ Id ] = alfa_x[ Id ] * ( u_star[ Id ] - c2 * ( eta( IDeast ) - eta( IDwest ) ) );
        }
      else
        {
          u[ Id ] = 0.;
        }


    }


  for ( const auto & Id : idStaggeredBoundaryVectWest )
    {

      const UInt i = Id / ( N_cols + 1 ),
      IDeast = Id - i + 1,
      IDwest = Id - i;

      const auto & H_interface = (H( IDwest ))*.5 + signum(-eta(IDeast)+eta(IDwest)) * (-H( IDwest ))*.5;
      if ( H_interface > H_min )
        {
//          u[ Id ] = isNonReflectingBC * alfa_x[ Id ] * ( u_star[ Id ] - c2 * ( orography[ IDeast ] - orography[ IDwest ] ) );
          u[ Id ] = isNonReflectingBC * alfa_x[ Id ] * ( u_star[ Id ] );
        }
      else
        {
          u[ Id ] = 0.;
        }



    }


  for ( const auto & Id : idStaggeredBoundaryVectEast )
    {

      const UInt i = Id / ( N_cols + 1 ),
      IDeast = Id - i - 1,
      IDwest = Id - i - 2;

      const auto & H_interface = (H( IDeast ))*.5 + signum(-eta(IDeast)+eta(IDwest)) * (H( IDeast ))*.5;
      if ( H_interface > H_min )
        {
//          u[ Id ] = isNonReflectingBC * alfa_x[ Id ] * ( u_star[ Id ] - c2 * ( orography[ IDeast ] - orography[ IDwest ] ) );
          u[ Id ] = isNonReflectingBC * alfa_x[ Id ] * ( u_star[ Id ] );
        }
      else
        {
          u[ Id ] = 0.;
        }

    }

}






Real
maxdt (const std::vector<Real>& u,
       const std::vector<Real>& v,
       const Real&              gravity,
       const Eigen::VectorXd&   H,
       const Real&              pixel_size)
{
  
    // +-----------------------------------------------+
    // |      Estimate vertical max Courant number     |
    // +-----------------------------------------------+
  
  const Real vel_max_y = std::max( *std::max_element( v.begin(), v.end() ), std::abs( *std::min_element( v.begin(), v.end() ) ) );
  
    // +-----------------------------------------------+
    // |    Estimate horizontal max Courant number     |
    // +-----------------------------------------------+
  
  const Real vel_max_x = std::max( *std::max_element( u.begin(), u.end() ), std::abs( *std::min_element( u.begin(), u.end() ) ) );
  
  const Real Co = 0.9; // 0.3
  const Real Co_cel = 1e4; // 10
  
  const Real cel = std::sqrt( H.maxCoeff() * gravity );
  
  Real dt_candidate = Co * pixel_size / (std::max( vel_max_x, vel_max_y )+std::numeric_limits<double>::epsilon());
  dt_candidate = std::min(dt_candidate, Co_cel * pixel_size / (cel+std::numeric_limits<double>::epsilon()));
  
  
  return( dt_candidate );
  
  
}


Real
maxCourant (const std::vector<Real>& u,
            const std::vector<Real>& v,
            const Real&              c1)
{

  // +-----------------------------------------------+
  // |      Estimate vertical max Courant number     |
  // +-----------------------------------------------+

  const Real Courant_y = std::max( *std::max_element( v.begin(), v.end() ), std::abs( *std::min_element( v.begin(), v.end() ) ) );

  // +-----------------------------------------------+
  // |    Estimate horizontal max Courant number     |
  // +-----------------------------------------------+

  const Real Courant_x = std::max( *std::max_element( u.begin(), u.end() ), std::abs( *std::min_element( u.begin(), u.end() ) ) );

  return( std::max( Courant_y, Courant_x ) * c1 );


}


Real
maxCourant (const Eigen::VectorXd& H,
            const Real&            gravity,
            const Real&            c1)
{
  const Real Courant_cel = std::sqrt( H.maxCoeff() * gravity );
  return( Courant_cel * c1 );
}




Real
compute_dt_sediment (const Real&              alpha,
                     const Real&              beta,
                     const Real&              S_x,
                     const Real&              S_y,
                     const std::vector<Real>& u,
                     const std::vector<Real>& v,
                     const Real&              pixel_size,
                     const Real&              dt_DSV,
                     UInt&              numberOfSteps)
{

  const Real max_courant_number = .95;
  // +-----------------------------------------------+
  // |      Estimate vertical max Courant number     |
  // +-----------------------------------------------+

  const Real dt_y = max_courant_number * pixel_size / ( alpha * std::pow( S_y, beta ) * std::max( *std::max_element( v.begin(), v.end() ), std::abs( *std::min_element( v.begin(), v.end() ) ) ) );

  // +-----------------------------------------------+
  // |    Estimate horizontal max Courant number     |
  // +-----------------------------------------------+

  const Real dt_x = max_courant_number * pixel_size / ( alpha * std::pow( S_x, beta ) * std::max( *std::max_element( u.begin(), u.end() ), std::abs( *std::min_element( u.begin(), u.end() ) ) ) );


  Real dt_sed = std::min( dt_y, dt_x );

  //std::cout << S_x << " " << std::pow( S_x, beta ) << " " << dt_y << " " << dt_x << " " << dt_DSV << std::endl;

  dt_sed = std::min( dt_DSV / std::floor( dt_DSV / dt_sed ), dt_DSV );

  numberOfSteps = std::floor( dt_DSV / dt_sed );

  return( dt_sed );


}


void
saveVector (const Eigen::VectorXd& b,
            const std::string& Name)
{
  std::ofstream ff( Name );

  for ( UInt k = 0; k < b.size(); k++ )
    {
      ff << b( k ) << " ";
      ff << std::endl;
    }

  ff.close();

}


void
saveMatrix (const SpMat& A,
            const std::string& Name)
{
  std::ofstream ff( Name );
  for ( UInt k = 0; k < A.outerSize(); ++k )
    {
      for ( SpMat::InnerIterator it( A, k ); it; ++it)
        {
          ff << it.row()+1 << " " << it.col()+1 << " " << it.value() << std::endl;   // row index
        }
    }
  ff.close();

}


void
saveSolution (const std::string& preName,
              const std::string& flag,
              const UInt& N_rows,
              const UInt& N_cols,
              const Real& xllcorner,
              const Real& yllcorner,
              const Real& cellsize,
              const Real& NODATA_value,
              const Eigen::VectorXd& H) // it is H or orography
{

  std::ofstream ff( preName + ".asc" );

  if ( flag == "u" )
    {
      ff << "ncols ";
      ff << N_cols + 1;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows;
      ff << std::endl;
    }
  else if ( flag == "v" )
    {
      ff << "ncols ";
      ff << N_cols;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows + 1;
      ff << std::endl;
    }
  else
    {
      ff << "ncols ";
      ff << N_cols;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows;
      ff << std::endl;
    }




  ff << "xllcorner ";
  ff << xllcorner;
  ff << std::endl;


  ff << "yllcorner ";
  ff << yllcorner;
  ff << std::endl;


  ff << "cellsize ";
  ff << cellsize;
  ff << std::endl;


  ff << "NODATA_value ";
  ff << NODATA_value;
  ff << std::endl;


  if ( flag == "u" )
    {

      for ( UInt i = 0; i < N_rows; i++ )
        {

          for ( UInt j = 0; j <= N_cols; j++ )
            {

              const auto Id = j + i * ( N_cols + 1 );

              ff << H[ Id ] << " ";

            }

          ff << std::endl;

        }


    }
  else if ( flag == "v" )
    {


      for ( UInt i = 0; i <= N_rows; i++ )
        {

          for ( UInt j = 0; j < N_cols; j++ )
            {

              const auto Id = j + i * N_cols;

              ff << H[ Id ] << " ";

            }

          ff << std::endl;

        }



    }
  else // H or orography
    {

      for ( UInt i = 0; i < N_rows; i++ )
        {

          for ( UInt j = 0; j < N_cols; j++ )
            {

              const auto k = j + i * N_cols;    // H

              ff << H[ k ] << " ";

            }

          ff << std::endl;

        }

    }


  ff.close();

}



void
saveSolution (const std::string& preName,
              const UInt& N_rows,
              const UInt& N_cols,
              const Real& xllcorner,
              const Real& yllcorner,
              const Real& cellsize,
              const Real& NODATA_value,
              const std::vector<std::tuple<bool,int>> excluded_ids) // excluded regions, high slopes I hope
{
  
  {
    std::ofstream ff ( preName + "_bool" + ".asc" );
    
    ff << "ncols ";
    ff << N_cols;
    ff << std::endl;
    
    
    ff << "nrows ";
    ff << N_rows;
    ff << std::endl;
    
    
    
    ff << "xllcorner ";
    ff << xllcorner;
    ff << std::endl;
    
    
    ff << "yllcorner ";
    ff << yllcorner;
    ff << std::endl;
    
    
    ff << "cellsize ";
    ff << cellsize;
    ff << std::endl;
    
    
    ff << "NODATA_value ";
    ff << NODATA_value;
    ff << std::endl;
    
    for ( UInt i = 0; i < N_rows; i++ )
    {
      
      for ( UInt j = 0; j < N_cols; j++ )
      {
        
        const auto k = j + i * N_cols;    // H
        
        ff << std::get<0> (excluded_ids[ k ]) << " ";
        
      }
      
      ff << std::endl;
      
    }
    
    
    ff.close();
  }
  
  
  std::ofstream ff( preName + "_pour" + ".asc" );
  
  ff << "ncols ";
  ff << N_cols;
  ff << std::endl;
  
  
  ff << "nrows ";
  ff << N_rows;
  ff << std::endl;
  
  
  
  ff << "xllcorner ";
  ff << xllcorner;
  ff << std::endl;
  
  
  ff << "yllcorner ";
  ff << yllcorner;
  ff << std::endl;
  
  
  ff << "cellsize ";
  ff << cellsize;
  ff << std::endl;
  
  
  ff << "NODATA_value ";
  ff << NODATA_value;
  ff << std::endl;
  
  for ( UInt i = 0; i < N_rows; i++ )
  {
    
    for ( UInt j = 0; j < N_cols; j++ )
    {
      
      const auto k = j + i * N_cols;    // H
      
      ff << std::get<1>(excluded_ids[ k ]) << " ";
      
    }
    
    ff << std::endl;
    
  }
  
  
  ff.close();
  
}



void
saveSolution (const std::string& preName,
              const std::string& flag,
              const UInt& N_rows,
              const UInt& N_cols,
              const Real& xllcorner,
              const Real& yllcorner,
              const Real& cellsize,
              const Real& NODATA_value,
              const std::vector<Real>& H) // it is H or orography
{

  std::ofstream ff( preName + ".asc" );

  if ( flag == "u" )
    {
      ff << "ncols ";
      ff << N_cols + 1;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows;
      ff << std::endl;
    }
  else if ( flag == "v" )
    {
      ff << "ncols ";
      ff << N_cols;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows + 1;
      ff << std::endl;
    }
  else
    {
      ff << "ncols ";
      ff << N_cols;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows;
      ff << std::endl;
    }




  ff << "xllcorner ";
  ff << xllcorner;
  ff << std::endl;


  ff << "yllcorner ";
  ff << yllcorner;
  ff << std::endl;


  ff << "cellsize ";
  ff << cellsize;
  ff << std::endl;


  ff << "NODATA_value ";
  ff << NODATA_value;
  ff << std::endl;


  if ( flag == "u" )
    {

      for ( UInt i = 0; i < N_rows; i++ )
        {

          for ( UInt j = 0; j <= N_cols; j++ )
            {

              const auto Id = j + i * ( N_cols + 1 );

              ff << H[ Id ] << " ";

            }

          ff << std::endl;

        }


    }
  else if ( flag == "v" )
    {


      for ( UInt i = 0; i <= N_rows; i++ )
        {

          for ( UInt j = 0; j < N_cols; j++ )
            {

              const auto Id = j + i * N_cols;

              ff << H[ Id ] << " ";

            }

          ff << std::endl;

        }



    }
  else // H or orography
    {

      for ( UInt i = 0; i < N_rows; i++ )
        {

          for ( UInt j = 0; j < N_cols; j++ )
            {

              const auto k = j + i * N_cols;    // H

              ff << H[ k ] << " ";

            }

          ff << std::endl;

        }

    }

  ff.close();



}


void
saveSolution (const std::string& preName,
              const std::string& flag,
              const UInt& N_rows,
              const UInt& N_cols,
              const Real& xllcorner,
              const Real& yllcorner,
              const Real& cellsize,
              const Real& NODATA_value,
              const std::vector<Int>& H) // it is H or orography
{

  std::ofstream ff( preName + ".asc" );

  if ( flag == "u" )
    {
      ff << "ncols ";
      ff << N_cols + 1;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows;
      ff << std::endl;
    }
  else if ( flag == "v" )
    {
      ff << "ncols ";
      ff << N_cols;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows + 1;
      ff << std::endl;
    }
  else
    {
      ff << "ncols ";
      ff << N_cols;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows;
      ff << std::endl;
    }




  ff << "xllcorner ";
  ff << xllcorner;
  ff << std::endl;


  ff << "yllcorner ";
  ff << yllcorner;
  ff << std::endl;


  ff << "cellsize ";
  ff << cellsize;
  ff << std::endl;


  ff << "NODATA_value ";
  ff << NODATA_value;
  ff << std::endl;


  if ( flag == "u" )
    {

      for ( UInt i = 0; i < N_rows; i++ )
        {

          for ( UInt j = 0; j <= N_cols; j++ )
            {

              const auto Id = j + i * ( N_cols + 1 );

              ff << H[ Id ] << " ";

            }

          ff << std::endl;

        }


    }
  else if ( flag == "v" )
    {


      for ( UInt i = 0; i <= N_rows; i++ )
        {

          for ( UInt j = 0; j < N_cols; j++ )
            {

              const auto Id = j + i * N_cols;

              ff << H[ Id ] << " ";

            }

          ff << std::endl;

        }



    }
  else // H or orography
    {

      for ( UInt i = 0; i < N_rows; i++ )
        {

          for ( UInt j = 0; j < N_cols; j++ )
            {

              const auto k = j + i * N_cols;    // H

              ff << H[ k ] << " ";

            }

          ff << std::endl;

        }

    }


  ff.close();

}





void
saveSolution (const std::string& preName,
              const std::string& flag,
              const UInt& N_rows,
              const UInt& N_cols,
              const Real& xllcorner,
              const Real& yllcorner,
              const Real& cellsize,
              const Real& NODATA_value,
              const UInt& n,
              const std::vector<Real>& u,
              const std::vector<Real>& v,
              const Eigen::VectorXd& H) // it is H or orography
{

  std::ofstream ff( preName + std::to_string( n ) + ".asc" );


  if ( flag == "u" )
    {
      ff << "ncols ";
      ff << N_cols + 1;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows;
      ff << std::endl;
    }
  else if ( flag == "v" )
    {
      ff << "ncols ";
      ff << N_cols;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows + 1;
      ff << std::endl;
    }
  else
    {
      ff << "ncols ";
      ff << N_cols;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows;
      ff << std::endl;
    }

  ff << "xllcorner ";
  ff << xllcorner;
  ff << std::endl;


  ff << "yllcorner ";
  ff << yllcorner;
  ff << std::endl;


  ff << "cellsize ";
  ff << cellsize;
  ff << std::endl;


  ff << "NODATA_value ";
  ff << NODATA_value;
  ff << std::endl;



  if ( flag == "u" )
    {

      for ( UInt i = 0; i < N_rows; i++ )
        {

          for ( UInt j = 0; j <= N_cols; j++ )
            {

              const auto Id = j + i * ( N_cols + 1 );

              ff << u[ Id ] << " ";

            }

          ff << std::endl;

        }


    }
  else if ( flag == "v" )
    {


      for ( UInt i = 0; i <= N_rows; i++ )
        {

          for ( UInt j = 0; j < N_cols; j++ )
            {

              const auto Id = j + i * N_cols;

              ff << v[ Id ] << " ";

            }

          ff << std::endl;

        }



    }
  else // H or orography
    {

      for ( UInt i = 0; i < N_rows; i++ )
        {

          for ( UInt j = 0; j < N_cols; j++ )
            {

              const auto k = j + i * N_cols;    // H

              ff << H( k ) << " ";

            }

          ff << std::endl;

        }

    }

  ff.close();

}



void
saveSolution (const std::string& preName,
              const std::string& flag,
              const UInt& N_rows,
              const UInt& N_cols,
              const Real& xllcorner,
              const Real& yllcorner,
              const Real& cellsize,
              const Real& NODATA_value,
              const UInt& n,
              const std::vector<Real>& u,
              const std::vector<Real>& v,
              const std::vector<Real>& H) // it is H or orography
{

  std::ofstream ff( preName + std::to_string( n ) + ".asc" );


  if ( flag == "u" )
    {
      ff << "ncols ";
      ff << N_cols + 1;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows;
      ff << std::endl;
    }
  else if ( flag == "v" )
    {
      ff << "ncols ";
      ff << N_cols;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows + 1;
      ff << std::endl;
    }
  else
    {
      ff << "ncols ";
      ff << N_cols;
      ff << std::endl;


      ff << "nrows ";
      ff << N_rows;
      ff << std::endl;
    }

  ff << "xllcorner ";
  ff << xllcorner;
  ff << std::endl;


  ff << "yllcorner ";
  ff << yllcorner;
  ff << std::endl;


  ff << "cellsize ";
  ff << cellsize;
  ff << std::endl;


  ff << "NODATA_value ";
  ff << NODATA_value;
  ff << std::endl;



  if ( flag == "u" )
    {

      for ( UInt i = 0; i < N_rows; i++ )
        {

          for ( UInt j = 0; j <= N_cols; j++ )
            {

              const auto Id = j + i * ( N_cols + 1 );

              ff << u[ Id ] << " ";

            }

          ff << std::endl;

        }


    }
  else if ( flag == "v" )
    {


      for ( UInt i = 0; i <= N_rows; i++ )
        {

          for ( UInt j = 0; j < N_cols; j++ )
            {

              const auto Id = j + i * N_cols;

              ff << v[ Id ] << " ";

            }

          ff << std::endl;

        }



    }
  else // H or orography
    {

      for ( UInt i = 0; i < N_rows; i++ )
        {

          for ( UInt j = 0; j < N_cols; j++ )
            {

              const auto k = j + i * N_cols;    // H

              ff << H[ k ] << " ";

            }

          ff << std::endl;

        }

    }

  ff.close();

}




void
saveTemporalSequence (const Vector2D&    X_gauges,
                      const Real&        time,
                      const std::string& preName,
                      const Real&        H)
{
  std::ofstream ff( preName + ".txt", std::ofstream::out | std::ofstream::app );
//  if ( time == 0 )
//    {
//      ff << X_gauges( 0 ) << " " << X_gauges( 1 ) << std::endl;
//      ff << dt << std::endl;
//    }
  ff << H << " " << time << std::endl;
  ff.close();
}



void
saveTemporalSequence (const Real&        time,
                      const std::string& preName,
                      const Real&        H)
{
  std::ofstream ff( preName + ".txt", std::ofstream::out | std::ofstream::app );
  ff << H << " " << time << std::endl;
  ff.close();
}



// For gravitational layer
void
computeResiduals (const std::vector<Real>& n_x,
                  const std::vector<Real>& n_y,
                  const UInt&              N_cols,
                  const UInt&              N_rows,
                  const std::vector<Real>& h,
                  const std::vector<Real>& coeff, // conducibilitÃ  idraulica
                  const std::vector<UInt>& idStaggeredInternalVectHorizontal,
                  const std::vector<UInt>& idStaggeredInternalVectVertical,
                  const std::vector<UInt>& idStaggeredBoundaryVectWest,
                  const std::vector<UInt>& idStaggeredBoundaryVectEast,
                  const std::vector<UInt>& idStaggeredBoundaryVectNorth,
                  const std::vector<UInt>& idStaggeredBoundaryVectSouth,
                  const std::vector<UInt>& idBasinVect,
                  std::vector<Real>& h_interface_x,
                  std::vector<Real>& h_interface_y,
                  std::vector<Real>& Res_x,
                  std::vector<Real>& Res_y)
{



  //std::cout << "Courant number for gravitational layer: " << coeff * c1 << std::endl;

  // +-----------------------------------------------+
  // |                  Horizontal                   |
  // +-----------------------------------------------+



  for ( const UInt & Id : idStaggeredInternalVectHorizontal )
    {

      const UInt i      = Id / ( N_cols + 1 ), // u
        IDeast = Id - i,                  // H
        IDwest = Id - i - 1;              // H


      const Real & h_left  = h[ IDwest ],
        & h_right = h[ IDeast ];

      const Real k_c_left  = coeff[ IDwest ],
        k_c_right = coeff[ IDeast ];

      h_interface_x[ Id ] = n_x[ Id ] * ( ( k_c_left * h_left + k_c_right * h_right ) + n_x[ Id ] * ( k_c_left * h_left - k_c_right * h_right ) ) * .5;


    }


  for ( const UInt & Id : idStaggeredBoundaryVectWest )
    {

      const UInt i = Id / ( N_cols + 1 );


      const Real h_left  = 0, //0,
        h_right = h[ Id - i ];

      const Real k_c_left  = 0.,
        k_c_right = coeff[ Id - i ];

      h_interface_x[ Id ] = n_x[ Id ] * ( ( k_c_left * h_left + k_c_right * h_right ) + n_x[ Id ] * ( k_c_left * h_left - k_c_right * h_right ) ) * .5;


    }


  for ( const UInt & Id : idStaggeredBoundaryVectEast )
    {

      const UInt i = Id / ( N_cols + 1 );


      const Real h_left  = h[ Id - i - 1 ],
        h_right = 0; //0;

      const Real k_c_left  = coeff[ Id - i - 1 ],
        k_c_right = 0.;


      h_interface_x[ Id ] = n_x[ Id ] * ( ( k_c_left * h_left + k_c_right * h_right ) + n_x[ Id ] * ( k_c_left * h_left - k_c_right * h_right ) ) * .5;

    }



  for ( const UInt & Id : idBasinVect )
    {
      const UInt i = Id / N_cols;

      Res_x[ Id ] = h_interface_x[ Id + 1 + i ] - h_interface_x[ Id + i ];
    }


  // +-----------------------------------------------+
  // |                   Vertical                    |
  // +-----------------------------------------------+


  for ( const UInt & Id : idStaggeredInternalVectVertical )
    {


      const UInt IDsouth = Id,              // H
        IDnorth = Id - N_cols;     // H


      const Real  h_left  = h[ IDnorth ],
        h_right = h[ IDsouth ];

      const Real k_c_left  = coeff[ IDnorth ],
        k_c_right = coeff[ IDsouth ];

      h_interface_y[ Id ] = n_y[ Id ] * ( ( k_c_left * h_left + k_c_right * h_right ) + n_y[ Id ] * ( k_c_left * h_left - k_c_right * h_right ) ) * .5;


    }


  for ( const UInt & Id : idStaggeredBoundaryVectNorth )
    {

      const Real h_left  = 0, // 0
        h_right = h[ Id ];

      const Real k_c_left  = 0.,
        k_c_right = coeff[ Id ];

      h_interface_y[ Id ] = n_y[ Id ] * ( ( k_c_left * h_left + k_c_right * h_right ) + n_y[ Id ] * ( k_c_left * h_left - k_c_right * h_right ) ) * .5;

    }



  for ( const UInt & Id : idStaggeredBoundaryVectSouth )
    {

      const Real h_left  = h[ Id - N_cols ],
        h_right = 0; //0;

      const Real k_c_left  = coeff[ Id - N_cols ],
        k_c_right = 0.;

      h_interface_y[ Id ] = n_y[ Id ] * ( ( k_c_left * h_left + k_c_right * h_right ) + n_y[ Id ] * ( k_c_left * h_left - k_c_right * h_right ) ) * .5;

    }


  for ( const UInt & Id : idBasinVect )
    {
      Res_y[ Id ] = h_interface_y[ Id + N_cols ] - h_interface_y[ Id ];
    }

}




// For sediment transport
void
computeResidualsTruncated (const std::vector<Real>&                 u,
                           const std::vector<Real>&                 v,
                           const UInt&                              N_cols,
                           const UInt&                              N_rows,
                           const UInt&                              N,
                           const Real&                              c1,
                           const std::vector<Real>&                 S_x,
                           const std::vector<Real>&                 S_y,
                           const Real&                              alpha,
                           const Real&                              beta,
                           const Real&                              gamma,
                           const std::vector<UInt>&                 idStaggeredInternalVectHorizontal,
                           const std::vector<UInt>&                 idStaggeredInternalVectVertical,
                           const std::vector<UInt>&                 idStaggeredBoundaryVectWest,
                           const std::vector<UInt>&                 idStaggeredBoundaryVectEast,
                           const std::vector<UInt>&                 idStaggeredBoundaryVectNorth,
                           const std::vector<UInt>&                 idStaggeredBoundaryVectSouth,
                           std::vector<std::array<Real, 2> >& Gamma_x,
                           std::vector<std::array<Real, 2> >& Gamma_y )
{


  for ( const auto & Id : idStaggeredInternalVectHorizontal )
    {
      const Real coeff_right = c1 * alpha * std::pow( std::abs( S_x[ Id ] ), beta ) * u[ Id ]
        * ( .5 - .5 * signum( u[ Id ] ) ),

        coeff_left = c1 * alpha * std::pow( std::abs( S_x[ Id ] ), beta ) * u[ Id ]
        * ( .5 + .5 * signum( u[ Id ] ) );

      Gamma_x[ Id ][ 0 ] = coeff_right;
      Gamma_x[ Id ][ 1 ] = coeff_left;
    }

  for ( const auto & Id : idStaggeredBoundaryVectWest )
    {
      const Real coeff_right = c1 * alpha * std::pow( std::abs( S_x[ Id ] ), beta ) * u[ Id ]
        * ( .5 - .5 * signum( u[ Id ] ) ),

        coeff_left = c1 * alpha * std::pow( std::abs( S_x[ Id ] ), beta ) * u[ Id ]
        * ( .5 + .5 * signum( u[ Id ] ) );

      Gamma_x[ Id ][ 0 ] = coeff_right;
      Gamma_x[ Id ][ 1 ] = coeff_left;
    }



  for ( const auto & Id : idStaggeredInternalVectVertical )
    {
      const Real coeff_right = c1 * alpha * std::pow( std::abs( S_y[ Id ] ), beta ) * v[ Id ]
        * ( .5 - .5 * signum( v[ Id ] ) ),

        coeff_left = c1 * alpha * std::pow( std::abs( S_y[ Id ] ), beta ) * v[ Id ]
        * ( .5 + .5 * signum( v[ Id ] ) );


      Gamma_y[ Id ] = std::array<Real,2>{{ coeff_right, coeff_left }};
    }

  for ( const auto & Id : idStaggeredBoundaryVectNorth )
    {
      const Real coeff_right = c1 * alpha * std::pow( std::abs( S_y[ Id ] ), beta ) * v[ Id ]
        * ( .5 - .5 * signum( v[ Id ] ) ),

        coeff_left = c1 * alpha * std::pow( std::abs( S_y[ Id ] ), beta ) * v[ Id ]
        * ( .5 + .5 * signum( v[ Id ] ) );


      Gamma_y[ Id ] = std::array<Real,2>{{ coeff_right, coeff_left }};
    }

  for ( const auto & Id : idStaggeredBoundaryVectSouth )
    {
      const Real coeff_right = c1 * alpha * std::pow( std::abs( S_y[ Id ] ), beta ) * v[ Id ]
        * ( .5 - .5 * signum( v[ Id ] ) ),

        coeff_left = c1 * alpha * std::pow( std::abs( S_y[ Id ] ), beta ) * v[ Id ]
        * ( .5 + .5 * signum( v[ Id ] ) );


      Gamma_y[ Id ] = std::array<Real,2>{{ coeff_right, coeff_left }};
    }


}



std::vector<Real>
compute_d_perc( const std::vector<Real>& clay, const std::vector<Real>& sand, const Real& perc )
{
  // linear interpolation in log10 x-scale
  std::vector<Real> d_perc( clay.size() );

  for ( UInt i = 0; i < clay.size(); i++ )
    {
      auto & d_perc_cell = d_perc[ i ];

      const auto & sand_cell = sand[ i ];
      const auto & clay_cell = clay[ i ];

      const auto Y_0 = ( 1 - sand_cell ) * 100;
      const auto Y_1 = clay_cell * 100; // ( 1 - sand_cell - silt_cell ) * 100
      const auto Y_2 = 100;

      if ( perc <= Y_0 )
        {
          const Real angular_coeff = ( Y_0 - Y_1 ) / ( std::log10( 25. ) ); // 50 \mu m / 2 \mu m
          const Real DY            = perc - Y_0;

          d_perc_cell = 50.e-6 * std::pow( 10, DY / angular_coeff );
        }
      else
        {
          const Real angular_coeff = ( Y_2 - Y_0 ) / ( std::log10( 40 ) ); // 2 mm / 50 \mu m
          const Real DY            = perc - Y_0;

          d_perc_cell = 50.e-6 * std::pow( 10, DY / angular_coeff );
        }

    }

  return d_perc;
}






