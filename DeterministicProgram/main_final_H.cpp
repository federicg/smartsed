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



//! std library
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <numeric>
#include <iosfwd>
#include <string>
#include <limits> 
#include <set>
#include <list>
#include <map>
#include <unordered_map>
#include <vector>
#include <array>
#include <chrono> 
#include <thread>
#include <ctime> 
#include <functional>


//! Parse library 
#include "GetPot.hpp" 


 
//! Eigen library    
#include "Eigen/Sparse"     
#include "Eigen/SparseCholesky"

//! Python interface --> correctly the QGIS path in Makefile
//#include "Python.h"

//! for parallelize with openmp
#include <omp.h>



//! Generic real data
typedef double Real;

typedef int8_t  int8_type;
typedef int16_t int16_type;
typedef int32_t int32_type;
typedef int64_t int64_type;

typedef uint8_t  uint8_type; 
typedef uint16_t uint16_type;
typedef uint32_t uint32_type;
typedef uint64_t uint64_type;

//! Generic short integer data
typedef int16_type  SInt;

//! Generic short unsigned integer data
typedef uint16_type SUInt;

//! Generic integer data
typedef int32_type  Int;

//! generic unsigned integer (used mainly for addressing)
typedef uint32_type UInt;



//! declares a COLUMN MAJOR sparse matrix type of Reals
typedef Eigen::SparseMatrix<Real> SpMat; 




#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif



using functionType = std::function < void (       std::vector<Real>&, 
                                                  std::vector<Real>&, 
                                            const std::vector<Real>&,
                                            const std::vector<Real>&,
                                            const std::vector<Real>&,
                                            const std::vector<Real>&,
                                            const std::vector<Real>&,
                                            const std::vector<Real>&,
                                            const std::vector<Real>&,
                                            const std::vector<Real>&,
                                            const Real&             ,
                                            const Real&             ,
                                            const Real&             ,
                                            const Real&             ,
                                            const Eigen::VectorXd&  ,
                                            const std::vector<Real>&,
                                            const std::vector<UInt>&,
                                            const std::vector<UInt>&,
                                            const std::vector<UInt>&,
                                            const std::vector<UInt>& ) >;





class Vector2D
{

public:

    //! Empty constructor (all components are set to zero)
    Vector2D()
    {}

    Vector2D( std::array<Real,2> const& indices ): // Note the use of Real indices
        M_coords( indices )
    {}
    
    //! Copy constructor
    Vector2D( Vector2D const& vector )
    {
        *this = vector;
    }    
    
        
    ~Vector2D() = default;
    
    //! Operator +=
    Vector2D& operator+= ( Vector2D const& vector )
    {
        for ( UInt i = 0; i < 2; i++ )
        {
            M_coords[ i ] += vector.M_coords[ i ];
        }
        return *this;
    }

    //! Assignment operator
    Vector2D& operator= ( Vector2D const& vector )
    {
        for ( UInt i = 0; i < 2; i++ )
        {
            M_coords[ i ] = vector.M_coords[ i ];
        }
        return *this;
    }


    Vector2D operator+ ( Vector2D const& vector ) const
    {
        Vector2D tmp( *this );
        return tmp += vector;
    }


    //! Operator -=
    Vector2D& operator-= ( Vector2D const& vector )
    {
        for ( UInt i = 0; i < 2; i++ )
        {
            M_coords[ i ] -= vector.M_coords[ i ];
        }
        return *this;
    }

    //! Operator -
    Vector2D operator- ( Vector2D const& vector ) const
    {
        Vector2D tmp( *this );
        return tmp -= vector;
    }
    
    //! Operator *= (multiplication by scalar)
    Vector2D& operator*= ( Real const& factor )
    {
        for ( UInt i = 0; i < 2; i++ )
        {
            M_coords[ i ] *= factor;
        }
        return *this;
    }

    //! Operator /= (division by scalar)
    Vector2D& operator/= ( Real const& factor )
    {
        *this *= 1. / factor;
        return *this;
    }

    //! Operator / (division by scalar)
    Vector2D operator/ ( Real const& factor ) const
    {
        Vector2D tmp( *this );
        return tmp /= factor;
    }
    
            
    Vector2D operator* ( Real const& factor )
    {
        Vector2D tmp( *this );
        return tmp *= factor;
    }

    
    Real dot ( Vector2D const& vector ) const
    {
        Real scalarProduct = 0.;
        for ( UInt i = 0; i < 2; i++ )
        {
            scalarProduct += M_coords[ i ] * vector.M_coords[ i ];
        }
        return scalarProduct;
    }    

    Real norm () const
    {
        return std::sqrt ( this->dot ( *this ) );
    }

     
    //! Operator ()
    Real const& operator() ( UInt const& i ) const
    {
        return M_coords[ i ];
    }

    //! Operator ()
    Real& operator() ( UInt const& i )
    {
        return M_coords[ i ];
    }    
    

private:

    std::array<Real,2> M_coords;

};


//! Operator * (multiplication by scalar on the right)
inline Vector2D operator* ( Vector2D const& vector, Real const& factor )
{
    Vector2D tmp( vector );
    return tmp *= factor;
}

//! Operator * (multiplication by scalar on the left)
inline Vector2D operator* ( Real const& factor, Vector2D const& vector )
{
    Vector2D tmp( vector );
    return tmp *= factor;
}


std::map<Int, std::array<Real,2> > createCN_map_Gav()
{
    std::map<Int, std::array<Real,2> > CN;
    
    CN[ 111 ] = std::array<Real,2>{{ 0.15,0.2 }};
    CN[ 112 ] = std::array<Real,2>{{ 0.35,0.2 }};
    CN[ 121 ] = std::array<Real,2>{{ 0.28,0.2 }};
    CN[ 122 ] = std::array<Real,2>{{ 1,0.2 }};
    CN[ 123 ] = std::array<Real,2>{{ 1,0.2 }};
    CN[ 124 ] = std::array<Real,2>{{ 1,0.2 }};
    CN[ 131 ] = std::array<Real,2>{{ 1,1.8 }};
    CN[ 132 ] = std::array<Real,2>{{ 1,1.8 }};
    CN[ 133 ] = std::array<Real,2>{{ 1,2 }};
    CN[ 141 ] = std::array<Real,2>{{ 0.5,1.5 }};
    CN[ 142 ] = std::array<Real,2>{{ 0.5,1.5 }};
    CN[ 211 ] = std::array<Real,2>{{ 0.8,1.8 }};
    CN[ 212 ] = std::array<Real,2>{{ 0.8,1.8 }};
    CN[ 213 ] = std::array<Real,2>{{ 0.8,1.8 }};
    CN[ 221 ] = std::array<Real,2>{{ 0.8,2 }};
    CN[ 222 ] = std::array<Real,2>{{ 0.8,1.8 }};
    CN[ 223 ] = std::array<Real,2>{{ 0.8,1.8 }};
    CN[ 231 ] = std::array<Real,2>{{ 0.6,1.6 }};
    CN[ 241 ] = std::array<Real,2>{{ 0.8,1.8 }};
    CN[ 242 ] = std::array<Real,2>{{ 0.8,1.8 }};
    CN[ 243 ] = std::array<Real,2>{{ 0.8,1.8 }};
    CN[ 244 ] = std::array<Real,2>{{ 0.8,1.6 }};
    CN[ 311 ] = std::array<Real,2>{{ 0.1,1.6 }};
    CN[ 312 ] = std::array<Real,2>{{ 0.05,1.6 }};
    CN[ 313 ] = std::array<Real,2>{{ 0.1,1.6 }};
    CN[ 321 ] = std::array<Real,2>{{ 0.5,1.6 }};
    CN[ 322 ] = std::array<Real,2>{{ 0.5,1.6 }};
    CN[ 323 ] = std::array<Real,2>{{ 0.4,1.6 }};
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
    
    return CN;
    
}


std::map<std::array<Int,2>, Int > createCN_map()
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





class Raster
{

public:


    Raster( const std::string& file )    
    {
    
        std::vector<Eigen::Triplet<Real> > cc;
    
    
    
    
        std::ifstream ff( file );
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
            std::cout << "Unable to open the file, check DEM_dir in the SMARTSED_input file" << std::endl;
            exit( -1. );
        }
        
                        
        
        Coords.resize( nrows, ncols );
                
         
        Coords.setFromTriplets( cc.begin(), cc.end() );

        

        
    }
    
    ~Raster() = default;    




    UInt ncols,
         nrows;
               
    Real xllcorner,
         yllcorner,
               
         cellsize,
               
         NODATA_value;      

    SpMat Coords;  // forse mettere una matrice densa


};


inline Real signum( const Real& x )
{ 
    return ( ( x > 0 ) - ( x < 0 ) ); 
}



class Rain // Previous interpolation not Linear 
{

public:

    Rain( const std::string& infiltrationModel,
          const UInt&        N,
          const bool&        isInitialLoss,
          const Real&        perc_initialLoss )
    {
        M_isInitialLoss = isInitialLoss;
        c = perc_initialLoss;
        
        if ( infiltrationModel == "None" )
        {
            M_infiltrationModel = false;
        }
        else if ( infiltrationModel  == "SCS-CN" )
        {
            M_infiltrationModel = true;
        }
        else
        {
            std::cout << "Insert a valid infiltration model!" << std::endl;
            exit( 1. );
        }
        
        DP_total.resize( N );
        DP_cumulative.resize( N );
        DP_infiltrated.resize( N );
        IDW_weights.resize( N );
        
    }
    
    ~Rain() = default;
    
    inline void constant_precipitation( const std::string&       file,
                                        const UInt&              ndata,
                                        const bool&              is_precipitation,
                                        const Real&              time_spacing )
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
                    if ( std::abs( ((i%rr)+1)*time_spacing - hour_full ) > .1 * ((i%rr)+1)*time_spacing )
                    {
                        std::cout << hour << std::endl;
                        std::cout << "Invalid rain file" << std::endl;
                        exit( 1. );
                    }
                    
                    
                    //value=1;
//                    std::cout << str1 << " " << hour << " " << value << " " << i+1 << std::endl;
                    
                    // mm/h  --> m/sec.
                    Hyetograph[ 0 ].push_back( value * 1.e-3 / ( time_spacing * 3600 ) );
                    //ff >> value; // second pluviometer data --> not used yet
                    
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
    
    
    inline void IDW_precipitation( const std::vector<std::string>& file_vect,
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
                    
                    
//                    std::cout << Id_sensor << " " << day << " " << i+1 << " " << hour << " " << value << std::endl;
                    
                    
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
//                        std::cout << ii << std::endl;
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
                //std::cout << IDW_weights[ IDcenter ][ ii ] << " " << *IDW_weights[ IDcenter ].begin() << " " << *IDW_weights[ IDcenter ].end();
                IDW_weights[ IDcenter ][ ii ] /= sum;
                //std::cout << " " << IDW_weights[ IDcenter ][ ii ] << std::endl;
                //IDW_weights[ IDcenter ][ ii ] = 1.;
            }
            
            
        }
        


    }
      


    
    inline void computePrecipitation( const UInt&              n,
                                      const UInt&              steps_per_hour,
                                      const std::vector<Real>& S,
                                      const std::vector<Real>& melt_mask,
                                      const std::vector<Real>& h_G,
                                      const Eigen::VectorXd&   H,
                                      const UInt&              N_rows,
                                      const UInt&              N_cols,
                                      const std::vector<UInt>& idBasinVect )
    {

        // SCS-CN method and Initial and Constant Loss Model
//        std::cout << Hyetograph.size() << " " << IDW_weights.size() << std::endl;
        for ( UInt Id = 0; Id < Hyetograph.size(); Id++ )
        {
            
//            std::cout << ( n - 1 ) / ( steps_per_hour * M_time_spacing_vect[ Id ] ) << " ietogramma" << std::endl;
            
            
            const UInt i_index = std::floor( ( n - 1 ) / ( steps_per_hour * M_time_spacing_vect[ Id ] ) );
            

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
//                    std::cout << H[ IDcenter ] << " " << h_G[ IDcenter ] << " " << c << " " << S[IDcenter ] << std::endl;
//                    std::cout << ( H[ IDcenter ] + h_G[ IDcenter ] ) << " " << ( c * S[ IDcenter ] ) << std::endl;
                    potential_runoff = 0.;
                    infiltrationRate = rainfall_intensity * melt_mask[ IDcenter ];
                }
                
                
                
                DP_total      [ IDcenter ] += rainfall_intensity;
                DP_cumulative [ IDcenter ] += potential_runoff; 
                DP_infiltrated[ IDcenter ] += infiltrationRate;
                
                //std::cout << DP_total[ IDcenter ] << " " << DP_cumulative[ IDcenter ] << " " << DP_infiltrated[ IDcenter ] << " " << potential_runoff << " " << infiltrationRate << std::endl;
                
                
                
            }
        
        }
//        exit(1);
                    
    }

    
    
    
    
    
    
    
    
    
//
//    // Compute the precipitation rate at a given time step n
//    inline void computePrecipitation( const UInt&              n,
//                                      const Real&              dt_DSV,
//                                      const UInt&              steps_per_hour,
//                                      const std::vector<Real>& basin_mask_Vec,
//                                      const std::vector<Real>& S,
//                                      const std::vector<Real>& melt_mask,
//                                      const std::vector<Real>& H_int_x,
//                                      const std::vector<Real>& H_int_y,
//                                      const std::vector<Real>& u,
//                                      const std::vector<Real>& v,
//                                      const std::vector<Real>& h_G,
//                                            std::vector<Real>& rhs_h_G,
//                                      const std::vector<Real>& S_coeff,
//                                      const std::vector<Real>& ET_ET_vec,
//                                      const std::vector<Real>& Res_x,
//                                      const std::vector<Real>& Res_y,
//                                      const Real&              pixel_size,
//                                      const UInt&              N_rows,
//                                      const UInt&              N_cols,
//                                      const std::vector<UInt>& idBasinVect )
//    {
//
//        // SCS-CN method
//
//        for ( UInt Id = 0; Id < Hyetograph.size(); Id++ )
//        {
//            const UInt i_index = std::floor( ( n - 1 ) / ( steps_per_hour * M_time_spacing_vect[ Id ] ) );
//            rainfall_intensity = Hyetograph[ Id ][ i_index ];
//
//            for ( const auto & k : idBasinVect )
//            {
//                rhs_h_G[ k ] = S_coeff[ k ] - ET_ET_vec[ k ] * ( 1 - ( rainfall_intensity > 0. ) ) - ( Res_x[ k ] + Res_y[ k ] ) / pixel_size;
//            }
//
//            for ( const auto & IDcenter : idBasinVect )
//            {
//
//
//
//                const UInt i        = IDcenter / N_cols,
//                IDsouth_ = IDcenter + N_cols, // v, v_star
//                IDnorth_ = IDcenter,          // v, v_star
//                IDeast_  = IDcenter + i + 1,  // u, u_star
//                IDwest_  = IDcenter + i;      // u, u_star
//
//
//
//                // define H at interfaces
//                const auto & H_interface_south = H_int_y[ IDsouth_ ],
//
//                & H_interface_north = H_int_y[ IDnorth_ ],
//
//                & H_interface_east  = H_int_x[ IDeast_  ],
//
//                & H_interface_west  = H_int_x[ IDwest_  ];
//
//
//                const Real deltaSoilMoisture = h_G[ IDcenter ] - S[ IDcenter ];
//                const Real exfiltration = std::max( deltaSoilMoisture, 0. ); // explicit discretization
//
//                Real weight = 0.;
//                if ( S[ IDcenter ] > 0 && deltaSoilMoisture < 0. )
//                {
//                    weight = std::pow( deltaSoilMoisture / S[ IDcenter ], 2. );
//                }
//
//
//
//                //weight=.5;
//                //std::cout << deltaSoilMoisture << " " << h_G[ IDcenter ] << std::endl;
//
//                //std::cout << exfiltration << " " << S[ IDcenter ] << " " << weight << std::endl;
//
//
//                if ( weight > 1. || h_G[ IDcenter ] < 0 )
//                {
//                    std::cout << "Error in weight infiltration model\n" <<
//                    "h_G = " << h_G[ IDcenter ] << "\nweight = " << weight <<
//                    "\nmax soil moisture ret. = " << S[ IDcenter ] << std::endl;
//                    exit( -1 );
//                }
//
//
//                const auto q_H = ( H_interface_east  * u[ IDeast_  ] - H_interface_west  * u[ IDwest_  ]
//                                  + H_interface_south * v[ IDsouth_ ] - H_interface_north * v[ IDnorth_ ] ) / pixel_size;
//
//
//                const auto infiltrationRate = ( weight   * ( rainfall_intensity * melt_mask[ IDcenter ] - q_H * is_modified_SCSCN_model ) + ( 1. - weight ) * is_modified_SCSCN_model * ( - rhs_h_G[ IDcenter ] ) );
//
//                //std::cout << weight << " " << " " << rainfall_intensity << " " << melt_mask[ IDcenter ] << std::endl;
//
//                /* // in case of modified SCS-CN happens that f is negative, so exfiltration? but creates problems with numerical stability of DSV
//                 if ( infiltrationRate < 0 )
//                 {
//                 std::cout << "negative Infiltration rate!" << std::endl;
//                 exit( -1. );
//                 }
//                 */
//
//                const auto potential_runoff = std::max( rainfall_intensity * melt_mask[ IDcenter ] - infiltrationRate, 0. );
//
//
//                DP_total      [ IDcenter ] = rainfall_intensity * dt_DSV;
//
//                // if precipitation rate exceeds infltration rate runoff occurs
//                DP_cumulative [ IDcenter ] = potential_runoff * dt_DSV + exfiltration * ( 1 - is_modified_SCSCN_model );
//                DP_infiltrated[ IDcenter ] = infiltrationRate * dt_DSV - exfiltration * ( 1 - is_modified_SCSCN_model );
//
//                //std::cout << DP_total[ IDcenter ] << " " << DP_cumulative[ IDcenter ] << " " << DP_infiltrated[ IDcenter ] << " " << potential_runoff << " " << infiltrationRate << std::endl;
//
//
//
//            }
//
//        }
//
//
//
//    }

    std::vector<Real> DP_total,
                      DP_cumulative,
                      DP_infiltrated;

private:
    
    std::vector<std::vector<Real> > Hyetograph, // # station times ndata
                                    IDW_weights;
    
    std::vector<Real> M_time_spacing_vect;

    bool M_infiltrationModel, M_isInitialLoss;


    Real rainfall_intensity = 0;

    Real c;

};



class Temperature
{

public:


    Temperature( const std::string&       file,
                 const UInt&              N,
                 const UInt&              max_Days,
                 const Real&              T_crit,
                 const std::vector<Real>& orography,
                 const UInt&              ndata,
                 const UInt&              steps_per_hour,
                 const Real&              time_spacing,
                 const Real&              height_thermometer,
                 const std::string        format_temp ):
    T_crit( T_crit ),
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
                    
                    //std::cout << str1 << std::endl;
                    if ( str1.length() != 10 )
                    {
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
                    
                    std::string hour_string  ( hour.begin(),   hour.begin()+2 ),
                                minute_string( hour.begin()+3, hour.begin()+5 ),
                                second_string( hour.begin()+6, hour.begin()+8 );
                    
                    
                    const Int hour_number   = std::stoi( hour_string   ),
                              minute_number = std::stoi( minute_string ),
                              second_number = std::stoi( second_string );
                    
                    const Real hour_full    = Real(hour_number) + Real(minute_number)/60 + Real(second_number)/3600;
                    
                    
                    Real value;  // temperature data
                    ff >> value;
                    

                    const int rr = int(std::round(24./time_spacing));
                    if ( std::abs( ((i%rr)+1)*time_spacing - hour_full ) > .1 * ((i%rr)+1)*time_spacing )
                    {
                        std::cout << hour << std::endl;
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
                    
                    
                    
                    std::cout << str1 << " " << hour << " " << value << " " << i+1 << std::endl;
                    
                    
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
                    
                    //std::cout << str1 << std::endl;
                    if ( str1.length() != 10 )
                    {
                        std::cout << "Wrong Temperature file format" << std::endl;
                        exit( 1. );
                    }
                    
                    std::string year_string ( str1.begin(),   str1.begin()+4 ),
                                month_string( str1.begin()+5, str1.begin()+7 ),
                                day_string  ( str1.begin()+8, str1.end()     );
                    
                    
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
                    
                    std::string hour_string  ( hour.begin(),   hour.begin()+2 ),
                                minute_string( hour.begin()+3, hour.begin()+5 );
                    
                    
                    const Int hour_number   = std::stoi( hour_string   ),
                              minute_number = std::stoi( minute_string );
                    
                    const Real hour_full    = Real(hour_number) + Real(minute_number)/60;
                    
                    
                    Real value;  // temperature data
                    ff >> value;
                    
                    const int rr = int(std::round(24./time_spacing));
                    if ( std::abs( (i%rr)*time_spacing - hour_full ) > .1 * (i%rr)*time_spacing )
                    {
                        std::cout << std::stoi( hour ) << std::endl;
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
                    
                    
                    
//                    std::cout << str1 << " " << hour << " " << value << " " << i+1 << std::endl;
                    
                    
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
                //std::cout << Temperature_Graph[ h ] << " " << h << std::endl;


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
            
            //std::cout << T_dailyMean[ n-1 ] << " " << T_dailyMin[n-1] << " " << T_dailyMax[n-1] << std::endl;
        }
        // --------------------------------------------- //


        // Now fill J starting from J_ndata
        for ( UInt n = 1; n <= max_Days; n++ )
        {
            const UInt i = std::floor( ( n - 1 ) * ( 24./time_spacing ) );
            J[ n - 1 ] = J_ndata[ i ];
        }

        
    }
    
    ~Temperature() = default;   


    inline void computeTemperature( const UInt&              n, 
                                    const UInt&              steps_per_hour,
                                    const Real&              time_spacing,
                                    const std::vector<Real>& orography,
                                    const std::vector<UInt>& idBasinVect )
    {

        const auto i = std::floor( ( n - 1 ) / ( steps_per_hour * time_spacing ) );

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

    std::vector<Real> T_raster,
                      melt_mask,
                      T_dailyMean,
                      T_dailyMin,
                      T_dailyMax,

                      J;
    
    const Real T_crit;

private:

    std::vector<Real> Temperature_Graph;           // length: ndata

    const Real Temp_diff = -6.5e-3;

    const Real height_th;

};




class evapoTranspiration
{

public:

    evapoTranspiration( const std::string&       ET_model,
                        const UInt&              N,
                        const std::vector<Real>& orography,
                        const std::vector<Real>& J,
                        const UInt&              max_Days,
                        const Real&              phi_rad,
                        const Real&              height_thermometer ):
    height_th( height_thermometer )
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
//            std::cout << Ra[ n - 1 ] << std::endl;
        }


    }


    ~evapoTranspiration() = default;

    inline void ET( const std::vector<Real>& T_mean, // lungo nstep: vettore delle temperature in °C
                    const std::vector<Real>& T_min,  // lungo nstep
                    const std::vector<Real>& T_max,  // lungo nstep
                    const Int&               n,       // time step number
                    const std::vector<UInt>& idBasinVect,
                    const std::vector<Real>& orography,
                    const Real&              steps_per_hour
                  )
    {

        switch ( M_evapoTranspiration_model )
        {
            case 0:

                break;

            case 1:

                const Int i = std::floor( ( n - 1 ) / ( steps_per_hour * 24 ) );
//                std::cout << "Evvv" << i << std::endl;
//                exit(1);
                for ( const auto & k : idBasinVect )
                {

                    const auto t_mean = T_mean[ i ] + Temp_diff * ( orography[ k ] - height_th ),
                               t_max  = T_max [ i ] + Temp_diff * ( orography[ k ] - height_th ),
                               t_min  = T_min [ i ] + Temp_diff * ( orography[ k ] - height_th );

                    // unity: mm/day --> m/sec.
                    ET_vec[ k ] = .0023 * Ra[ i ] * ( t_mean + 17.8 ) * std::pow( ( t_max - t_min ), .5 ) * ( 1.e-3/(24*3600) ); // più che altro mettere la T del giorno media per poter calcolare T_min e T_max
                    
//                    if ( std::isnan(ET_vec[k]) )
//                    {
//                        std::cout << "nan in evapotranspiration" << std::endl;
//                        exit( -1 );
//                    }

                } 

                //exit(1);

                break;

        }

    }

    std::vector<Real> ET_vec;

private:

    std::vector<Real> Ra;

    UInt M_evapoTranspiration_model;

    const Real M_Gsc = .082; // Solar constant

    const Real height_th;
    const Real Temp_diff = -6.5e-3;


};




class frictionClass
{

public:

    frictionClass( const std::string& friction_model, 
                   const Real& n_manning, 
                   const Real& dt_DSV,
                   const std::vector<Real>& d_90,
                   const Real& r,
                   const Real& H_min, 
                   const UInt& N_rows,
                   const UInt& N_cols,
                   const std::vector<Real>& S_x,
                   const std::vector<Real>& S_y ):
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
                    
                    const auto d_90_cell = r * d_90[ IDcell ];
                    
                    M_fc0_greater_x[ IDleft ] = std::pow( d_90_cell, .45 ) / ( .56  * std::pow( M_g, .44 ) );
                    M_fc0_lower_x  [ IDleft ] = std::pow( d_90_cell, .24 ) / ( 2.73 * std::pow( M_g, .49 ) );
                    M_fc0_greater_y[ IDup ] = std::pow( d_90_cell, .45 ) / ( .56  * std::pow( M_g, .44 ) );
                    M_fc0_lower_y  [ IDup ] = std::pow( d_90_cell, .24 ) / ( 2.73 * std::pow( M_g, .49 ) );
                    
                    M_fc0_greater_x[ IDright ] = std::pow( d_90_cell, .45 ) / ( .56  * std::pow( M_g, .44 ) );
                    M_fc0_lower_x  [ IDright ] = std::pow( d_90_cell, .24 ) / ( 2.73 * std::pow( M_g, .49 ) );
                    M_fc0_greater_y[ IDdown ] = std::pow( d_90_cell, .45 ) / ( .56  * std::pow( M_g, .44 ) );
                    M_fc0_lower_y  [ IDdown ] = std::pow( d_90_cell, .24 ) / ( 2.73 * std::pow( M_g, .49 ) );
                }
            }
            
            
            

            alfa_x.resize( S_x.size() );
            alfa_y.resize( S_y.size() );

            if ( friction_model == "None" )
            {
                M_frictionModel = 0;
            }
            else if ( friction_model == "Manning" )
            {
                M_frictionModel = 1;
                

                M_gamma_dt_DSV_x = M_dt_DSV * M_g * std::pow( M_n_manning, 2 );
                M_gamma_dt_DSV_y = M_dt_DSV * M_g * std::pow( M_n_manning, 2 );
                
            }
            else if ( friction_model == "Rickenmann" )
            {
                M_frictionModel = 2;
                
                M_gamma_dt_DSV_x = M_dt_DSV * M_g * std::pow( M_n_manning, 2 );
                M_gamma_dt_DSV_y = M_dt_DSV * M_g * std::pow( M_n_manning, 2 );

                
                
                M_expo_r_x_vect.resize( S_x.size() );
                M_gamma_dt_DSV_x_.resize  ( S_x.size() );
                for ( UInt k = 0; k < S_x.size(); k++ )
                {
                    M_expo_r_x_vect[ k ] = M_expo_r1 * ( std::abs( S_x[ k ] ) >   .006 ) +
                                           M_expo_r2 * ( std::abs( S_x[ k ] ) <=  .006 );
                    
                    
                    const auto M_Rick_x = ( M_fc0_greater_x[ k ] * std::pow( std::abs( S_x[ k ] ), .33 ) ) * ( std::abs( S_x[ k ] ) >  .006 ) +
                                          ( M_fc0_lower_x[ k ]   * std::pow( std::abs( S_x[ k ] ), .08 ) ) * ( std::abs( S_x[ k ] ) <= .006 );
                    M_gamma_dt_DSV_x_[ k ] = M_dt_DSV * M_g * std::pow( M_Rick_x, 2 );
                }


                M_expo_r_y_vect.resize( S_y.size() );
                M_gamma_dt_DSV_y_.resize  ( S_y.size() );
                for ( UInt k = 0; k < S_y.size(); k++ )
                {
                    M_expo_r_y_vect[ k ] = M_expo_r1 * ( std::abs( S_y[ k ] ) >   .006 ) +
                                           M_expo_r2 * ( std::abs( S_y[ k ] ) <=  .006 );
                    
                    const auto M_Rick_y = ( M_fc0_greater_y[ k ] * std::pow( std::abs( S_y[ k ] ), .33 ) ) * ( std::abs( S_y[ k ] ) >  .006 ) +
                                          ( M_fc0_lower_y[ k ]   * std::pow( std::abs( S_y[ k ] ), .08 ) ) * ( std::abs( S_y[ k ] ) <= .006 );
                    M_gamma_dt_DSV_y_ [ k ] = M_dt_DSV * M_g * std::pow( M_Rick_y, 2 );
                }

            }
            else
            {
                std::cout << "No friction model inserted!!" << std::endl;
                exit( -1. );
            }


        }

    ~frictionClass() = default;  

    inline void f_x( const std::vector<Real>& H_interface,
                     const std::vector<Real>& u,
                     const std::vector<UInt>& idStaggeredInternalVectHorizontal,
                     const std::vector<UInt>& idStaggeredBoundaryVectWest,
                     const std::vector<UInt>& idStaggeredBoundaryVectEast )
    {
        
        switch ( M_frictionModel )
        {
            case 0:

                

                for ( const auto & Id : idStaggeredInternalVectHorizontal )
                {
                    alfa_x[ Id ] = 1;
                }
                

                
                for ( const auto & Id : idStaggeredBoundaryVectWest )
                {
                    alfa_x[ Id ] = 1;
                }
                
                
                for ( const auto & Id : idStaggeredBoundaryVectEast )
                {
                    alfa_x[ Id ] = 1;
                }
                


                break;

            case 1:
    
                
                for ( const auto & Id : idStaggeredInternalVectHorizontal )
                {
                    
                    Real alfa = 1.;
                    const Real den = std::pow( H_interface[ Id ], M_expo );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1. + M_gamma_dt_DSV_x * std::abs( u[ Id ] ) / den );
                    }
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( u[ Id ] ) << " " << std::pow( H_interface[ Id ], M_expo ) << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }

                    alfa_x[ Id ] = alfa;
                }
                

                
                for ( const auto & Id : idStaggeredBoundaryVectWest )
                {
                
                    Real alfa = 1.;
                    const Real den = std::pow( H_interface[ Id ], M_expo );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1. + M_gamma_dt_DSV_x * std::abs( u[ Id ] ) / den );
                    }
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( u[ Id ] ) << " " << std::pow( H_interface[ Id ], M_expo ) << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }
                    
                    alfa_x[ Id ] = alfa;
            
                }
                
                
                for ( const auto & Id : idStaggeredBoundaryVectEast )
                {

                    Real alfa = 1.;
                    const Real den = std::pow( H_interface[ Id ], M_expo );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1. + M_gamma_dt_DSV_x * std::abs( u[ Id ] ) / den );
                    }
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( u[ Id ] ) << " " << std::pow( H_interface[ Id ], M_expo ) << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }


                    alfa_x[ Id ] = alfa;

            
                }
                
                
                
                


                break;

            case 2:
                
                


                for ( const auto & Id : idStaggeredInternalVectHorizontal )
                {
                    
                    Real alfa = 1.;
                    
                    const Real den = std::pow( H_interface[ Id ], M_expo + M_expo_r_x_vect[ Id ] );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1. + M_gamma_dt_DSV_x_[ Id ] * std::pow( std::abs( u[ Id ] ), 1. - M_expo_r_x_vect[ Id ] ) / den );
                        alfa = std::min( alfa, 1. / ( 1. + M_gamma_dt_DSV_x * std::abs( u[ Id ] ) / std::pow( H_interface[ Id ], M_expo ) ) );
                    }
                    
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( u[ Id ] ) << " " << M_gamma_dt_DSV_x << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }

                    alfa_x[ Id ] = alfa;
                }
                

                
                for ( const auto & Id : idStaggeredBoundaryVectWest )
                {
                
                     
                    Real alfa = 1.;

                    const Real den = std::pow( H_interface[ Id ], M_expo + M_expo_r_x_vect[ Id ] );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1. + M_gamma_dt_DSV_x_[ Id ] * std::pow( std::abs( u[ Id ] ), 1. - M_expo_r_x_vect[ Id ] ) / den );
                        alfa = std::min( alfa, 1. / ( 1. + M_gamma_dt_DSV_x * std::abs( u[ Id ] ) / std::pow( H_interface[ Id ], M_expo ) ) );
                    }
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( u[ Id ] ) << " " << M_gamma_dt_DSV_x << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }
                    
                    alfa_x[ Id ] = alfa;
            
                }
                
                
                for ( const auto & Id : idStaggeredBoundaryVectEast )
                {
                    
                    Real alfa = 1.;
                    

                    const Real den = std::pow( H_interface[ Id ], M_expo + M_expo_r_x_vect[ Id ] );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1. + M_gamma_dt_DSV_x_[ Id ] * std::pow( std::abs( u[ Id ] ), 1. - M_expo_r_x_vect[ Id ] ) / den );
                        alfa = std::min( alfa, 1. / ( 1. + M_gamma_dt_DSV_x * std::abs( u[ Id ] ) / std::pow( H_interface[ Id ], M_expo ) ) );
                    }
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( u[ Id ] ) << " " << M_gamma_dt_DSV_x << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }
                    
                   

                    alfa_x[ Id ] = alfa;
            
                }
                
                
                
                



                break;

            default:
                std::cout << "Error in frictionClass" << std::endl;
                exit( -1. );
        }

    }


    inline void f_y( const std::vector<Real>& H_interface,
                     const std::vector<Real>& v,
                     const std::vector<UInt>& idStaggeredInternalVectVertical,
                     const std::vector<UInt>& idStaggeredBoundaryVectNorth,
                     const std::vector<UInt>& idStaggeredBoundaryVectSouth )
    {

        switch ( M_frictionModel )
        {
            case 0:


                for ( const auto & Id : idStaggeredInternalVectVertical )
                {
                    alfa_y[ Id ] = 1.;
                }

                
                for ( const auto & Id : idStaggeredBoundaryVectNorth )
                {
                    alfa_y[ Id ] = 1.;
                }
            

                for ( const auto & Id : idStaggeredBoundaryVectSouth )
                {
                    alfa_y[ Id ] = 1.;
                }
                


                break;

            case 1: 

                for ( const auto & Id : idStaggeredInternalVectVertical )
                {
                    Real alfa = 1.;
                    const Real den = std::pow( H_interface[ Id ], M_expo );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1 + M_gamma_dt_DSV_y * std::abs( v[ Id ] ) / den );
                    }
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( v[ Id ] ) << " " << std::pow( H_interface[ Id ], M_expo ) << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }


                    alfa_y[ Id ] = alfa;
                }

                
                for ( const auto & Id : idStaggeredBoundaryVectNorth )
                {
                    Real alfa = 1.;
                    const Real den = std::pow( H_interface[ Id ], M_expo );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1 + M_gamma_dt_DSV_y * std::abs( v[ Id ] ) / den );
                    }
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( v[ Id ] ) << " " << std::pow( H_interface[ Id ], M_expo ) << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }


                    alfa_y[ Id ] = alfa;
                }
            

                for ( const auto & Id : idStaggeredBoundaryVectSouth )
                {
                    Real alfa = 1.;
                    const Real den = std::pow( H_interface[ Id ], M_expo );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1 + M_gamma_dt_DSV_y * std::abs( v[ Id ] ) / den );
                    }
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( v[ Id ] ) << " " << std::pow( H_interface[ Id ], M_expo ) << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }


                    alfa_y[ Id ] = alfa;
                }
                
                
                


                break;

            case 2:

                for ( const auto & Id : idStaggeredInternalVectVertical )
                {
                    
                    Real alfa = 1.;
                    
                    const Real den = std::pow( H_interface[ Id ], M_expo + M_expo_r_y_vect[ Id ] );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1 + M_gamma_dt_DSV_y_[ Id ] * std::pow( std::abs( v[ Id ] ), 1. - M_expo_r_y_vect[ Id ] ) / den );
                        alfa = std::min( alfa, 1. / ( 1. + M_gamma_dt_DSV_y * std::abs( v[ Id ] ) / std::pow( H_interface[ Id ], M_expo ) ) );
                    }
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( v[ Id ] ) << " " << M_gamma_dt_DSV_y << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }
                    
                    
                    
                    alfa_y[ Id ] = alfa;
                }

                
                for ( const auto & Id : idStaggeredBoundaryVectNorth )
                {
                    
                    Real alfa = 1.;
                    
                    const Real den = std::pow( H_interface[ Id ], M_expo + M_expo_r_y_vect[ Id ] );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1 + M_gamma_dt_DSV_y_[ Id ] * std::pow( std::abs( v[ Id ] ), 1. - M_expo_r_y_vect[ Id ] ) / den );
                        alfa = std::min( alfa, 1. / ( 1. + M_gamma_dt_DSV_y * std::abs( v[ Id ] ) / std::pow( H_interface[ Id ], M_expo ) ) );
                    }
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( v[ Id ] ) << " " << M_gamma_dt_DSV_y << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }
                    
                    alfa_y[ Id ] = alfa;
                }
            
                
                for ( const auto & Id : idStaggeredBoundaryVectSouth )
                {
                    
                    Real alfa = 1.;
                    
                    const Real den = std::pow( H_interface[ Id ], M_expo + M_expo_r_y_vect[ Id ] );
                    if ( den > M_H_min )
                    {
                        alfa = 1. / ( 1 + M_gamma_dt_DSV_y_[ Id ] * std::pow( std::abs( v[ Id ] ), 1. - M_expo_r_y_vect[ Id ] ) / den );
                        alfa = std::min( alfa, 1. / ( 1. + M_gamma_dt_DSV_y * std::abs( v[ Id ] ) / std::pow( H_interface[ Id ], M_expo ) ) );
                    }
                    
                    if ( std::isnan( alfa ) )
                    {
                        std::cout << std::abs( v[ Id ] ) << " " << M_gamma_dt_DSV_y << std::endl;
                        std::cout << "nan in frictionClass" << std::endl;
                        exit( -1. );
                    }
                    
                    alfa_y[ Id ] = alfa;
                }
            
                

                   

                break;

            default:

                std::cout << "Error in frictionClass" << std::endl;
                exit( -1. );

        }

    }

    



    std::vector<Real> alfa_x,
                      alfa_y;

private:

          UInt        M_frictionModel;
    
    const Real        M_n_manning,
                      M_dt_DSV;
    
    Real              M_gamma_dt_DSV_x,
                      M_gamma_dt_DSV_y;
    
    Real M_H_min;

    const UInt N_rows,
               N_cols;

    const Real M_expo    = 4./3.;
    const Real M_expo_r1 = .11 * 2;
    const Real M_expo_r2 = .03 * 2;
    const Real M_g       = 9.81;
    const Real M_toll    = 1.e-4;

          std::vector<Real> M_expo_r_x_vect,
                            M_expo_r_y_vect,
                            M_gamma_dt_DSV_x_,
                            M_gamma_dt_DSV_y_;
    
};







class upwind
{

public:

    upwind( const std::vector<Real>& u, const std::vector<Real>& v, const Real& N_rows, const Real& N_cols ):
    N_cols( N_cols ),
    N_rows( N_rows )
    {
        horizontal.resize( u.size() );
        vertical.  resize( v.size() );
    }

    ~upwind() = default;    

    inline void computeHorizontal( const Eigen::VectorXd&   H,
                                   const std::vector<Real>& u,
                                   const std::vector<UInt>& idStaggeredInternalVectHorizontal,
                                   const std::vector<UInt>& idStaggeredBoundaryVectWest,
                                   const std::vector<UInt>& idStaggeredBoundaryVectEast )
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
                                        

            horizontal[ Id ] = ( H_left + H_right ) * .5 + signum( u[ Id + 1 ] ) * ( H_left - H_right ) * .5; // messo il segno della velocità interna di fianco
    
        }
        
        
        for ( const auto & Id : idStaggeredBoundaryVectEast )
        {

            const UInt i = Id / ( N_cols + 1 );
          


            const Real & H_left  = H( Id - i - 1 ),
                         H_right = 0;
                                        

            horizontal[ Id ] = ( H_left + H_right ) * .5 + signum( u[ Id - 1 ] ) * ( H_left - H_right ) * .5; // messo il segno della velocità interna di fianco

    
        }
        

 

    }


    inline void computeVertical( const Eigen::VectorXd&   H,
                                 const std::vector<Real>& v,
                                 const std::vector<UInt>& idStaggeredInternalVectVertical,
                                 const std::vector<UInt>& idStaggeredBoundaryVectNorth,
                                 const std::vector<UInt>& idStaggeredBoundaryVectSouth )
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


    std::vector<Real> horizontal,
                      vertical;

private:

    UInt N_cols,
         N_rows;

};



inline void bilinearInterpolation( const std::vector<Real>& u,
                                   const std::vector<Real>& v,
                                         std::vector<Real>& u_star,
                                         std::vector<Real>& v_star,
                                   const UInt&              nrows,
                                   const UInt&              ncols,
                                   const Real&              dt,
                                   const Real&              pixel_size
                                 )
{


    // +-----------------------------------------------+
    // |              Horizontal Velocity              |
    // +-----------------------------------------------+

    
    for ( UInt i = 0; i < nrows; i++ )
    {
        
        for ( UInt j = 1; j < ncols; j++ )
        {
                
            const auto Id    = j + i * ( ncols + 1 ), // u

                       // ID della velocità
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

                       // ID della velocità
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



inline void bilinearInterpolation( const std::vector<Real>& u,
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
                                   const std::vector<UInt>& idStaggeredBoundaryVectSouth
                                 )
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




inline void computeAdjacencies( const std::vector<Real>& basin_mask_Vec,
                               
                                      std::vector<UInt>& idStaggeredBoundaryVectSouth, 
                                      std::vector<UInt>& idStaggeredBoundaryVectNorth,
                                      std::vector<UInt>& idStaggeredBoundaryVectWest,
                                      std::vector<UInt>& idStaggeredBoundaryVectEast,
                               
                                      std::vector<UInt>& idStaggeredInternalVectHorizontal,
                                      std::vector<UInt>& idStaggeredInternalVectVertical,
                               
                                      std::vector<UInt>& idBasinVect,
                                      std::vector<UInt>& idBasinVectReIndex,
                                
                                const UInt&              N_rows,
                                const UInt&              N_cols )
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






inline void buildMatrix( const std::vector<Real>& H_int_x,
                         const std::vector<Real>& H_int_y,
                         const std::vector<Real>& orography,
                         const std::vector<Real>& u_star,
                         const std::vector<Real>& v_star,
                         const Eigen::VectorXd&   H,
                         const UInt&              N_cols,
                         const UInt&              N_rows,
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
                         const std::vector<UInt>& idBasinVectReIndex,
                         const bool&              isNonReflectingBC,
                         
                               std::vector<Eigen::Triplet<Real> >& coefficients,
                               Eigen::VectorXd&                    rhs
                   
                )


{

    coefficients.reserve( idBasinVect.size() +
                      4 * idStaggeredInternalVectHorizontal.size() +
                      4 * idStaggeredInternalVectVertical.size() );

    
    

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
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDleftReIndex, - c3 * coeff_m ) );
            
            
            coefficients.push_back( Eigen::Triplet<Real>( IDleftReIndex,  IDleftReIndex,   c3 * coeff_m ) );
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDrightReIndex,  c3 * coeff_m ) );
            

         
            rhs( IDleftReIndex )  += - c1 * ( + coeff_m  * u_star[ Id ] ) - ( c3 * coeff_m * orography[ IDleft ]
                                                                            - c3 * coeff_m * orography[ IDright ] );
            
            rhs( IDrightReIndex ) += - c1 * ( - coeff_m  * u_star[ Id ] ) - ( c3 * coeff_m * orography[ IDright ]
                                                                            - c3 * coeff_m * orography[ IDleft ] );
             
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
                
            rhs( IDrightReIndex ) += isNonReflectingBC *
                                   ( - c1 * ( - coeff_m  * u_star[ Id ] ) - ( c3 * coeff_m * orography[ IDrightright ]
                                                                            - c3 * coeff_m * orography[ IDright ] ) );
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
        
            rhs( IDleftReIndex ) += isNonReflectingBC *
                                  ( - c1 * ( + coeff_m  * u_star[ Id ] ) - ( c3 * coeff_m * orography[ IDleftleft ]
                                                                           - c3 * coeff_m * orography[ IDleft ] ) );
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
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDleftReIndex, - c3 * coeff_m ) );
            
            
            coefficients.push_back( Eigen::Triplet<Real>( IDleftReIndex,  IDleftReIndex,   c3 * coeff_m ) );
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDrightReIndex,  c3 * coeff_m ) );
            
            
            rhs( IDleftReIndex )  += - c1 * ( + coeff_m  * v_star[ Id ] ) - ( c3 * coeff_m * orography[ IDleft ]
                                                                            - c3 * coeff_m * orography[ IDright ] );
            
            rhs( IDrightReIndex ) += - c1 * ( - coeff_m  * v_star[ Id ] ) - ( c3 * coeff_m * orography[ IDright ]
                                                                            - c3 * coeff_m * orography[ IDleft ] );
             
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

            rhs( IDrightReIndex ) += isNonReflectingBC *
                                   ( - c1 * ( - coeff_m  * v_star[ Id ] ) - ( c3 * coeff_m * orography[ IDrightright ]
                                                                            - c3 * coeff_m * orography[ IDright ] ) );
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

            rhs( IDleftReIndex )  += isNonReflectingBC *
                                   ( - c1 * ( + coeff_m  * v_star[ Id ] ) - ( c3 * coeff_m * orography[ IDleftleft ]
                                                                            - c3 * coeff_m * orography[ IDleft ] ) );
        }
        

    
    
    }
    
    

}





inline void buildMatrix( const std::vector<Real>& H_int_x,
                         const std::vector<Real>& H_int_y,
                         const std::vector<Real>& orography,
                         const std::vector<Real>& u_star,
                         const std::vector<Real>& v_star,
                         const Eigen::VectorXd&   H,
                         const UInt&              N_cols,
                         const UInt&              N_rows,
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
                         const std::vector<UInt>& idBasinVectReIndex,
                         const bool&              isNonReflectingBC,
                         const Real&              alp,  // 0 < alp <= 1
                         
                               std::vector<Eigen::Triplet<Real> >& coefficients,
                               Eigen::VectorXd&                    rhs
                   
                )


{

    coefficients.reserve( idBasinVect.size() +
                      4 * idStaggeredInternalVectHorizontal.size() +
                      4 * idStaggeredInternalVectVertical.size() );

    
    

    for ( const auto & Id : idBasinVect )
    {
        const auto IDreIndex = idBasinVectReIndex[ Id ];
        coefficients.push_back( Eigen::Triplet<Real>( IDreIndex,  IDreIndex,  1. ) );
        rhs( IDreIndex ) = ( H( Id ) + precipitation[ Id ] * dt_DSV ) * alp;
//        if (precipitation[ Id ] != 0) std::cout << precipitation[ Id ] << std::endl;
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
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDleftReIndex, - c3 * coeff_m ) );
            
            
            coefficients.push_back( Eigen::Triplet<Real>( IDleftReIndex,  IDleftReIndex,   c3 * coeff_m ) );
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDrightReIndex,  c3 * coeff_m ) );
            

         
            rhs( IDleftReIndex )  += - c1 * ( + coeff_m  * u_star[ Id ] ) - ( orography[ IDleft ] - orography[ IDright ] ) * c3 * coeff_m * alp;
            
            rhs( IDrightReIndex ) += - c1 * ( - coeff_m  * u_star[ Id ] ) - ( orography[ IDright ] - orography[ IDleft ] ) * c3 * coeff_m * alp;
             
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
                
            rhs( IDrightReIndex ) += isNonReflectingBC *
                                   ( - c1 * ( - coeff_m  * u_star[ Id ] ) - ( orography[ IDrightright ] - orography[ IDright ] ) * c3 * coeff_m * alp );
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
        
            rhs( IDleftReIndex ) += isNonReflectingBC *
                                  ( - c1 * ( + coeff_m  * u_star[ Id ] ) - ( orography[ IDleftleft ] - orography[ IDleft ] ) * c3 * coeff_m * alp );
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
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDleftReIndex, - c3 * coeff_m ) );
            
            
            coefficients.push_back( Eigen::Triplet<Real>( IDleftReIndex,  IDleftReIndex,   c3 * coeff_m ) );
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDrightReIndex,  c3 * coeff_m ) );
            
            
            rhs( IDleftReIndex )  += - c1 * ( + coeff_m  * v_star[ Id ] ) - ( orography[ IDleft ] - orography[ IDright ] ) * c3 * coeff_m * alp;
            
            rhs( IDrightReIndex ) += - c1 * ( - coeff_m  * v_star[ Id ] ) - ( orography[ IDright ] - orography[ IDleft ] ) * c3 * coeff_m * alp;
             
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

            rhs( IDrightReIndex ) += isNonReflectingBC *
                                   ( - c1 * ( - coeff_m  * v_star[ Id ] ) - ( orography[ IDrightright ] - orography[ IDright ] ) * c3 * coeff_m * alp );
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

            rhs( IDleftReIndex )  += isNonReflectingBC *
                                   ( - c1 * ( + coeff_m  * v_star[ Id ] ) - ( orography[ IDleftleft ] - orography[ IDleft ] ) * c3 * coeff_m * alp );
        }
        

    
    
    }
    
    

}





/*
inline void buildMatrix( const std::vector<Real>& H_int_x,
                         const std::vector<Real>& H_int_y,
                         const std::vector<Real>& orography,
                         const std::vector<Real>& u_star,
                         const std::vector<Real>& v_star,
                         const Eigen::VectorXd&   eta,
                         const UInt&              N_cols,
                         const UInt&              N_rows,
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
                         const std::vector<UInt>& idBasinVectReIndex,
                         const bool&              isNonReflectingBC,
                         
                               std::vector<Eigen::Triplet<Real> >& coefficients,
                               Eigen::VectorXd&                    rhs
                   
                )


{

    coefficients.reserve( idBasinVect.size() +
                      4 * idStaggeredInternalVectHorizontal.size() +
                      4 * idStaggeredInternalVectVertical.size() );

    
    

    for ( const auto & Id : idBasinVect )
    {
        const auto IDreIndex = idBasinVectReIndex[ Id ];
        coefficients.push_back( Eigen::Triplet<Real>( IDreIndex,  IDreIndex,  1. ) );
        rhs( IDreIndex ) = eta( Id ) + precipitation[ Id ] * dt_DSV;
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
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDleftReIndex, - c3 * coeff_m ) );
            
            
            coefficients.push_back( Eigen::Triplet<Real>( IDleftReIndex,  IDleftReIndex,   c3 * coeff_m ) );
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDrightReIndex,  c3 * coeff_m ) );
            

         
            rhs( IDleftReIndex )  += - c1 * ( + coeff_m  * u_star[ Id ] );
            
            rhs( IDrightReIndex ) += - c1 * ( - coeff_m  * u_star[ Id ] );
             
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
                
            rhs( IDrightReIndex ) += isNonReflectingBC *
                                   ( - c1 * ( - coeff_m  * u_star[ Id ] ) - ( c3 * coeff_m * orography[ IDrightright ]
                                                                            - c3 * coeff_m * orography[ IDright ] ) );
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
        
            rhs( IDleftReIndex ) += isNonReflectingBC *
                                  ( - c1 * ( + coeff_m  * u_star[ Id ] ) - ( c3 * coeff_m * orography[ IDleftleft ]
                                                                           - c3 * coeff_m * orography[ IDleft ] ) );
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
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDleftReIndex, - c3 * coeff_m ) );
            
            
            coefficients.push_back( Eigen::Triplet<Real>( IDleftReIndex,  IDleftReIndex,   c3 * coeff_m ) );
            coefficients.push_back( Eigen::Triplet<Real>( IDrightReIndex, IDrightReIndex,  c3 * coeff_m ) );
            
            
            rhs( IDleftReIndex )  += - c1 * ( + coeff_m  * v_star[ Id ] );
            
            rhs( IDrightReIndex ) += - c1 * ( - coeff_m  * v_star[ Id ] );
             
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

            rhs( IDrightReIndex ) += isNonReflectingBC *
                                   ( - c1 * ( - coeff_m  * v_star[ Id ] ) - ( c3 * coeff_m * orography[ IDrightright ]
                                                                            - c3 * coeff_m * orography[ IDright ] ) );
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

            rhs( IDleftReIndex )  += isNonReflectingBC *
                                   ( - c1 * ( + coeff_m  * v_star[ Id ] ) - ( c3 * coeff_m * orography[ IDleftleft ]
                                                                            - c3 * coeff_m * orography[ IDleft ] ) );
        }
        

    
    
    }
    
    

}
*/









inline void updateVel(       std::vector<Real>& u,
                             std::vector<Real>& v,
                       const std::vector<Real>& u_star,
                       const std::vector<Real>& v_star,
                       const std::vector<Real>& H_int_x,
                       const std::vector<Real>& H_int_y,
                       const std::vector<Real>& alfa_x,
                       const std::vector<Real>& alfa_y,
                       const Real&              N_rows,
                       const Real&              N_cols,
                       const Real&              c2,
                       const Real&              H_min,
                       const Eigen::VectorXd&   eta,
                       const std::vector<Real>& orography,
                       const std::vector<UInt>& idStaggeredInternalVectHorizontal,
                       const std::vector<UInt>& idStaggeredInternalVectVertical,
                       const std::vector<UInt>& idStaggeredBoundaryVectWest,
                       const std::vector<UInt>& idStaggeredBoundaryVectEast,
                       const std::vector<UInt>& idStaggeredBoundaryVectNorth,
                       const std::vector<UInt>& idStaggeredBoundaryVectSouth,
                       const bool&              isNonReflectingBC
                     )
{

    
    // +-----------------------------------------------+
    // |              Update Vertical Velocity         |
    // +-----------------------------------------------+

    for ( const auto & Id : idStaggeredInternalVectVertical )
    {
        
        const UInt IDsouth = Id,
                   IDnorth = Id - N_cols;



        const auto & H_interface = H_int_y[ Id ];
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
    
        const auto & H_interface = H_int_y[ Id ];
        if ( H_interface > H_min )
        {
//            std::cout << H_interface << std::endl;
            v[ Id ] = isNonReflectingBC * alfa_y[ Id ] * ( v_star[ Id ] - c2 * ( orography[ Id + N_cols ] - orography[ Id ] ) );
            
        }
        else
        {
            v[ Id ] = 0.;
        }
        

           
    
    }


    // last row
    for ( const auto & Id : idStaggeredBoundaryVectSouth )
    {
        
        const auto & H_interface = H_int_y[ Id ];
        if ( H_interface > H_min )
        {
            v[ Id ] = isNonReflectingBC * alfa_y[ Id ] * ( v_star[ Id ] - c2 * ( orography[ Id - N_cols ] - orography[ Id - 2*N_cols ] ) );
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
        

        const auto & H_interface = H_int_x[ Id ];
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
        
        const UInt i = Id / ( N_cols + 1 );

        const auto & H_interface = H_int_x[ Id ];
        if ( H_interface > H_min )
        {
            u[ Id ] = isNonReflectingBC * alfa_x[ Id ] * ( u_star[ Id ] - c2 * ( orography[ Id - i + 1 ] - orography[ Id - i ] ) );
        }
        else
        {
            u[ Id ] = 0.;
        }



    }
    
    
    for ( const auto & Id : idStaggeredBoundaryVectEast )
    {

        const UInt i = Id / ( N_cols + 1 );

        
        const auto & H_interface = H_int_x[ Id ];
        if ( H_interface > H_min )
        {
            u[ Id ] = isNonReflectingBC * alfa_x[ Id ] * ( u_star[ Id ] - c2 * ( orography[ Id - i - 1 ] - orography[ Id - i - 2 ] ) );
        }
        else
        {
            u[ Id ] = 0.;
        }

    }
    
    
    
}









inline Real maxCourant( const std::vector<Real>& u, 
                        const std::vector<Real>& v,
                        const Real&              c1 )
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


inline Real maxCourant( const Eigen::VectorXd& H,
                        const Real&            gravity,
                        const Real&            c1 )
{

    const Real Courant_cel = std::sqrt( H.maxCoeff() * gravity );
    return( Courant_cel * c1 );
}




inline Real compute_dt_sediment( const Real&              alpha,
                                 const Real&              beta,
                                 const Real&              S_x,
                                 const Real&              S_y,
                                 const std::vector<Real>& u, 
                                 const std::vector<Real>& v, 
                                 const Real&              pixel_size,
                                 const Real&              dt_DSV,
                                       UInt&              numberOfSteps )
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

    //std::cout << dt_sed << " " << std::floor( dt_DSV / dt_sed ) << " " << dt_DSV / dt_sed << " " << dt_DSV / std::floor( dt_DSV / dt_sed ) << std::endl;
    
    dt_sed = std::min( dt_DSV / std::floor( dt_DSV / dt_sed ), dt_DSV );

    numberOfSteps = std::floor( dt_DSV / dt_sed );

    return( dt_sed );

    
}


inline void saveVector( const Eigen::VectorXd& b,
                        const std::string& Name )
{
    std::ofstream ff( Name );
    
    for ( UInt k = 0; k < b.size(); k++ )
    {
        ff << b( k ) << " ";
        ff << std::endl;
    }
        
    ff.close();
    
}


inline void saveMatrix( const SpMat& A,
                        const std::string& Name )
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


inline void saveSolution( const std::string& preName,
                          const std::string& flag,
                          const UInt& N_rows,
                          const UInt& N_cols,
                          const Real& xllcorner,
                          const Real& yllcorner,
                          const Real& cellsize,
                          const Real& NODATA_value,
                          const Eigen::VectorXd& H ) // it is H or orography
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



inline void saveSolution( const std::string& preName,
                          const std::string& flag,
                          const UInt& N_rows,
                          const UInt& N_cols,
                          const Real& xllcorner,
                          const Real& yllcorner,
                          const Real& cellsize,
                          const Real& NODATA_value,
                          const std::vector<Real>& H ) // it is H or orography
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


inline void saveSolution( const std::string& preName,
                          const std::string& flag,
                          const UInt& N_rows,
                          const UInt& N_cols,
                          const Real& xllcorner,
                          const Real& yllcorner,
                          const Real& cellsize,
                          const Real& NODATA_value,
                          const std::vector<Int>& H ) // it is H or orography
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





inline void saveSolution( const std::string& preName, 
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
                          const Eigen::VectorXd& H ) // it is H or orography
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



inline void saveSolution( const std::string& preName,
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
                          const std::vector<Real>& H ) // it is H or orography
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




inline void saveTemporalSequence( const Vector2D&    X_gauges,
                                  const Real&        dt,
                                  const UInt&        n,
                                  const std::string& preName,
                                  const Real&        H ) 
{
    std::ofstream ff( preName + ".txt", std::ofstream::out | std::ofstream::app );
    if ( n == 1 )
    {
        ff << X_gauges( 0 ) << " " << X_gauges( 1 ) << std::endl;
        ff << dt << std::endl;
    }
    ff << H << std::endl;
    ff.close();
}



inline void saveTemporalSequence( const Real&        dt,
                                  const UInt&        n,
                                  const std::string& preName,
                                  const Real&        H )
{
    std::ofstream ff( preName + ".txt", std::ofstream::out | std::ofstream::app );
    if ( n == 1 )
    {
        ff << dt << std::endl;
    }
    ff << H << std::endl;        
    ff.close();
}



// For gravitational layer
inline void computeResiduals( const std::vector<Real>& n_x,
                              const std::vector<Real>& n_y,
                              const UInt&              N_cols,
                              const UInt&              N_rows,
                              const std::vector<Real>& h,
                              const std::vector<Real>& coeff, // conducibilità idraulica
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
                                    std::vector<Real>& Res_y )
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
inline void computeResidualsTruncated( const std::vector<Real>&                 u,
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
    
    /*
    for ( const auto & Id : idStaggeredBoundaryVectEast )
    {
        const Real coeff_right = c1 * alpha * std::pow( std::abs( S_x[ Id ] ), beta ) * u[ Id ]
                                              * ( .5 - .5 * signum( u[ Id ] ) ),
        
                   coeff_left = c1 * alpha * std::pow( std::abs( S_x[ Id ] ), beta ) * u[ Id ]
                                              * ( .5 + .5 * signum( u[ Id ] ) );
                                
        
        Gamma_x[ Id ] = std::array<Real,2>{{ coeff_right, coeff_left }};
    }*/
    
    
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
    
    
    /*
    for ( UInt Id = 0; Id < N_rows * ( N_cols + 1 ); Id++ )
    {
        const Real coeff_right = c1 * alpha * std::pow( S_x[ Id ], beta ) * u[ Id ]
                                              * ( .5 - .5 * signum( u[ Id ] ) ),
        
                   coeff_left = c1 * alpha * std::pow( S_x[ Id ], beta ) * u[ Id ]
                                              * ( .5 + .5 * signum( u[ Id ] ) );
                                
        
        Gamma_x[ Id ] = std::array<Real,2>{{ coeff_right, coeff_left }};
    }
    

    
    for ( UInt Id = 0; Id < N_cols * ( N_rows + 1 ); Id++ )
    {
        const Real coeff_right = c1 * alpha * std::pow( S_y[ Id ], beta ) * v[ Id ]
                                              * ( .5 - .5 * signum( v[ Id ] ) ),
        
                   coeff_left = c1 * alpha * std::pow( S_y[ Id ], beta ) * v[ Id ]
                                              * ( .5 + .5 * signum( v[ Id ] ) );
                                
        
        Gamma_y[ Id ] = std::array<Real,2>{{ coeff_right, coeff_left }};
    }*/

    


}



inline std::vector<Real> compute_d_perc( const std::vector<Real>& clay, const std::vector<Real>& sand, const Real& perc )
{
    // linear interpolation in log10 x-scale
    std::vector<Real> d_perc( clay.size() );
    
    for ( UInt i = 0; i < clay.size(); i++ )
    {
        auto & d_perc_cell = d_perc[ i ];
        
        const auto & sand_cell = sand[ i ];
        const auto & clay_cell = clay[ i ];
        
        const auto Y_0 = ( 1 - sand_cell ) * 100;
        const auto Y_1 = clay_cell * 100;
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











