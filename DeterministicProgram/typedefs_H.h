#ifndef TYPDEFS_H_H
#  define TYPDEFS_H_H

  //! std library
#include <cstdint>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>


//! Eigen library    
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

//! IML++ CG template
#include "cg.hpp"

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
#  define M_PI 3.14159265358979323846
#endif

using functionType = std::function <void (std::vector<Real>&, 
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
                                          const std::vector<UInt>&)>;

#endif
