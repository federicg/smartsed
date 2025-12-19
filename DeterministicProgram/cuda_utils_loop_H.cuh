#pragma once

//! set CUDA headers
#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <cuda_runtime_api.h>

//! Thrust
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>


#define CHECK_CUDA(func)                                                       \
{                                                                              \
    cudaError_t status = (func);                                               \
    if (status != cudaSuccess) {                                               \
        printf("CUDA API failed at line %d with error: %s (%d)\n",             \
               __LINE__, cudaGetErrorString(status), status);                  \
        return EXIT_FAILURE;                                                   \
    }                                                                          \
}

#define CHECK_CUSPARSE(func)                                                   \
{                                                                              \
    cusparseStatus_t status = (func);                                          \
    if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
        printf("cuSPARSE API failed at line %d with error: %s (%d)\n",         \
               __LINE__, cusparseGetErrorString(status), status);              \
        return EXIT_FAILURE;                                                   \
    }                                                                          \
}

#define CHECK_CUBLAS(func)                                                     \
{                                                                              \
    cublasStatus_t status = (func);                                            \
    if (status != CUBLAS_STATUS_SUCCESS) {                                     \
        printf("CUBLAS API failed at line %d with error: %d\n",                \
               __LINE__, status);                                              \
        return EXIT_FAILURE;                                                   \
    }                                                                          \
}

#if defined(NDEBUG)
#   define PRINT_INFO(var)
#else
#   define PRINT_INFO(var) printf("  " #var ": %f\n", var);
#endif

//==============================================================================
#ifdef __CUDACC__
template <typename Kernel, typename... Args>
inline void launch_kernel(Kernel kernel, int n, Args&&... args) {
    constexpr int blockSize = 256;
    int numBlocks = (n + blockSize - 1) / blockSize;
    kernel<<<numBlocks, blockSize>>>(std::forward<Args>(args)..., n);
}
#endif
//==============================================================================
// Horizontal upwind
__global__ void computeHorizontalInternalKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* u,
    double* horizontal,
    unsigned int N_cols,
    unsigned int n);

__global__ void computeHorizontalWestKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* u,
    double* horizontal,
    unsigned int N_cols,
    unsigned int n);

__global__ void computeHorizontalEastKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* u,
    double* horizontal,
    unsigned int N_cols,
    unsigned int n);

void compute_horizontal_interface_wrapper(const thrust::device_vector<unsigned int>& idStaggeredInternalVectHorizontal, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectWest, const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectEast,
		const thrust::device_vector<double>& H, const thrust::device_vector<double>& u, thrust::device_vector<double>& horizontal, const unsigned int N_cols);

// Vertical upwind
__global__ void computeVerticalInternalKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* v,
    double* vertical,
    unsigned int N_cols,
    unsigned int n);

__global__ void computeVerticalNorthKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* v,
    double* vertical,
    unsigned int N_cols,
    unsigned int n);

__global__ void computeVerticalSouthKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* v,
    double* vertical,
    unsigned int N_cols,
    unsigned int n);

void compute_vertical_interface_wrapper(const thrust::device_vector<unsigned int>& idStaggeredInternalVectVertical, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectNorth, const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectSouth,
		const thrust::device_vector<double>& H, const thrust::device_vector<double>& v, thrust::device_vector<double>& vertical, const unsigned int N_cols);

//==============================================================================
// Horizontal friction
__global__ void computeHorizontalInternalKernel_friction(
    const unsigned int* ids,
    const double* H_interface_horizontal,
    const double* u,
    const double* M_expo_r_x_vect,
    double* alfa_x,
    const double* M_gamma_dt_DSV_x_,
    double M_dt_DSV, double M_coeff, double M_H_min,
    double M_expo, unsigned int M_frictionModel,
    unsigned int n);

__global__ void computeHorizontalWestKernel_friction(
    const unsigned int* ids,
    const double* H_interface_horizontal,
    const double* u,
    const double* M_expo_r_x_vect,
    double* alfa_x,
    const double* M_gamma_dt_DSV_x_,
    double M_dt_DSV, double M_coeff, double M_H_min,
    double M_expo, unsigned int M_frictionModel,
    unsigned int n);

__global__ void computeHorizontalEastKernel_friction(
    const unsigned int* ids,
    const double* H,
    const double* u,
    double* horizontal,
    unsigned int N_cols,
    unsigned int n);

void compute_horizontal_friction_wrapper(const thrust::device_vector<unsigned int>& idStaggeredInternalVectHorizontal, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectWest, const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectEast,
		const thrust::device_vector<double>& H_interface_horizontal, const thrust::device_vector<double>& u, const thrust::device_vector<double>& M_expo_r_x_vect, 
		      thrust::device_vector<double>& alfa_x, const thrust::device_vector<double>& M_gamma_dt_DSV_x_, double M_dt_DSV, double M_coeff, double M_H_min, 
                double M_expo, unsigned int M_frictionModel); 

// Vertical friction
__global__ void computeVerticalInternalKernel_friction(
    const unsigned int* ids,
    const double* H_interface_vertical,
    const double* v,
    const double* M_expo_r_y_vect,
    double* alfa_y,
    const double* M_gamma_dt_DSV_y_,
    double M_dt_DSV, double M_coeff, double M_H_min,
    double M_expo, unsigned int M_frictionModel,
    unsigned int n);

__global__ void computeVerticalNorthKernel_friction(
    const unsigned int* ids,
    const double* H_interface_vertical,
    const double* v,
    const double* M_expo_r_y_vect,
    double* alfa_y,
    const double* M_gamma_dt_DSV_y_,
    double M_dt_DSV, double M_coeff, double M_H_min,
    double M_expo, unsigned int M_frictionModel,
    unsigned int n);

__global__ void computeVerticalSouthKernel_friction(
    const unsigned int* ids,
    const double* H_interface_vertical,
    const double* v,
    const double* M_expo_r_y_vect,
    double* alfa_y,
    const double* M_gamma_dt_DSV_y_,
    double M_dt_DSV, double M_coeff, double M_H_min,
    double M_expo, unsigned int M_frictionModel,
    unsigned int n);

void compute_vertical_friction_wrapper(const thrust::device_vector<unsigned int>& idStaggeredInternalVectVertical, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectNorth, const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectSouth,
		const thrust::device_vector<double>& H_interface_vertical, const thrust::device_vector<double>& v, const thrust::device_vector<double>& M_expo_r_y_vect, 
		      thrust::device_vector<double>& alfa_y, const thrust::device_vector<double>& M_gamma_dt_DSV_y_, double M_dt_DSV, double M_coeff, double M_H_min, 
                double M_expo, unsigned int M_frictionModel); 

//==============================================================================

typedef struct VecStruct {
    cusparseDnVecDescr_t vec;
    double*              ptr;
} Vec;

//==============================================================================

/// A 5-point Laplacian on a g x g grid with Dirichlet boundary conditions.
/// This code allocates. The caller must free.
void make_laplace_matrix(int * n_out,
                         int **row_offsets_out, 
                         int **columns_out, 
                         double **values_out);

//==============================================================================

int gpu_CG(cublasHandle_t       cublasHandle,
           cusparseHandle_t     cusparseHandle,
           int                  m,
           cusparseSpMatDescr_t matA,
           cusparseSpMatDescr_t matL,
           Vec                  d_B,
           Vec                  d_X,
           Vec                  d_R,
           Vec                  d_R_aux,
           Vec                  d_P,
           Vec                  d_T,
           Vec                  d_tmp,
           void*                d_bufferMV,
           int                  maxIterations,
           double               tolerance);

//==============================================================================
//==============================================================================





