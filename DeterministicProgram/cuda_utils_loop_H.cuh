#pragma once

#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <cublas_v2.h>

#include "utils_H.h"

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

static constexpr int blockSize = 256;

#define LAUNCH_KERNEL(kernel, n, ...)                    \
  do {                                                   \
    int numBlocks = ((n) + blockSize - 1) / blockSize;  \
    kernel<<<numBlocks, blockSize>>>(__VA_ARGS__, n);   \
  } while (0)

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

//==============================================================================
// Horizontal friction
__global__ void computeHorizontalInternalKernel_friction(
    const unsigned int* ids,
    const double* H,
    const double* u,
    double* horizontal,
    unsigned int N_cols,
    unsigned int n);

__global__ void computeHorizontalWestKernel_friction(
    const unsigned int* ids,
    const double* H,
    const double* u,
    double* horizontal,
    unsigned int N_cols,
    unsigned int n);

__global__ void computeHorizontalEastKernel_friction(
    const unsigned int* ids,
    const double* H,
    const double* u,
    double* horizontal,
    unsigned int N_cols,
    unsigned int n);

// Vertical friciton
__global__ void computeVerticalInternalKernel_friction(
    const unsigned int* ids,
    const double* H,
    const double* v,
    double* vertical,
    unsigned int N_cols,
    unsigned int n);

__global__ void computeVerticalNorthKernel_friction(
    const unsigned int* ids,
    const double* H,
    const double* v,
    double* vertical,
    unsigned int N_cols,
    unsigned int n);

__global__ void computeVerticalSouthKernel_friction(
    const unsigned int* ids,
    const double* H,
    const double* v,
    double* vertical,
    unsigned int N_cols,
    unsigned int n);

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





