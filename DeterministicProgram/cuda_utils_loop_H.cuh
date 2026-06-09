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



double deviceMax(const thrust::device_vector<double> &v);
double deviceMin(const thrust::device_vector<double> &v);
double deviceSum(const thrust::device_vector<double> &v);

struct MinMaxResult {
  double min_val;
  double max_val;
};

MinMaxResult deviceMinMax(const thrust::device_vector<double> &v);

//==============================================================================

template <typename Kernel, typename... Args>
inline void launch_kernel(Kernel kernel, std::size_t n, cudaStream_t stream, Args&&... args) {
    constexpr int blockSize = 256;
    int numBlocks = (static_cast<int>(n) + blockSize - 1) / blockSize;
    kernel<<<numBlocks, blockSize, 0, stream>>>(std::forward<Args>(args)..., static_cast<int>(n));
}

//==============================================================================
__device__ __forceinline__
double M_gamma_dt_DSV(double dt, double coeff);

__device__ inline int signum_dev(double val);

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

void compute_horizontal_interface_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectHorizontal, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectWest, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectEast,
		const thrust::device_vector<double>& H, 
		const thrust::device_vector<double>& u, 
		thrust::device_vector<double>& horizontal, const unsigned int N_cols,
		cudaStream_t stream);

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

void compute_vertical_interface_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectVertical, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectNorth, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectSouth,
		const thrust::device_vector<double>& H, 
		const thrust::device_vector<double>& v, 
		thrust::device_vector<double>& vertical, const unsigned int N_cols,
		cudaStream_t stream);

//==============================================================================

__global__ void computeKernel_friction(
    const unsigned int* ids,
    const double* H_interface_dir,
    const double* vel,
    const double* M_expo_r_dir_vect,
    double* alfa_dir,
    const double* M_gamma_dt_DSV_dir_,
    double M_dt_DSV, double M_coeff, double M_H_min,
    double M_expo, unsigned int M_frictionModel,
    unsigned int n);

// Horizontal friction
void compute_horizontal_friction_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectHorizontal, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectWest, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectEast,
		const thrust::device_vector<double>& H_interface_horizontal, 
		const thrust::device_vector<double>& u, const thrust::device_vector<double>& M_expo_r_x_vect, 
		      thrust::device_vector<double>& alfa_x, 
		const thrust::device_vector<double>& M_gamma_dt_DSV_x_, 
		double M_dt_DSV, double M_coeff, double M_H_min, 
                double M_expo, unsigned int M_frictionModel, cudaStream_t stream); 

// Vertical friction
void compute_vertical_friction_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectVertical, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectNorth, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectSouth,
		const thrust::device_vector<double>& H_interface_vertical, 
		const thrust::device_vector<double>& v, 
		const thrust::device_vector<double>& M_expo_r_y_vect, 
		      thrust::device_vector<double>& alfa_y, 
		const thrust::device_vector<double>& M_gamma_dt_DSV_y_, 
		double M_dt_DSV, double M_coeff, double M_H_min, 
                double M_expo, unsigned int M_frictionModel, cudaStream_t stream); 

//==============================================================================

__global__ void computeTemperatureKernel(
    const unsigned int* idBasinVect, 
    double* T_raster, double* melt_mask,
    const double* orography,
    const double T, const double Temp_diff,
    const double height_th, const double T_crit,
    unsigned int n);

void computeTemperature_wrapper(
    const thrust::device_vector<unsigned int> &idBasinVect, 
    thrust::device_vector<double>& T_raster, thrust::device_vector<double>& melt_mask,
    const thrust::device_vector<double>& orography,
    const double T, const double Temp_diff,
    const double height_th, const double T_crit, cudaStream_t stream);

//==============================================================================

__global__ void computeETKernel(
    const unsigned int* idBasinVect,
    double* ET_vec,
    const double* orography,
    const double Ra,
    const double t_mean_base,
    const double t_max_base,
    const double t_min_base,
    const double Temp_diff,
    const double height_th,
    const unsigned int ET_model,
    const unsigned int n);

void computeET_wrapper(
    const thrust::device_vector<unsigned int> &idBasinVect,
    thrust::device_vector<double> &ET_vec,
    const thrust::device_vector<double> &orography,
    const double Ra,
    const double t_mean_base,
    const double t_max_base,
    const double t_min_base,
    const double Temp_diff,
    const double height_th,
    const unsigned int ET_model,
    cudaStream_t stream);

//==============================================================================

__global__ void computeKernel_sediment(
    const unsigned int* ids,
    double* Gamma_dir_1,
    double* Gamma_dir_2,
    const double* S_dir_mod,
    const double* vel,
    const double c1,
    unsigned int n);

// Horizontal sediment
void computeResidualsTruncatedHorizontal_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectHorizontal, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectWest, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectEast,
		thrust::device_vector<double>& Gamma_x_1, thrust::device_vector<double>& Gamma_x_2, 
		const thrust::device_vector<double>& S_x_mod, 
		const thrust::device_vector<double>& u, 
                const double c1, cudaStream_t stream); 

// Vertical sediment
void computeResidualsTruncatedVertical_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectVertical, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectNorth, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectSouth,
                thrust::device_vector<double>& Gamma_y_1, thrust::device_vector<double>& Gamma_y_2,
		const thrust::device_vector<double>& S_y_mod,
		const thrust::device_vector<double>& v,
		const double c1, cudaStream_t stream);

//==============================================================================

__global__ void computeKernelHorizontalInternal_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_x,
    const double* h,
    double* h_interface_x,
    const unsigned int N_cols,
    unsigned int n);

__global__ void computeKernelHorizontalWest_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_x,
    const double* h,
    double* h_interface_x,
    const unsigned int N_cols,
    unsigned int n);

__global__ void computeKernelHorizontalEast_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_x,
    const double* h,
    double* h_interface_x,
    const unsigned int N_cols,
    unsigned int n);

// Horizontal gravitational layer
void computeResidualsHorizontal_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectHorizontal, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectWest, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectEast,
		const thrust::device_vector<double>& coeff,
		const thrust::device_vector<double>& n_x,
		const thrust::device_vector<double>& h, 
		thrust::device_vector<double>& h_interface_x, 
                const unsigned int N_cols, cudaStream_t stream);

__global__ void computeKernelVerticalInternal_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_y,
    const double* h,
    double* h_interface_y,
    const unsigned int N_cols,
    unsigned int n);

__global__ void computeKernelVerticalNorth_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_y,
    const double* h,
    double* h_interface_y,
    const unsigned int N_cols,
    unsigned int n);

__global__ void computeKernelVerticalSouth_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_y,
    const double* h,
    double* h_interface_y,
    const unsigned int N_cols,
    unsigned int n);

// Vertical gravitational layer
void computeResidualsVertical_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectVertical, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectNorth, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectSouth,
		const thrust::device_vector<double>& coeff,
		const thrust::device_vector<double>& n_y,
		const thrust::device_vector<double>& h, 
		thrust::device_vector<double>& h_interface_y, 
                const unsigned int N_cols, cudaStream_t stream);

//==============================================================================

__global__ void computeupdateSnowGravLayers(
    const unsigned int* ids,
    const double* T_raster,
    const double* melt_mask,
    double* h_sn,
    double* h_G,
    const double* ET_vec,
    const double* DP_infiltrated,
    const double* DP_total,
    const double c1_min,
    const double dt_min,
    const double T_thr,
    const unsigned int N_cols,
    const double* h_interface_x,
    const double* h_interface_y,
    unsigned int n);

// Snow and Gravitational layer updates
void updateSnowGravLayers_wrapper(
        const thrust::device_vector<unsigned int>& idBasinVect, 
	const thrust::device_vector<double>& T_raster, 
	const thrust::device_vector<double>& melt_mask, 
	thrust::device_vector<double>& h_sn, 
	thrust::device_vector<double>& h_G, 
	const thrust::device_vector<double>& ET_vec, 
	const thrust::device_vector<double>& DP_infiltrated, 
	const thrust::device_vector<double>& DP_total, 
	const double c1_min, 
	const double dt_min, 
	const double T_thr,
	const unsigned int N_cols, 
	const thrust::device_vector<double>& h_interface_x,
	const thrust::device_vector<double>& h_interface_y, cudaStream_t stream);

//==============================================================================

__global__ void computePrecipitation(
    const unsigned int* ids,
    const double time,
    const double c,
    const bool M_isInitialLoss,
    const double* M_time_spacing_vect,
    const double* S,
    const double* melt_mask,
    const double* Hyetograph,
    const double* IDW_weights,
    double* DP_total,
    double* DP_cumulative,
    double* DP_infiltrated,
    const double* h_G,
    const double* H,
    const unsigned int N_cols,
    const unsigned int* offset_hy,
    const unsigned int Hyetograph_size,
    unsigned int n);

// Compute the precipitation in the computational domain
void computePrecipitation_wrapper(
        const thrust::device_vector<unsigned int>& idBasinVect, 
        const double time,
	const double c, const bool M_isInitialLoss,
	const thrust::device_vector<double>& M_time_spacing_vect,
        const thrust::device_vector<double>& S,
	const thrust::device_vector<double>& melt_mask,
	const thrust::device_vector<double>& Hyetograph,
        const thrust::device_vector<double>& IDW_weights,
	thrust::device_vector<double>& DP_total, 
	thrust::device_vector<double>& DP_cumulative, 
	thrust::device_vector<double>& DP_infiltrated,
	const thrust::device_vector<double>& h_G, 
	const thrust::device_vector<double>& H, 
	const unsigned int N_cols, 
	const thrust::device_vector<unsigned int>& offset_hy,
	const unsigned int Hyetograph_size,
	cudaStream_t stream);

//==============================================================================

__global__ void bilinearInterpolationHorizontal(
		const unsigned int* __restrict__ ids,
		const double* __restrict__ u,
		const double* __restrict__ v,
		double* __restrict__ u_star,
		const double scale,
		const unsigned int nrows,
		const unsigned int ncols,
		unsigned int n);

__global__ void bilinearInterpolationHorizontalWest(
		const unsigned int* ids,
                const double* u,
		const double* v,
                double* u_star,
		const double scale,
		const unsigned int N_rows, const unsigned int N_cols,
		unsigned int n);

__global__ void bilinearInterpolationHorizontalEast(
		const unsigned int* ids,
                const double* u,
		const double* v,
                double* u_star,
		const double scale,
		const unsigned int N_rows, const unsigned int N_cols,
		unsigned int n);


void bilinearInterpolationHorizontal_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectHorizontal,
                const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectWest, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectEast,
                const thrust::device_vector<double>& u, 
		const thrust::device_vector<double>& v, 
		thrust::device_vector<double>& u_star, 
		const double scale, 
		const unsigned int nrows, const unsigned int ncols, cudaStream_t stream);


__global__ void bilinearInterpolationVertical(
		const unsigned int* __restrict__ ids,
		const double* __restrict__ u,
		const double* __restrict__ v,
		double* __restrict__ v_star,
		const double scale,
		const unsigned int nrows,
		const unsigned int ncols,
		unsigned int n);

__global__ void bilinearInterpolationVerticalNorth(
		const unsigned int* ids,
                const double* u,
		const double* v,
                double* v_star,
		const double scale,
		const unsigned int N_rows, const unsigned int N_cols,
		unsigned int n);

__global__ void bilinearInterpolationVerticalSouth(
		const unsigned int* ids,
                const double* u,
		const double* v,
                double* v_star,
		const double scale,
		const unsigned int N_rows, const unsigned int N_cols,
		unsigned int n);


void bilinearInterpolationVertical_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectVertical,
                const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectNorth, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectSouth,
                const thrust::device_vector<double>& u, 
		const thrust::device_vector<double>& v, 
		thrust::device_vector<double>& v_star, 
		const double scale,
		const unsigned int nrows, const unsigned int ncols, cudaStream_t stream);

//==============================================================================

typedef struct VecStruct {
    cusparseDnVecDescr_t vec;
    double*              ptr;
} Vec;

//==============================================================================

__device__ int findPosition(const int* d_A_rows, const int* d_A_columns,
                             int row, int col);


__global__ void buildMatrix_cell_center(
              const unsigned int* ids,
              const double* H,
	      const double* precipitation,
	      const unsigned int* idBasinVectReIndex,
	      const double dt_DSV,
	      double* rhs,
	      const int* d_A_rows,
	      const int* d_A_columns,
	      double* d_A_values,
	      unsigned int n);


__global__ void buildMatrix_horizontal_internal(
		const unsigned int* ids,
                const double* H_int_x,
		const double* alfa_x,
                const unsigned int* idBasinVectReIndex,
		const double* orography,
		const double* u_star,
		const unsigned int N_cols,
		const double c1, 
		const double c3,
		const double H_min,
		double* rhs,
		const int* d_A_rows,
		const int* d_A_columns,
		double* d_A_values,
		unsigned int n);

__global__ void buildMatrix_horizontal_West(
		const unsigned int* ids,
                const double* H_int_x,
		const double* alfa_x,
                const unsigned int* idBasinVectReIndex,
		const double* u_star,
		const unsigned int N_cols,
		const double c1, 
		const double H_min,
		const bool isNonReflectingBC,
		double* rhs,
		const int* d_A_rows,
		const int* d_A_columns,
		double* d_A_values,
		unsigned int n);
 
__global__ void buildMatrix_horizontal_East(
		const unsigned int* ids,
                const double* H_int_x,
		const double* alfa_x,
                const unsigned int* idBasinVectReIndex,
		const double* u_star,
		const unsigned int N_cols,
		const double c1, 
		const double H_min,
		const bool isNonReflectingBC,
		double* rhs,
		const int* d_A_rows,
		const int* d_A_columns,
		double* d_A_values,
		unsigned int n);



__global__ void buildMatrix_vertical_internal(
		const unsigned int* ids,
                const double* H_int_y,
		const double* alfa_y,
                const unsigned int* idBasinVectReIndex,
		const double* orography,
		const double* v_star,
		const unsigned int N_cols,
		const double c1, 
		const double c3,
		const double H_min,
		double* rhs,
		const int* d_A_rows,
		const int* d_A_columns,
		double* d_A_values,
		unsigned int n);

__global__ void buildMatrix_vertical_North(
		const unsigned int* ids,
                const double* H_int_y,
		const double* alfa_y,
                const unsigned int* idBasinVectReIndex,
		const double* v_star,
		const unsigned int N_cols,
		const double c1, 
		const double H_min,
		const bool isNonReflectingBC,
		double* rhs,
		const int* d_A_rows,
		const int* d_A_columns,
		double* d_A_values,
		unsigned int n);

__global__ void buildMatrix_vertical_South(
		const unsigned int* ids,
                const double* H_int_y,
		const double* alfa_y,
                const unsigned int* idBasinVectReIndex,
		const double* v_star,
		const unsigned int N_cols,
		const double c1, 
		const double H_min,
		const bool isNonReflectingBC,
		double* rhs,
		const int* d_A_rows,
		const int* d_A_columns,
		double* d_A_values,
		unsigned int n);

void buildMatrix_wrapper(const thrust::device_vector<double>& H_int_x, 
		const thrust::device_vector<double>& H_int_y, 
		const thrust::device_vector<double>& orography, 
		const thrust::device_vector<double>& u_star, 
		const thrust::device_vector<double>& v_star, 
		const thrust::device_vector<double>& u, 
		const thrust::device_vector<double>& v, 
		const thrust::device_vector<double>& H, 
		const unsigned int N_cols, 
		const double c1, 
		const double c3, 
		const double H_min, 
		const thrust::device_vector<double>& precipitation,
		const double dt_DSV, 
		const thrust::device_vector<double>& alfa_x, 
		const thrust::device_vector<double>& alfa_y, 
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectHorizontal, 
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectVertical, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectWest, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectEast, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectNorth, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectSouth, 
		const thrust::device_vector<unsigned int>& idBasinVect,
		const thrust::device_vector<unsigned int>& idBasinVectReIndex,
		const bool isNonReflectingBC, 
		const int nnz,
		const int* d_A_rows, 
		const int* d_A_columns,
	       	double* d_A_values, 
		Vec& rhs, cudaStream_t stream);

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
