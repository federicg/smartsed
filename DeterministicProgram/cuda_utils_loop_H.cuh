#ifndef CUDA_UTILS_LOOP_H_CUH
#define CUDA_UTILS_LOOP_H_CUH

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


double deviceMax(const double* raw_ptr, size_t n);
double deviceMax(const thrust::device_vector<double> &v);
double deviceMin(const double* raw_ptr, size_t n);
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
// Merged internal+West+East (0/1/2 tag), single launch.
__global__ void computeHorizontalMergedKernel_interface(
    const unsigned int* ids,
    const unsigned int* tag,
    const double* H,
    const double* u,
    double* horizontal,
    unsigned int N_cols,
    unsigned int n);

void compute_horizontal_interface_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		const thrust::device_vector<unsigned int>& tag,
		const thrust::device_vector<double>& H,
		const thrust::device_vector<double>& u,
		thrust::device_vector<double>& horizontal, const unsigned int N_cols,
		cudaStream_t stream);

// Vertical upwind, merged internal+North+South (0/1/2 tag).
__global__ void computeVerticalMergedKernel_interface(
    const unsigned int* ids,
    const unsigned int* tag,
    const double* H,
    const double* v,
    double* vertical,
    unsigned int N_cols,
    unsigned int n);

void compute_vertical_interface_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		const thrust::device_vector<unsigned int>& tag,
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

// Merged internal+boundary friction: single launch over the concatenated id
// list, used for both directions (computeKernel_friction's formula doesn't
// distinguish face type).
void compute_friction_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		const thrust::device_vector<double>& H_interface_dir,
		const thrust::device_vector<double>& vel,
		const thrust::device_vector<double>& M_expo_r_dir_vect,
		      thrust::device_vector<double>& alfa_dir,
		const thrust::device_vector<double>& M_gamma_dt_DSV_dir_,
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

// Sediment residuals, one direction: single launch over the concatenated
// internal+boundary id list (computeKernel_sediment's formula doesn't
// distinguish face type).
void computeResidualsTruncated_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		thrust::device_vector<double>& Gamma_dir_1, thrust::device_vector<double>& Gamma_dir_2,
		const thrust::device_vector<double>& S_dir_mod,
		const thrust::device_vector<double>& vel,
                const double c1, cudaStream_t stream);

// Gravitational-layer sediment production W_Gav (base EPM term only; the CPU
// per-pour-point additional_source_term from the static excluded-subbasin
// approximation is not carried on the GPU). pi_1e3 = 1.e-3 * M_PI folded on
// the host to keep the multiply order identical to the CPU reference.
__global__ void computeKernel_WGav(
    const unsigned int* ids,
    const double* Z_Gav,
    const double* T_raster,
    const double* melt_mask,
    const double* DP_total,
    double* W_Gav,
    const double pi_1e3,
    const double dt_sed,
    unsigned int n);

void computeWGav_wrapper(
    const thrust::device_vector<unsigned int>& idBasinVect,
    const thrust::device_vector<double>& Z_Gav,
    const thrust::device_vector<double>& T_raster,
    const thrust::device_vector<double>& melt_mask,
    const thrust::device_vector<double>& DP_total,
          thrust::device_vector<double>& W_Gav,
    const double pi_1e3, const double dt_sed, cudaStream_t stream);

//==============================================================================

// Gravitational layer, merged internal+West+East / internal+North+South
// (0/1/2 tag), single launch per direction.
__global__ void computeHorizontalMergedKernel_gravitational(
    const unsigned int* ids,
    const unsigned int* tag,
    const double* coeff,
    const double* n_x,
    const double* h,
    double* h_interface_x,
    const unsigned int N_cols,
    unsigned int n);

// Horizontal gravitational layer
void computeResidualsHorizontal_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		const thrust::device_vector<unsigned int>& tag,
		const thrust::device_vector<double>& coeff,
		const thrust::device_vector<double>& n_x,
		const thrust::device_vector<double>& h,
		thrust::device_vector<double>& h_interface_x,
                const unsigned int N_cols, cudaStream_t stream);

__global__ void computeVerticalMergedKernel_gravitational(
    const unsigned int* ids,
    const unsigned int* tag,
    const double* coeff,
    const double* n_y,
    const double* h,
    double* h_interface_y,
    const unsigned int N_cols,
    unsigned int n);

// Vertical gravitational layer
void computeResidualsVertical_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		const thrust::device_vector<unsigned int>& tag,
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

// Bilinear back-trace, split as internal (tag==0) then merged boundary
// (tag 1/2). Both kernels are launched over the whole concatenated
// internal+West+East / internal+North+South list; the boundary pass must
// stay a separate launch because it reads u_star/v_star the internal pass writes.
__global__ void bilinearInterpolationHorizontal(
		const unsigned int* __restrict__ ids,
		const unsigned int* __restrict__ tag,
		const double* __restrict__ u,
		const double* __restrict__ v,
		double* __restrict__ u_star,
		const double scale,
		const unsigned int nrows,
		const unsigned int ncols,
		unsigned int n);

__global__ void bilinearInterpolationHorizontalBoundary(
		const unsigned int* ids,
                const unsigned int* tag,
                const double* u,
		const double* v,
                double* u_star,
		const double scale,
		const unsigned int N_rows, const unsigned int N_cols,
		unsigned int n);


void bilinearInterpolationHorizontal_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
                const thrust::device_vector<unsigned int>& tag,
                const thrust::device_vector<double>& u,
		const thrust::device_vector<double>& v,
		thrust::device_vector<double>& u_star,
		const double scale,
		const unsigned int nrows, const unsigned int ncols, cudaStream_t stream);


__global__ void bilinearInterpolationVertical(
		const unsigned int* __restrict__ ids,
		const unsigned int* __restrict__ tag,
		const double* __restrict__ u,
		const double* __restrict__ v,
		double* __restrict__ v_star,
		const double scale,
		const unsigned int nrows,
		const unsigned int ncols,
		unsigned int n);

__global__ void bilinearInterpolationVerticalBoundary(
		const unsigned int* ids,
                const unsigned int* tag,
                const double* u,
		const double* v,
                double* v_star,
		const double scale,
		const unsigned int N_rows, const unsigned int N_cols,
		unsigned int n);


void bilinearInterpolationVertical_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
                const thrust::device_vector<unsigned int>& tag,
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


// Cell-centered gather replacing the old face-centered scatter kernels
// (buildMatrix_horizontal_internal + buildMatrix_vertical_internal). See the
// .cu definition for why this needs no atomics.
__global__ void buildMatrix_cell_gather(
		const unsigned int* ids,
                const double* H_int_x,
                const double* H_int_y,
		const double* alfa_x,
		const double* alfa_y,
		const double* orography,
		const double* u_star,
		const double* v_star,
		const unsigned int* basin_mask,
                const unsigned int* idBasinVectReIndex,
		const unsigned int N_cols,
		const unsigned int N,
		const double c1,
		const double c3,
		const double H_min,
		double* rhs,
		const int* d_A_rows,
		const int* d_A_columns,
		double* d_A_values,
		unsigned int n);

// Merged boundary rhs contribution; launched once per direction (horizontal
// list + _x arrays + isHorizontal=1, vertical list + _y arrays + isHorizontal=0).
// Replaces buildMatrix_{horizontal,vertical}_{West,East,North,South}.
__global__ void buildMatrix_boundary_merged(
		const unsigned int* ids,
                const unsigned int* tag,
                const double* H_int,
		const double* alfa,
		const double* vel_star,
                const unsigned int* idBasinVectReIndex,
		const unsigned int N_cols,
		const double c1,
		const double H_min,
		const bool isNonReflectingBC,
		const int isHorizontal,
		double* rhs,
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
		const thrust::device_vector<unsigned int>& idHorizontalAll,
		const thrust::device_vector<unsigned int>& horizontalTag,
		const thrust::device_vector<unsigned int>& idVerticalAll,
		const thrust::device_vector<unsigned int>& verticalTag,
		const thrust::device_vector<unsigned int>& idBasinVect,
		const thrust::device_vector<unsigned int>& idBasinVectReIndex,
		const thrust::device_vector<unsigned int>& basin_mask,
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
	   cusparseSpSVDescr_t  spsvDescrL,
           cusparseSpSVDescr_t  spsvDescrLT,
           int                  maxIterations,
           double               tolerance,
	   const bool           use_preconditioner);

//==============================================================================
// Scatter CG solution d_X back into H and update eta = H + orography
__global__ void updateH_kernel(
    const unsigned int* ids,
    const unsigned int* idBasinVectReIndex,
    double* H,
    double* eta,
    const double* orography,
    const double* d_X,
    unsigned int n);

void updateH_wrapper(
    const thrust::device_vector<unsigned int>& idBasinVect,
    const thrust::device_vector<unsigned int>& idBasinVectReIndex,
    thrust::device_vector<double>& H,
    thrust::device_vector<double>& eta,
    const thrust::device_vector<double>& orography,
    const double* d_X,
    cudaStream_t stream);

//==============================================================================
// updateVel: merged internal+boundary kernels (2 launches instead of 6; see
// the tag-construction comment in main_final_H.cpp and the kernel comments
// in cuda_utils_loop_H.cu for why membership is tagged at setup rather than
// re-derived per face).
__global__ void updateVelHorizontalMerged(
    const unsigned int* ids,
    const unsigned int* tag,
    double* u,
    const double* u_star,
    const double* alfa_x,
    const double* H,
    const double* eta,
    double c2, double H_min,
    int isNonReflectingBC,
    unsigned int N_cols,
    unsigned int n);

__global__ void updateVelVerticalMerged(
    const unsigned int* ids,
    const unsigned int* tag,
    double* v,
    const double* v_star,
    const double* alfa_y,
    const double* H,
    const double* eta,
    double c2, double H_min,
    int isNonReflectingBC,
    unsigned int N_cols,
    unsigned int n);

void updateVel_wrapper(
    thrust::device_vector<double>& u,
    thrust::device_vector<double>& v,
    const thrust::device_vector<double>& u_star,
    const thrust::device_vector<double>& v_star,
    const thrust::device_vector<double>& alfa_x,
    const thrust::device_vector<double>& alfa_y,
    unsigned int N_rows, unsigned int N_cols,
    double c2, double H_min,
    const thrust::device_vector<double>& eta,
    const thrust::device_vector<double>& H,
    const thrust::device_vector<double>& orography,
    const thrust::device_vector<unsigned int>& idHorizontalAll,
    const thrust::device_vector<unsigned int>& horizontalTag,
    const thrust::device_vector<unsigned int>& idVerticalAll,
    const thrust::device_vector<unsigned int>& verticalTag,
    int isNonReflectingBC,
    cudaStream_t stream);

//==============================================================================
#endif
