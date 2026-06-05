#include <thrust/extrema.h>      // for thrust::max_element, min_element, minmax_element
#include <thrust/reduce.h>       // for thrust::reduce
#include <thrust/device_vector.h>
#include "cuda_utils_loop_H.cuh"


double deviceMax(const thrust::device_vector<double> &v) {
  return *thrust::max_element(v.begin(), v.end());
}

double deviceMin(const thrust::device_vector<double> &v) {
  return *thrust::min_element(v.begin(), v.end());
}

double deviceSum(const thrust::device_vector<double> &v) {
  return thrust::reduce(v.begin(), v.end(), 0.0, thrust::plus<double>());
}

MinMaxResult deviceMinMax(const thrust::device_vector<double> &v) {
  auto mm = thrust::minmax_element(v.begin(), v.end());
  return { static_cast<double>(*mm.first), static_cast<double>(*mm.second) };
}

//==============================================================================

__device__ __forceinline__
double M_gamma_dt_DSV(double dt, double coeff) { return dt * coeff; }

__device__ inline int signum_dev(double val) {
    return (0.0 < val) - (val < 0.0);
}

//==============================================================================
// Horizontal upwind
__global__ void computeHorizontalInternalKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* u,
    double* horizontal,
    unsigned int N_cols,
    unsigned int n)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    unsigned int Id = ids[i];
    unsigned int ii = Id / (N_cols + 1);

    unsigned int IDeast = Id - ii;
    unsigned int IDwest = Id - ii - 1;

    double H_left  = H[IDwest];
    double H_right = H[IDeast];

    double s = (u[Id] > 0.0) - (u[Id] < 0.0);

    horizontal[Id] =
        0.5 * (H_left + H_right + s * (H_left - H_right));
}

__global__ void computeHorizontalWestKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* u,
    double* horizontal,
    unsigned int N_cols,
    unsigned int n)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    unsigned int Id = ids[i];
    unsigned int ii = Id / (N_cols + 1);

    double H_left  = 0.0;
    double H_right = H[Id - ii];

    double s = (u[Id + 1] > 0.0) - (u[Id + 1] < 0.0);

    horizontal[Id] =
        0.5 * (H_left + H_right + s * (H_left - H_right));
}

__global__ void computeHorizontalEastKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* u,
    double* horizontal,
    unsigned int N_cols,
    unsigned int n)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    unsigned int Id = ids[i];
    unsigned int ii = Id / (N_cols + 1);

    double H_left  = H[Id - ii - 1];
    double H_right = 0.0;

    double s = (u[Id - 1] > 0.0) - (u[Id - 1] < 0.0);

    horizontal[Id] =
        0.5 * (H_left + H_right + s * (H_left - H_right));
}

void compute_horizontal_interface_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectHorizontal, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectWest, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectEast,
		const thrust::device_vector<double>& H, 
		const thrust::device_vector<double>& u, 
		thrust::device_vector<double>& horizontal, const unsigned int N_cols,
		cudaStream_t stream) {

    launch_kernel(
      computeHorizontalInternalKernel_interface,
      idStaggeredInternalVectHorizontal.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredInternalVectHorizontal.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(horizontal.data()),
      N_cols
    );

    launch_kernel(
      computeHorizontalWestKernel_interface,
      idStaggeredBoundaryVectWest.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectWest.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(horizontal.data()),
      N_cols
    );

    launch_kernel(
      computeHorizontalEastKernel_interface,
      idStaggeredBoundaryVectEast.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectEast.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(horizontal.data()),
      N_cols
    );

}

// Vertical upwind
__global__ void computeVerticalInternalKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* v,
    double* vertical,
    unsigned int N_cols,
    unsigned int n)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  unsigned int Id = ids[i];
  unsigned int IDsouth = Id, IDnorth = Id - N_cols;
  
  double H_left = H[IDnorth]; 
  double H_right = H[IDsouth];

  double v_cc = v[Id];
  
  double s = (v_cc > 0.0) - (v_cc < 0.0);

  vertical[Id] =
          (H_left + H_right) * .5 + s * (H_left - H_right) * .5;

}

__global__ void computeVerticalNorthKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* v,
    double* vertical,
    unsigned int N_cols,
    unsigned int n)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  unsigned int Id = ids[i];

 
  double H_left = 0;
  double H_right = H[Id];

  double v_right = v[Id + N_cols];

  double s = (v_right > 0.0) - (v_right < 0.0);

  vertical[Id] = (H_left + H_right) * .5 +
                     s * (H_left - H_right) * .5;

}

__global__ void computeVerticalSouthKernel_interface(
    const unsigned int* ids,
    const double* H,
    const double* v,
    double* vertical,
    unsigned int N_cols,
    unsigned int n)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
    
  unsigned int Id = ids[i];
  
  double H_left = H[Id - N_cols]; 
  double H_right = 0;

  double v_left = v[Id - N_cols];
  
  double s = (v_left > 0.0) - (v_left < 0.0);

  vertical[Id] = (H_left + H_right) * .5 +
                     s * (H_left - H_right) * .5;

}

void compute_vertical_interface_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectVertical, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectNorth, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectSouth,
		const thrust::device_vector<double>& H, 
		const thrust::device_vector<double>& v, 
		thrust::device_vector<double>& vertical, const unsigned int N_cols,
		cudaStream_t stream) {

    launch_kernel(
      computeVerticalInternalKernel_interface,
      idStaggeredInternalVectVertical.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredInternalVectVertical.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(vertical.data()),
      N_cols
    );

    launch_kernel(
      computeVerticalNorthKernel_interface,
      idStaggeredBoundaryVectNorth.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectNorth.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(vertical.data()),
      N_cols
    );

    launch_kernel(
      computeVerticalSouthKernel_interface,
      idStaggeredBoundaryVectSouth.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectSouth.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(vertical.data()),
      N_cols
    );

}

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
    unsigned int n)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const unsigned int Id = ids[i];

  const float H_int = fmaxf((float)H_interface_dir[Id], (float)M_H_min);
  const float exponent = (float)M_expo_r_dir_vect[Id];

  const float expo1 = (float)M_expo + exponent * (M_frictionModel == 2);
  const float den = __expf(expo1 * __logf(H_int));
  double alfa = 1.0;

  if (den > (float)M_H_min) {
    const float vel_abs = fmaxf(fabsf((float)vel[Id]), 1e-12f);
    const float expo2 = 1.0f - exponent * (M_frictionModel == 2);

    const float pow_vel = __expf(expo2 * __logf(vel_abs));

    double coeff =
            M_gamma_dt_DSV(M_dt_DSV, M_coeff) * vel_abs / den *
            (M_frictionModel > 0);

    coeff = fmax(
            coeff,
            M_dt_DSV * M_gamma_dt_DSV_dir_[Id] *
            pow_vel / den);

    alfa = 1.0 / (1.0 + coeff);
  }
  alfa_dir[Id] = alfa;
}

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
                double M_expo, unsigned int M_frictionModel, cudaStream_t stream) {

    launch_kernel(
      computeKernel_friction,
      idStaggeredInternalVectHorizontal.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredInternalVectHorizontal.data()),
      thrust::raw_pointer_cast(H_interface_horizontal.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(M_expo_r_x_vect.data()),
      thrust::raw_pointer_cast(alfa_x.data()),
      thrust::raw_pointer_cast(M_gamma_dt_DSV_x_.data()),
      M_dt_DSV, M_coeff, M_H_min, 
      M_expo, M_frictionModel
    );

    launch_kernel(
      computeKernel_friction,
      idStaggeredBoundaryVectWest.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectWest.data()),
      thrust::raw_pointer_cast(H_interface_horizontal.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(M_expo_r_x_vect.data()),
      thrust::raw_pointer_cast(alfa_x.data()),
      thrust::raw_pointer_cast(M_gamma_dt_DSV_x_.data()),
      M_dt_DSV, M_coeff, M_H_min, 
      M_expo, M_frictionModel
    );

    launch_kernel(
      computeKernel_friction,
      idStaggeredBoundaryVectEast.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectEast.data()),
      thrust::raw_pointer_cast(H_interface_horizontal.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(M_expo_r_x_vect.data()),
      thrust::raw_pointer_cast(alfa_x.data()),
      thrust::raw_pointer_cast(M_gamma_dt_DSV_x_.data()),
      M_dt_DSV, M_coeff, M_H_min, 
      M_expo, M_frictionModel
    );

}

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
                double M_expo, unsigned int M_frictionModel,
		cudaStream_t stream) {

    launch_kernel(
      computeKernel_friction,
      idStaggeredInternalVectVertical.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredInternalVectVertical.data()),
      thrust::raw_pointer_cast(H_interface_vertical.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(M_expo_r_y_vect.data()),
      thrust::raw_pointer_cast(alfa_y.data()),
      thrust::raw_pointer_cast(M_gamma_dt_DSV_y_.data()),
      M_dt_DSV, M_coeff, M_H_min, 
      M_expo, M_frictionModel
    );

    launch_kernel(
      computeKernel_friction,
      idStaggeredBoundaryVectNorth.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectNorth.data()),
      thrust::raw_pointer_cast(H_interface_vertical.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(M_expo_r_y_vect.data()),
      thrust::raw_pointer_cast(alfa_y.data()),
      thrust::raw_pointer_cast(M_gamma_dt_DSV_y_.data()),
      M_dt_DSV, M_coeff, M_H_min, 
      M_expo, M_frictionModel
    );

    launch_kernel(
      computeKernel_friction,
      idStaggeredBoundaryVectSouth.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectSouth.data()),
      thrust::raw_pointer_cast(H_interface_vertical.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(M_expo_r_y_vect.data()),
      thrust::raw_pointer_cast(alfa_y.data()),
      thrust::raw_pointer_cast(M_gamma_dt_DSV_y_.data()),
      M_dt_DSV, M_coeff, M_H_min, 
      M_expo, M_frictionModel
    );

}

//==============================================================================

__global__ void computeTemperatureKernel(
    const unsigned int* idBasinVect, 
    double* T_raster, double* melt_mask,
    const double* orography,
    const double T, const double Temp_diff,
    const double height_th, const double T_crit,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const unsigned int j = idBasinVect[i];
  T_raster[j] = T + Temp_diff * (orography[j] - height_th);
  melt_mask[j] = (T_raster[j] > T_crit) ? 1.0 : 0.0;
}

void computeTemperature_wrapper(
    const thrust::device_vector<unsigned int> &idBasinVect,
    thrust::device_vector<double> &T_raster,
    thrust::device_vector<double> &melt_mask,
    const thrust::device_vector<double> &orography,
    const double T, const double Temp_diff,
    const double height_th, const double T_crit,
    cudaStream_t stream) {

  launch_kernel(computeTemperatureKernel, idBasinVect.size(), stream,
    thrust::raw_pointer_cast(idBasinVect.data()),
    thrust::raw_pointer_cast(T_raster.data()),
    thrust::raw_pointer_cast(melt_mask.data()),
    thrust::raw_pointer_cast(orography.data()),
    T, Temp_diff, height_th, T_crit);
}

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
    const unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const unsigned int k = idBasinVect[i];

  const double t_mean = t_mean_base + Temp_diff * (orography[k] - height_th);
  const double t_max  = t_max_base  + Temp_diff * (orography[k] - height_th);
  const double t_min  = t_min_base  + Temp_diff * (orography[k] - height_th);

  ET_vec[k] = .0023 * Ra * (t_mean + 17.8) *
              sqrt(t_max - t_min) * (1.e-3 / (24 * 3600));
}

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
    cudaStream_t stream) {

  if (ET_model == 0) return;

  launch_kernel(computeETKernel,
    idBasinVect.size(),
    stream,
    thrust::raw_pointer_cast(idBasinVect.data()),
    thrust::raw_pointer_cast(ET_vec.data()),
    thrust::raw_pointer_cast(orography.data()),
    Ra, t_mean_base, t_max_base, t_min_base,
    Temp_diff, height_th, ET_model);
}

//==============================================================================

__global__ void computeKernelHorizontalInternal_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_x,
    const double* h,
    double* h_interface_x,
    const unsigned int N_cols,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
 
  const unsigned int ii = Id / (N_cols + 1), // u
	IDeast = Id - ii,                      // H
	IDwest = Id - ii - 1;                  // H

  const double &h_left = h[IDwest], &h_right = h[IDeast];

  const double k_c_left = coeff[IDwest], k_c_right = coeff[IDeast];


  h_interface_x[Id] = n_x[Id] *
	  ((k_c_left * h_left + k_c_right * h_right) +
	   n_x[Id] * (k_c_left * h_left - k_c_right * h_right)) *
	  .5;
}

__global__ void computeKernelHorizontalWest_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_x,
    const double* h,
    double* h_interface_x,
    const unsigned int N_cols,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
 
  const unsigned int ii = Id / (N_cols + 1);
  const double h_left = 0, h_right = h[Id - ii];
  const double k_c_left = 0., k_c_right = coeff[Id - ii];

  h_interface_x[Id] = n_x[Id] *
	  ((k_c_left * h_left + k_c_right * h_right) +
	   n_x[Id] * (k_c_left * h_left - k_c_right * h_right)) *
	  .5;
}

__global__ void computeKernelHorizontalEast_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_x,
    const double* h,
    double* h_interface_x,
    const unsigned int N_cols,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
 
  const unsigned int ii = Id / (N_cols + 1);
  const double h_left = h[Id - ii - 1], h_right = 0;
  const double k_c_left = coeff[Id - ii - 1], k_c_right = 0.;

  h_interface_x[Id] = n_x[Id] *
	  ((k_c_left * h_left + k_c_right * h_right) +
	   n_x[Id] * (k_c_left * h_left - k_c_right * h_right)) *
	  .5;
}


// Horizontal gravitational layer
void computeResidualsHorizontal_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectHorizontal, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectWest, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectEast,
		const thrust::device_vector<double>& coeff,
		const thrust::device_vector<double>& n_x,
		const thrust::device_vector<double>& h, 
		thrust::device_vector<double>& h_interface_x, 
                const unsigned int N_cols, cudaStream_t stream) {

    launch_kernel(
      computeKernelHorizontalInternal_gravitational,
      idStaggeredInternalVectHorizontal.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredInternalVectHorizontal.data()),
      thrust::raw_pointer_cast(coeff.data()),
      thrust::raw_pointer_cast(n_x.data()),
      thrust::raw_pointer_cast(h.data()),
      thrust::raw_pointer_cast(h_interface_x.data()),
      N_cols
    );

    launch_kernel(
      computeKernelHorizontalWest_gravitational,
      idStaggeredBoundaryVectWest.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectWest.data()),
      thrust::raw_pointer_cast(coeff.data()),
      thrust::raw_pointer_cast(n_x.data()),
      thrust::raw_pointer_cast(h.data()),
      thrust::raw_pointer_cast(h_interface_x.data()),
      N_cols
    );

    launch_kernel(
      computeKernelHorizontalEast_gravitational,
      idStaggeredBoundaryVectEast.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectEast.data()),
      thrust::raw_pointer_cast(coeff.data()),
      thrust::raw_pointer_cast(n_x.data()),
      thrust::raw_pointer_cast(h.data()),
      thrust::raw_pointer_cast(h_interface_x.data()),
      N_cols
    );

}

__global__ void computeKernelVerticalInternal_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_y,
    const double* h,
    double* h_interface_y,
    const unsigned int N_cols,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
 
  const unsigned int IDsouth = Id, // H
	IDnorth = Id - N_cols;       // H

  const double h_left = h[IDnorth], h_right = h[IDsouth];

  const double k_c_left = coeff[IDnorth], k_c_right = coeff[IDsouth];

  h_interface_y[Id] = n_y[Id] *
	  ((k_c_left * h_left + k_c_right * h_right) +
	   n_y[Id] * (k_c_left * h_left - k_c_right * h_right)) *
	  .5;
}

__global__ void computeKernelVerticalNorth_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_y,
    const double* h,
    double* h_interface_y,
    const unsigned int N_cols,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
  const double h_left = 0, h_right = h[Id];

  const double k_c_left = 0., k_c_right = coeff[Id];

  h_interface_y[Id] = n_y[Id] *
	  ((k_c_left * h_left + k_c_right * h_right) +
	   n_y[Id] * (k_c_left * h_left - k_c_right * h_right)) *
	  .5;

}

__global__ void computeKernelVerticalSouth_gravitational(
    const unsigned int* ids,
    const double* coeff,
    const double* n_y,
    const double* h,
    double* h_interface_y,
    const unsigned int N_cols,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
  const double h_left = h[Id - N_cols], h_right = 0;

  const double k_c_left = coeff[Id - N_cols], k_c_right = 0.;

  h_interface_y[Id] = n_y[Id] *
	  ((k_c_left * h_left + k_c_right * h_right) +
	   n_y[Id] * (k_c_left * h_left - k_c_right * h_right)) *
	  .5;

}

// Vertical gravitational layer
void computeResidualsVertical_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectVertical, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectNorth, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectSouth,
		const thrust::device_vector<double>& coeff,
		const thrust::device_vector<double>& n_y,
		const thrust::device_vector<double>& h, 
		thrust::device_vector<double>& h_interface_y, 
                const unsigned int N_cols, cudaStream_t stream) {

    launch_kernel(
      computeKernelVerticalInternal_gravitational,
      idStaggeredInternalVectVertical.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredInternalVectVertical.data()),
      thrust::raw_pointer_cast(coeff.data()),
      thrust::raw_pointer_cast(n_y.data()),
      thrust::raw_pointer_cast(h.data()),
      thrust::raw_pointer_cast(h_interface_y.data()),
      N_cols
    );

    launch_kernel(
      computeKernelVerticalNorth_gravitational,
      idStaggeredBoundaryVectNorth.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectNorth.data()),
      thrust::raw_pointer_cast(coeff.data()),
      thrust::raw_pointer_cast(n_y.data()),
      thrust::raw_pointer_cast(h.data()),
      thrust::raw_pointer_cast(h_interface_y.data()),
      N_cols
    );

    launch_kernel(
      computeKernelVerticalSouth_gravitational,
      idStaggeredBoundaryVectSouth.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectSouth.data()),
      thrust::raw_pointer_cast(coeff.data()),
      thrust::raw_pointer_cast(n_y.data()),
      thrust::raw_pointer_cast(h.data()),
      thrust::raw_pointer_cast(h_interface_y.data()),
      N_cols
    );

}

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
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto k = ids[i];
    
  const unsigned int ii = k / N_cols;	
  const double S_coeff = 4.62e-10 * h_sn[k] * (T_raster[k] - T_thr) * melt_mask[k];
  const double Res_x_cell = h_interface_x[k + 1 + ii] - h_interface_x[k + ii];
  const double Res_y_cell = h_interface_y[k + N_cols] - h_interface_y[k];	  
  h_G[k] = fmax(h_G[k] + (S_coeff - ET_vec[k]) * dt_min + 
		  DP_infiltrated[k] * dt_min - c1_min * (Res_x_cell + Res_y_cell), 0.);
  h_sn[k] += DP_total[k] * (1. - melt_mask[k]) * dt_min - S_coeff * dt_min;
}

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
	const thrust::device_vector<double>& h_interface_y, cudaStream_t stream) {
 
     launch_kernel(
      computeupdateSnowGravLayers,
      idBasinVect.size(),
      stream,
      thrust::raw_pointer_cast(idBasinVect.data()),
      thrust::raw_pointer_cast(T_raster.data()),
      thrust::raw_pointer_cast(melt_mask.data()),
      thrust::raw_pointer_cast(h_sn.data()),
      thrust::raw_pointer_cast(h_G.data()),
      thrust::raw_pointer_cast(ET_vec.data()),
      thrust::raw_pointer_cast(DP_infiltrated.data()),
      thrust::raw_pointer_cast(DP_total.data()),
      c1_min, dt_min, T_thr, 
      N_cols, 
      thrust::raw_pointer_cast(h_interface_x.data()), 
      thrust::raw_pointer_cast(h_interface_y.data())
    );

}

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
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto IDcenter = ids[i];

  const auto H_cell = H[IDcenter];
  const auto h_G_cell = h_G[IDcenter];
  const auto S_cell = S[IDcenter];
  const auto melt_mask_cell = melt_mask[IDcenter];

  double DP_total_cell = 0, DP_cumulative_cell = 0, DP_infiltrated_cell = 0;

  const double deltaSoilMoisture = h_G_cell - S_cell;
  const bool cond_cc = (((H_cell + h_G_cell) < (c * S_cell)) && M_isInitialLoss);
	
  const double weight = (S_cell > 0 && deltaSoilMoisture < 0.) ? 
	  (deltaSoilMoisture/S_cell)*(deltaSoilMoisture/S_cell) : 0.;

  // SCS-CN method and Initial and Constant Loss Model
  for (unsigned int Id = 0; Id < Hyetograph_size; Id++) {

      const unsigned int i_index = floor(time / (M_time_spacing_vect[Id] * 3600));

      const double rainfall_intensity = Hyetograph[offset_hy[Id] + i_index] * 
	      IDW_weights[IDcenter*Hyetograph_size + Id];
  
      const double infiltrationRate = cond_cc ? rainfall_intensity*melt_mask_cell
	      : weight*rainfall_intensity*melt_mask_cell;
      const double potential_runoff = cond_cc ? 0. 
	      : fmax(rainfall_intensity*melt_mask_cell-infiltrationRate, 0.);

      DP_total_cell += rainfall_intensity;
      DP_cumulative_cell += potential_runoff;
      DP_infiltrated_cell += infiltrationRate;

  }

  DP_total[IDcenter] = DP_total_cell;
  DP_cumulative[IDcenter] = DP_cumulative_cell;
  DP_infiltrated[IDcenter] = DP_infiltrated_cell;

}

// Compute the precipitation in the computational domain
void computePrecipitation_wrapper(
        const thrust::device_vector<unsigned int>& idBasinVect,
        const double time, const double c,
        const bool M_isInitialLoss,
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
	cudaStream_t stream) {
 
     launch_kernel(
      computePrecipitation,
      idBasinVect.size(),
      stream,
      thrust::raw_pointer_cast(idBasinVect.data()),
      time, c, M_isInitialLoss,
      thrust::raw_pointer_cast(M_time_spacing_vect.data()),
      thrust::raw_pointer_cast(S.data()),
      thrust::raw_pointer_cast(melt_mask.data()),
      thrust::raw_pointer_cast(Hyetograph.data()),
      thrust::raw_pointer_cast(IDW_weights.data()),
      thrust::raw_pointer_cast(DP_total.data()),
      thrust::raw_pointer_cast(DP_cumulative.data()),
      thrust::raw_pointer_cast(DP_infiltrated.data()),
      thrust::raw_pointer_cast(h_G.data()),
      thrust::raw_pointer_cast(H.data()),
      N_cols,
      thrust::raw_pointer_cast(offset_hy.data()),
      Hyetograph_size
    );

}

//==============================================================================

__global__ void computeKernel_sediment(
    const unsigned int* ids,
    double* Gamma_dir_1,
    double* Gamma_dir_2,
    const double* S_dir_mod,
    const double* vel,
    const double c1,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];

  const auto vel_Id = vel[Id];

  const double coeff_right = c1 * S_dir_mod[Id] * vel_Id * (.5 - .5 * signum_dev(vel_Id));
  const double coeff_left  = c1 * S_dir_mod[Id] * vel_Id * (.5 + .5 * signum_dev(vel_Id));

  Gamma_dir_1[Id] = coeff_right;
  Gamma_dir_2[Id] = coeff_left;
}

// Horizontal sediment
void computeResidualsTruncatedHorizontal_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectHorizontal, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectWest, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectEast,
		thrust::device_vector<double>& Gamma_x_1, thrust::device_vector<double>& Gamma_x_2, 
		const thrust::device_vector<double>& S_x_mod, 
		const thrust::device_vector<double>& u, 
                const double c1, cudaStream_t stream) {

    launch_kernel(
      computeKernel_sediment,
      idStaggeredInternalVectHorizontal.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredInternalVectHorizontal.data()),
      thrust::raw_pointer_cast(Gamma_x_1.data()),
      thrust::raw_pointer_cast(Gamma_x_2.data()),
      thrust::raw_pointer_cast(S_x_mod.data()),
      thrust::raw_pointer_cast(u.data()),
      c1
    );

    launch_kernel(
      computeKernel_sediment,
      idStaggeredBoundaryVectWest.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectWest.data()),
      thrust::raw_pointer_cast(Gamma_x_1.data()),
      thrust::raw_pointer_cast(Gamma_x_2.data()),
      thrust::raw_pointer_cast(S_x_mod.data()),
      thrust::raw_pointer_cast(u.data()),
      c1
    );

    launch_kernel(
      computeKernel_sediment,
      idStaggeredBoundaryVectEast.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectEast.data()),
      thrust::raw_pointer_cast(Gamma_x_1.data()),
      thrust::raw_pointer_cast(Gamma_x_2.data()),
      thrust::raw_pointer_cast(S_x_mod.data()),
      thrust::raw_pointer_cast(u.data()),
      c1
    );

}

// Vertical sediment
void computeResidualsTruncatedVertical_wrapper(
		const thrust::device_vector<unsigned int>& idStaggeredInternalVectVertical, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectNorth, 
		const thrust::device_vector<unsigned int>& idStaggeredBoundaryVectSouth,
                thrust::device_vector<double>& Gamma_y_1, thrust::device_vector<double>& Gamma_y_2,
		const thrust::device_vector<double>& S_y_mod,
		const thrust::device_vector<double>& v,
		const double c1, cudaStream_t stream) {

    launch_kernel(
      computeKernel_sediment,
      idStaggeredInternalVectVertical.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredInternalVectVertical.data()),
      thrust::raw_pointer_cast(Gamma_y_1.data()),
      thrust::raw_pointer_cast(Gamma_y_2.data()),
      thrust::raw_pointer_cast(S_y_mod.data()),
      thrust::raw_pointer_cast(v.data()),
      c1
    );

    launch_kernel(
      computeKernel_sediment,
      idStaggeredBoundaryVectNorth.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectNorth.data()),
      thrust::raw_pointer_cast(Gamma_y_1.data()),
      thrust::raw_pointer_cast(Gamma_y_2.data()),
      thrust::raw_pointer_cast(S_y_mod.data()),
      thrust::raw_pointer_cast(v.data()),
      c1
    );

    launch_kernel(
      computeKernel_sediment,
      idStaggeredBoundaryVectSouth.size(),
      stream,
      thrust::raw_pointer_cast(idStaggeredBoundaryVectSouth.data()),
      thrust::raw_pointer_cast(Gamma_y_1.data()),
      thrust::raw_pointer_cast(Gamma_y_2.data()),
      thrust::raw_pointer_cast(S_y_mod.data()),
      thrust::raw_pointer_cast(v.data()),
      c1
    );

}



//==============================================================================
/// A 5-point Laplacian on a g x g grid with Dirichlet boundary conditions.
/// This code allocates. The caller must free.
void make_laplace_matrix(int * n_out,
                         int **row_offsets_out, 
                         int **columns_out, 
                         double **values_out) {
    int grid = 700; // grid resolution

    int n = grid * grid;
    *n_out = n;
    // vertices have 5 neighbors, 
    // but each vertex on the boundary loses 1. corners lose 2.
    int nnz = 5 * n - 4 * grid;

    printf("Creating 5-point time-dependent diffusion matrix.\n"
           " grid size: %d x %d\n"
           " matrix rows:   %d\n"
           " matrix cols:   %d\n"
           " nnz:         %d\n",
           grid, grid, n, n, nnz);

    int* row_offsets = *row_offsets_out = (int*)malloc((n + 1) * sizeof(int));
    int* columns     = *columns_out     = (int*)malloc(nnz * sizeof(int));
    double* values   = *values_out      = (double*)malloc(nnz * sizeof(double));
    assert(row_offsets);
    assert(columns);
    assert(values);

    // The Laplacian stencil looks like [-1;-1,4,-1;-1].
    // ICHOL doesn't work great with that stencil.
    // ICHOL is better suited when there's some more mass on the diagonal.
    double mass = 0.04;

    int it = 0; // next unused index into `columns`/`values`

#define INSERT(u,v, x)                    \
    if(0<=(u) && (u)<grid &&              \
       0<=(v) && (v)<grid)                \
    {                                     \
        columns[it] = ((u) * grid + (v)); \
        values[it] = x;                   \
        ++it;                             \
    }

    int row = 0;
    row_offsets[row] = 0;
    for (int i = 0; i < grid; ++i) {
        for (int j = 0; j < grid; ++j)
        {
            INSERT(i - 1, j    , -1.0);
            INSERT(i    , j - 1, -1.0);
            INSERT(i    , j    ,  4.0 + mass);
            INSERT(i    , j + 1, -1.0);
            INSERT(i + 1, j    , -1.0);
            row_offsets[++row] = it;
        }
    }
    assert(it == nnz);
#undef INSERT
}

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
           double               tolerance) {
    const double zero      = 0.0;
    const double one       = 1.0;
    const double minus_one = -1.0;
    //--------------------------------------------------------------------------
    // ### 1 ### R0 = b - A * X0 (using initial guess in X)
    //    (a) copy b in R0
    CHECK_CUDA( cudaMemcpy(d_R.ptr, d_B.ptr, m * sizeof(double),
                           cudaMemcpyDeviceToDevice) )
    //    (b) compute R = -A * X0 + R
    CHECK_CUSPARSE( cusparseSpMV(cusparseHandle,
                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &minus_one, matA, d_X.vec, &one, d_R.vec,
                                 CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT,
                                 d_bufferMV) )
    //--------------------------------------------------------------------------
    // ### 2 ### R_i_aux = L^-1 L^-T R_i
    size_t              bufferSizeL, bufferSizeLT;
    void*               d_bufferL, *d_bufferLT;
    cusparseSpSVDescr_t spsvDescrL, spsvDescrLT;
    //    (a) L^-1 tmp => R_i_aux    (triangular solver)
    CHECK_CUSPARSE( cusparseSpSV_createDescr(&spsvDescrL) )
    CHECK_CUSPARSE( cusparseSpSV_bufferSize(
                        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &one, matL, d_R.vec, d_tmp.vec, CUDA_R_64F,
                        CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL, &bufferSizeL) )
    CHECK_CUDA( cudaMalloc(&d_bufferL, bufferSizeL) )
    CHECK_CUSPARSE( cusparseSpSV_analysis(
                        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &one, matL, d_R.vec, d_tmp.vec, CUDA_R_64F,
                        CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL, d_bufferL) )
    CHECK_CUDA( cudaMemset(d_tmp.ptr, 0x0, m * sizeof(double)) )
    CHECK_CUSPARSE( cusparseSpSV_solve(
                        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &one, matL, d_R.vec, d_tmp.vec, CUDA_R_64F,
                        CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL) )

    //    (b) L^-T R_i => tmp    (triangular solver)
    CHECK_CUSPARSE( cusparseSpSV_createDescr(&spsvDescrLT) )
    CHECK_CUSPARSE( cusparseSpSV_bufferSize(
                        cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE,
                        &one, matL, d_tmp.vec, d_R_aux.vec, CUDA_R_64F,
                        CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrLT, &bufferSizeLT) )
    CHECK_CUDA( cudaMalloc(&d_bufferLT, bufferSizeLT) )
    CHECK_CUSPARSE( cusparseSpSV_analysis(
                        cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE,
                        &one, matL, d_tmp.vec, d_R_aux.vec, CUDA_R_64F,
                        CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrLT, d_bufferLT) )
    CHECK_CUDA( cudaMemset(d_R_aux.ptr, 0x0, m * sizeof(double)) )
    CHECK_CUSPARSE( cusparseSpSV_solve(
                        cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE,
                        &one, matL, d_tmp.vec, d_R_aux.vec, CUDA_R_64F,
                        CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrLT) )
    //--------------------------------------------------------------------------
    // ### 3 ### P0 = R0_aux
    CHECK_CUDA( cudaMemcpy(d_P.ptr, d_R_aux.ptr, m * sizeof(double),
                           cudaMemcpyDeviceToDevice) )
    //--------------------------------------------------------------------------
    // nrm_R0 = ||R||
    double nrm_R;
    CHECK_CUBLAS( cublasDnrm2(cublasHandle, m, d_R.ptr, 1, &nrm_R) )
    double threshold = tolerance * nrm_R;
    printf("  Initial Residual: Norm %e' threshold %e\n", nrm_R, threshold);
    //--------------------------------------------------------------------------
    double delta;
    CHECK_CUBLAS( cublasDdot(cublasHandle, m, d_R.ptr, 1, d_R_aux.ptr, 1, &delta) )
    //--------------------------------------------------------------------------
    // ### 4 ### repeat until convergence based on max iterations and
    //           and relative residual
    for (int i = 0; i < maxIterations; i++) {
        printf("  Iteration = %d; Error Norm = %e\n", i, nrm_R);
        //----------------------------------------------------------------------
        // ### 5 ### alpha = (R_i, R_aux_i) / (A * P_i, P_i)
        //     (a) T  = A * P_i
        CHECK_CUSPARSE( cusparseSpMV(cusparseHandle,
                                     CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                                     matA, d_P.vec, &zero, d_T.vec,
                                     CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT,
                                     d_bufferMV) )
        //     (b) denominator = (T, P_i)
        double denominator;
        CHECK_CUBLAS( cublasDdot(cublasHandle, m, d_T.ptr, 1, d_P.ptr, 1,
                                 &denominator) )
        //     (c) alpha = delta / denominator
        double alpha = delta / denominator;
        PRINT_INFO(delta)
        PRINT_INFO(denominator)
        PRINT_INFO(alpha)
        //----------------------------------------------------------------------
        // ### 6 ###  X_i+1 = X_i + alpha * P
        //    (a) X_i+1 = -alpha * T + X_i
        CHECK_CUBLAS( cublasDaxpy(cublasHandle, m, &alpha, d_P.ptr, 1,
                                  d_X.ptr, 1) )
        //----------------------------------------------------------------------
        // ### 7 ###  R_i+1 = R_i - alpha * (A * P)
        //    (a) R_i+1 = -alpha * T + R_i
        double minus_alpha = -alpha;
        CHECK_CUBLAS( cublasDaxpy(cublasHandle, m, &minus_alpha, d_T.ptr, 1,
                                  d_R.ptr, 1) )
        //----------------------------------------------------------------------
        // ### 8 ###  check ||R_i+1|| < threshold
        CHECK_CUBLAS( cublasDnrm2(cublasHandle, m, d_R.ptr, 1, &nrm_R) )
        PRINT_INFO(nrm_R)
        if (nrm_R < threshold)
            break;
        //----------------------------------------------------------------------
        // ### 9 ### R_aux_i+1 = L^-1 L^-T R_i+1
        //    (a) L^-1 R_i+1 => tmp    (triangular solver)
        CHECK_CUDA( cudaMemset(d_tmp.ptr,   0x0, m * sizeof(double)) )
        CHECK_CUDA( cudaMemset(d_R_aux.ptr, 0x0, m * sizeof(double)) )
        CHECK_CUSPARSE( cusparseSpSV_solve(cusparseHandle,
                                           CUSPARSE_OPERATION_NON_TRANSPOSE,
                                           &one, matL, d_R.vec, d_tmp.vec,
                                           CUDA_R_64F,
                                           CUSPARSE_SPSV_ALG_DEFAULT,
                                           spsvDescrL) )
        //    (b) L^-T tmp => R_aux_i+1    (triangular solver)
        CHECK_CUSPARSE( cusparseSpSV_solve(cusparseHandle,
                                           CUSPARSE_OPERATION_TRANSPOSE,
                                           &one, matL, d_tmp.vec,
                                           d_R_aux.vec, CUDA_R_64F,
                                           CUSPARSE_SPSV_ALG_DEFAULT,
                                           spsvDescrLT) )
        //----------------------------------------------------------------------
        // ### 10 ### beta = (R_i+1, R_aux_i+1) / (R_i, R_aux_i)
        //    (a) delta_new => (R_i+1, R_aux_i+1)
        double delta_new;
        CHECK_CUBLAS( cublasDdot(cublasHandle, m, d_R.ptr, 1, d_R_aux.ptr, 1,
                                 &delta_new) )
        //    (b) beta => delta_new / delta
        double beta = delta_new / delta;
        PRINT_INFO(delta_new)
        PRINT_INFO(beta)
        delta       = delta_new;
        //----------------------------------------------------------------------
        // ### 11 ###  P_i+1 = R_aux_i+1 + beta * P_i
        //    (a) P = beta * P
        CHECK_CUBLAS(cublasDscal(cublasHandle, m, &beta, d_P.ptr, 1))
        //    (b) P = R_aux + P
        CHECK_CUBLAS(
            cublasDaxpy(cublasHandle, m, &one, d_R_aux.ptr, 1, d_P.ptr, 1))
    }
    //--------------------------------------------------------------------------
    printf("Check Solution\n"); // ||R = b - A * X||
    //    (a) copy b in R
    CHECK_CUDA( cudaMemcpy(d_R.ptr, d_B.ptr, m * sizeof(double),
                           cudaMemcpyDeviceToDevice) )
    // R = -A * X + R
    CHECK_CUSPARSE( cusparseSpMV(cusparseHandle,
                                 CUSPARSE_OPERATION_NON_TRANSPOSE, &minus_one,
                                 matA, d_X.vec, &one, d_R.vec, CUDA_R_64F,
                                 CUSPARSE_SPMV_ALG_DEFAULT, d_bufferMV) )
    // check ||R||
    CHECK_CUBLAS( cublasDnrm2(cublasHandle, m, d_R.ptr, 1, &nrm_R) )
    printf("Final error norm = %e\n", nrm_R);
    //--------------------------------------------------------------------------
    CHECK_CUSPARSE( cusparseSpSV_destroyDescr(spsvDescrL) )
    CHECK_CUSPARSE( cusparseSpSV_destroyDescr(spsvDescrLT) )
    CHECK_CUDA( cudaFree(d_bufferL) )
    CHECK_CUDA( cudaFree(d_bufferLT) )
    return EXIT_SUCCESS;
}

//==============================================================================
//==============================================================================

#if 0 == 1
int main(int argc, char** argv) {
    const int    maxIterations = 10000;
    const double tolerance     = 1e-8f;
    if (argc != 1) {
        printf("Wrong number of command line arguments. cg_example accepts no arguments.\n");
        return EXIT_FAILURE;
    }
    int     base        = 0;
    int     m           = -1;
    int*    h_A_rows    = NULL;
    int*    h_A_columns = NULL;
    double* h_A_values  = NULL;
    make_laplace_matrix(&m, &h_A_rows, &h_A_columns, &h_A_values);
    int num_offsets = m + 1;
    int nnz = h_A_rows[m];
    double* h_X = (double*)malloc(m * sizeof(double));

    printf("Testing CG\n");
    for (int i = 0; i < m; i++)
        h_X[i] = 1.0;
    //--------------------------------------------------------------------------
    // ### Device memory management ###
    int*    d_A_rows, *d_A_columns;
    double* d_A_values, *d_L_values;
    Vec     d_B, d_X, d_R, d_R_aux, d_P, d_T, d_tmp;

    // allocate device memory for CSR matrices
    CHECK_CUDA( cudaMalloc((void**) &d_A_rows,    num_offsets * sizeof(int)) )
    CHECK_CUDA( cudaMalloc((void**) &d_A_columns, nnz * sizeof(int)) )
    CHECK_CUDA( cudaMalloc((void**) &d_A_values,  nnz * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_L_values,  nnz * sizeof(double)) )

    CHECK_CUDA( cudaMalloc((void**) &d_B.ptr,     m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_X.ptr,     m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_R.ptr,     m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_R_aux.ptr, m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_P.ptr,     m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_T.ptr,     m * sizeof(double)) )
    CHECK_CUDA( cudaMalloc((void**) &d_tmp.ptr,   m * sizeof(double)) )

    // copy the CSR matrices and vectors into device memory
    CHECK_CUDA( cudaMemcpy(d_A_rows, h_A_rows, num_offsets * sizeof(int),
                           cudaMemcpyHostToDevice) )
    CHECK_CUDA( cudaMemcpy(d_A_columns, h_A_columns, nnz *  sizeof(int),
                           cudaMemcpyHostToDevice) )
    CHECK_CUDA( cudaMemcpy(d_A_values, h_A_values, nnz * sizeof(double),
                           cudaMemcpyHostToDevice) )
    CHECK_CUDA( cudaMemcpy(d_L_values, h_A_values, nnz * sizeof(double),
                           cudaMemcpyHostToDevice) )
    CHECK_CUDA( cudaMemcpy(d_X.ptr, h_X, m * sizeof(double),
                           cudaMemcpyHostToDevice) )
    //--------------------------------------------------------------------------
    // ### cuSPARSE Handle and descriptors initialization ###
    // create the test matrix on the host
    cublasHandle_t   cublasHandle   = NULL;
    cusparseHandle_t cusparseHandle = NULL;
    CHECK_CUBLAS( cublasCreate(&cublasHandle) )
    CHECK_CUSPARSE( cusparseCreate(&cusparseHandle) )
    // Create dense vectors
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_B.vec,     m, d_B.ptr, CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_X.vec,     m, d_X.ptr, CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_R.vec,     m, d_R.ptr, CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_R_aux.vec, m, d_R_aux.ptr,
                                        CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_P.vec,   m, d_P.ptr,   CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_T.vec,   m, d_T.ptr,   CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseCreateDnVec(&d_tmp.vec, m, d_tmp.ptr, CUDA_R_64F) )

    cusparseIndexBase_t  baseIdx = CUSPARSE_INDEX_BASE_ZERO;
    cusparseSpMatDescr_t matA, matL;
    int*                 d_L_rows      = d_A_rows;
    int*                 d_L_columns   = d_A_columns;
    cusparseFillMode_t   fill_lower    = CUSPARSE_FILL_MODE_LOWER;
    cusparseDiagType_t   diag_non_unit = CUSPARSE_DIAG_TYPE_NON_UNIT;
    // A
    CHECK_CUSPARSE( cusparseCreateCsr(&matA, m, m, nnz, d_A_rows,
                                      d_A_columns, d_A_values,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      baseIdx, CUDA_R_64F) )
    // L
    CHECK_CUSPARSE( cusparseCreateCsr(&matL, m, m, nnz, d_L_rows,
                                      d_L_columns, d_L_values,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                      baseIdx, CUDA_R_64F) )
    CHECK_CUSPARSE( cusparseSpMatSetAttribute(matL,
                                              CUSPARSE_SPMAT_FILL_MODE,
                                              &fill_lower, sizeof(fill_lower)) )
    CHECK_CUSPARSE( cusparseSpMatSetAttribute(matL,
                                              CUSPARSE_SPMAT_DIAG_TYPE,
                                              &diag_non_unit,
                                              sizeof(diag_non_unit)) )
    //--------------------------------------------------------------------------
    // ### Preparation ### b = A * X
    const double alpha = 0.75;
    size_t       bufferSizeMV;
    void*        d_bufferMV;
    double       beta = 0.0;
    CHECK_CUSPARSE( cusparseSpMV_bufferSize(
                        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha, matA, d_X.vec, &beta, d_B.vec, CUDA_R_64F,
                        CUSPARSE_SPMV_ALG_DEFAULT, &bufferSizeMV) )
    CHECK_CUDA( cudaMalloc(&d_bufferMV, bufferSizeMV) )

    CHECK_CUSPARSE( cusparseSpMV(
                        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha, matA, d_X.vec, &beta, d_B.vec, CUDA_R_64F,
                        CUSPARSE_SPMV_ALG_DEFAULT, d_bufferMV) )
    // X0 = 0
    CHECK_CUDA( cudaMemset(d_X.ptr, 0x0, m * sizeof(double)) )
    //--------------------------------------------------------------------------
    // Perform Incomplete-Cholesky factorization of A (csric0) -> L, L^T
    cusparseMatDescr_t descrM;
    csric02Info_t      infoM        = NULL;
    int                bufferSizeIC = 0;
    void*              d_bufferIC;
    CHECK_CUSPARSE( cusparseCreateMatDescr(&descrM) )
    CHECK_CUSPARSE( cusparseSetMatIndexBase(descrM, baseIdx) )
    CHECK_CUSPARSE( cusparseSetMatType(descrM, CUSPARSE_MATRIX_TYPE_GENERAL) )
    CHECK_CUSPARSE( cusparseSetMatFillMode(descrM, CUSPARSE_FILL_MODE_LOWER) )
    CHECK_CUSPARSE( cusparseSetMatDiagType(descrM,
                                           CUSPARSE_DIAG_TYPE_NON_UNIT) )
    CHECK_CUSPARSE( cusparseCreateCsric02Info(&infoM) )

    CHECK_CUSPARSE( cusparseDcsric02_bufferSize(
                        cusparseHandle, m, nnz, descrM, d_L_values,
                        d_A_rows, d_A_columns, infoM, &bufferSizeIC) )
    CHECK_CUDA( cudaMalloc(&d_bufferIC, bufferSizeIC) )
    CHECK_CUSPARSE( cusparseDcsric02_analysis(
                        cusparseHandle, m, nnz, descrM, d_L_values,
                        d_A_rows, d_A_columns, infoM,
                        CUSPARSE_SOLVE_POLICY_NO_LEVEL, d_bufferIC) )
    int structural_zero;
    CHECK_CUSPARSE( cusparseXcsric02_zeroPivot(cusparseHandle, infoM,
                                               &structural_zero) )
    // M = L * L^T
    CHECK_CUSPARSE( cusparseDcsric02(
                        cusparseHandle, m, nnz, descrM, d_L_values,
                        d_A_rows, d_A_columns, infoM,
                        CUSPARSE_SOLVE_POLICY_NO_LEVEL, d_bufferIC) )
    // Find numerical zero
    int numerical_zero;
    CHECK_CUSPARSE( cusparseXcsric02_zeroPivot(cusparseHandle, infoM,
                                               &numerical_zero) )

    CHECK_CUSPARSE( cusparseDestroyCsric02Info(infoM) )
    CHECK_CUSPARSE( cusparseDestroyMatDescr(descrM) )
    CHECK_CUDA( cudaFree(d_bufferIC) )
    //--------------------------------------------------------------------------
    // ### Run CG computation ###
    printf("CG loop:\n");
    gpu_CG(cublasHandle, cusparseHandle, m,
           matA, matL, d_B, d_X, d_R, d_R_aux, d_P, d_T,
           d_tmp, d_bufferMV, maxIterations, tolerance);
    //--------------------------------------------------------------------------
    // ### Free resources ###
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_B.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_X.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_R.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_R_aux.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_P.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_T.vec) )
    CHECK_CUSPARSE( cusparseDestroyDnVec(d_tmp.vec) )
    CHECK_CUSPARSE( cusparseDestroySpMat(matA) )
    CHECK_CUSPARSE( cusparseDestroySpMat(matL) )
    CHECK_CUSPARSE( cusparseDestroy(cusparseHandle) )
    CHECK_CUBLAS( cublasDestroy(cublasHandle) )

    free(h_A_rows);
    free(h_A_columns);
    free(h_A_values);
    free(h_X);

    CHECK_CUDA( cudaFree(d_X.ptr) )
    CHECK_CUDA( cudaFree(d_B.ptr) )
    CHECK_CUDA( cudaFree(d_R.ptr) )
    CHECK_CUDA( cudaFree(d_R_aux.ptr) )
    CHECK_CUDA( cudaFree(d_P.ptr) )
    CHECK_CUDA( cudaFree(d_T.ptr) )
    CHECK_CUDA( cudaFree(d_tmp.ptr) )
    CHECK_CUDA( cudaFree(d_A_values) )
    CHECK_CUDA( cudaFree(d_A_columns) )
    CHECK_CUDA( cudaFree(d_A_rows) )
    CHECK_CUDA( cudaFree(d_L_values) )
    CHECK_CUDA( cudaFree(d_bufferMV) )
    return EXIT_SUCCESS;
}
#endif

//==============================================================================
