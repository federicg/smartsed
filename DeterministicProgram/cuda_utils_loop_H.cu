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

double deviceMin(const double* raw_ptr, size_t n) {
    thrust::device_ptr<const double> dev_ptr(raw_ptr);
    return *thrust::min_element(dev_ptr, dev_ptr + n);
}

double deviceMax(const double* raw_ptr, size_t n) {
    thrust::device_ptr<const double> dev_ptr(raw_ptr);
    return *thrust::max_element(dev_ptr, dev_ptr + n);
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
// Horizontal upwind, merged internal+West+East: single launch, tag-dispatched
// (0=internal, 1=West, 2=East). See build_merged_ids_and_tags in utils_H.h.
__global__ void computeHorizontalMergedKernel_interface(
    const unsigned int* ids,
    const unsigned int* tag,
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
    unsigned int t  = tag[i];

    double H_left, H_right, s;
    if (t == 0) {
        unsigned int IDeast = Id - ii;
        unsigned int IDwest = Id - ii - 1;
        H_left  = H[IDwest];
        H_right = H[IDeast];
        s = (u[Id] > 0.0) - (u[Id] < 0.0);
    } else if (t == 1) {
        H_left  = 0.0;
        H_right = H[Id - ii];
        s = (u[Id + 1] > 0.0) - (u[Id + 1] < 0.0);
    } else {
        H_left  = H[Id - ii - 1];
        H_right = 0.0;
        s = (u[Id - 1] > 0.0) - (u[Id - 1] < 0.0);
    }

    horizontal[Id] = 0.5 * (H_left + H_right + s * (H_left - H_right));
}

void compute_horizontal_interface_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		const thrust::device_vector<unsigned int>& tag,
		const thrust::device_vector<double>& H,
		const thrust::device_vector<double>& u,
		thrust::device_vector<double>& horizontal, const unsigned int N_cols,
		cudaStream_t stream) {

    launch_kernel(
      computeHorizontalMergedKernel_interface,
      idAll.size(),
      stream,
      thrust::raw_pointer_cast(idAll.data()),
      thrust::raw_pointer_cast(tag.data()),
      thrust::raw_pointer_cast(H.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(horizontal.data()),
      N_cols
    );
}

// Vertical upwind, merged internal+North+South (0=internal, 1=North, 2=South).
__global__ void computeVerticalMergedKernel_interface(
    const unsigned int* ids,
    const unsigned int* tag,
    const double* H,
    const double* v,
    double* vertical,
    unsigned int N_cols,
    unsigned int n)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  unsigned int Id = ids[i];
  unsigned int t  = tag[i];

  double H_left, H_right, s;
  if (t == 0) {
    unsigned int IDsouth = Id, IDnorth = Id - N_cols;
    H_left = H[IDnorth];
    H_right = H[IDsouth];
    double v_cc = v[Id];
    s = (v_cc > 0.0) - (v_cc < 0.0);
  } else if (t == 1) {
    H_left = 0;
    H_right = H[Id];
    double v_right = v[Id + N_cols];
    s = (v_right > 0.0) - (v_right < 0.0);
  } else {
    H_left = H[Id - N_cols];
    H_right = 0;
    double v_left = v[Id - N_cols];
    s = (v_left > 0.0) - (v_left < 0.0);
  }

  vertical[Id] = (H_left + H_right) * .5 + s * (H_left - H_right) * .5;
}

void compute_vertical_interface_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		const thrust::device_vector<unsigned int>& tag,
		const thrust::device_vector<double>& H,
		const thrust::device_vector<double>& v,
		thrust::device_vector<double>& vertical, const unsigned int N_cols,
		cudaStream_t stream) {

    launch_kernel(
      computeVerticalMergedKernel_interface,
      idAll.size(),
      stream,
      thrust::raw_pointer_cast(idAll.data()),
      thrust::raw_pointer_cast(tag.data()),
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

// computeKernel_friction's formula doesn't distinguish internal from
// boundary faces (it depends only on already-fully-populated per-face
// arrays), so a single launch over the concatenated id list is equivalent
// to the three separate internal/West-or-North/East-or-South launches.
void compute_friction_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		const thrust::device_vector<double>& H_interface_dir,
		const thrust::device_vector<double>& vel,
		const thrust::device_vector<double>& M_expo_r_dir_vect,
		      thrust::device_vector<double>& alfa_dir,
		const thrust::device_vector<double>& M_gamma_dt_DSV_dir_,
		double M_dt_DSV, double M_coeff, double M_H_min,
                double M_expo, unsigned int M_frictionModel, cudaStream_t stream) {

    launch_kernel(
      computeKernel_friction,
      idAll.size(),
      stream,
      thrust::raw_pointer_cast(idAll.data()),
      thrust::raw_pointer_cast(H_interface_dir.data()),
      thrust::raw_pointer_cast(vel.data()),
      thrust::raw_pointer_cast(M_expo_r_dir_vect.data()),
      thrust::raw_pointer_cast(alfa_dir.data()),
      thrust::raw_pointer_cast(M_gamma_dt_DSV_dir_.data()),
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

// Horizontal gravitational layer, merged internal+West+East: single launch,
// tag-dispatched (0=internal, 1=West, 2=East). Each branch only sets
// h_left/h_right/k_c_left/k_c_right; the closing expression is the verbatim
// formula shared by the three kernels this replaces.
__global__ void computeHorizontalMergedKernel_gravitational(
    const unsigned int* ids,
    const unsigned int* tag,
    const double* coeff,
    const double* n_x,
    const double* h,
    double* h_interface_x,
    const unsigned int N_cols,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
  const unsigned int t  = tag[i];
  const unsigned int ii = Id / (N_cols + 1);

  double h_left, h_right, k_c_left, k_c_right;
  if (t == 0) {
    const unsigned int IDeast = Id - ii, IDwest = Id - ii - 1;
    h_left  = h[IDwest];      h_right  = h[IDeast];
    k_c_left = coeff[IDwest]; k_c_right = coeff[IDeast];
  } else if (t == 1) {
    h_left  = 0.;            h_right  = h[Id - ii];
    k_c_left = 0.;           k_c_right = coeff[Id - ii];
  } else {
    h_left  = h[Id - ii - 1]; h_right  = 0.;
    k_c_left = coeff[Id - ii - 1]; k_c_right = 0.;
  }

  h_interface_x[Id] = n_x[Id] *
	  ((k_c_left * h_left + k_c_right * h_right) +
	   n_x[Id] * (k_c_left * h_left - k_c_right * h_right)) *
	  .5;
}


// Horizontal gravitational layer
void computeResidualsHorizontal_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		const thrust::device_vector<unsigned int>& tag,
		const thrust::device_vector<double>& coeff,
		const thrust::device_vector<double>& n_x,
		const thrust::device_vector<double>& h,
		thrust::device_vector<double>& h_interface_x,
                const unsigned int N_cols, cudaStream_t stream) {

    launch_kernel(
      computeHorizontalMergedKernel_gravitational,
      idAll.size(),
      stream,
      thrust::raw_pointer_cast(idAll.data()),
      thrust::raw_pointer_cast(tag.data()),
      thrust::raw_pointer_cast(coeff.data()),
      thrust::raw_pointer_cast(n_x.data()),
      thrust::raw_pointer_cast(h.data()),
      thrust::raw_pointer_cast(h_interface_x.data()),
      N_cols
    );
}

// Vertical gravitational layer, merged internal+North+South (0/1/2 tag).
__global__ void computeVerticalMergedKernel_gravitational(
    const unsigned int* ids,
    const unsigned int* tag,
    const double* coeff,
    const double* n_y,
    const double* h,
    double* h_interface_y,
    const unsigned int N_cols,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
  const unsigned int t = tag[i];

  double h_left, h_right, k_c_left, k_c_right;
  if (t == 0) {
    const unsigned int IDsouth = Id, IDnorth = Id - N_cols;
    h_left  = h[IDnorth];      h_right  = h[IDsouth];
    k_c_left = coeff[IDnorth]; k_c_right = coeff[IDsouth];
  } else if (t == 1) {
    h_left  = 0.;             h_right  = h[Id];
    k_c_left = 0.;            k_c_right = coeff[Id];
  } else {
    h_left  = h[Id - N_cols]; h_right  = 0.;
    k_c_left = coeff[Id - N_cols]; k_c_right = 0.;
  }

  h_interface_y[Id] = n_y[Id] *
	  ((k_c_left * h_left + k_c_right * h_right) +
	   n_y[Id] * (k_c_left * h_left - k_c_right * h_right)) *
	  .5;
}

// Vertical gravitational layer
void computeResidualsVertical_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		const thrust::device_vector<unsigned int>& tag,
		const thrust::device_vector<double>& coeff,
		const thrust::device_vector<double>& n_y,
		const thrust::device_vector<double>& h,
		thrust::device_vector<double>& h_interface_y,
                const unsigned int N_cols, cudaStream_t stream) {

    launch_kernel(
      computeVerticalMergedKernel_gravitational,
      idAll.size(),
      stream,
      thrust::raw_pointer_cast(idAll.data()),
      thrust::raw_pointer_cast(tag.data()),
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

// Internal-face back-trace. `ids`/`tag` are the concatenated internal+West+East
// list: this kernel is launched over the whole list but only processes tag==0
// (internal) entries, so the West/East kernel -- which reads u_star values this
// kernel writes -- can safely run as a second launch on the same stream over
// the same list. (Merging all three into one launch would race.)
__global__ void bilinearInterpolationHorizontal(
    const unsigned int* __restrict__ ids,
    const unsigned int* __restrict__ tag,
    const double* __restrict__ u,
    const double* __restrict__ v,
    double* __restrict__ u_star,
    const double scale,
    const unsigned int nrows,
    const unsigned int ncols,
    unsigned int n) {

  const unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= n) return;
  if (__ldg(&tag[tid]) != 0u) return;

  const unsigned int Id  = __ldg(&ids[tid]);
  const unsigned int row = Id / (ncols + 1);
  const unsigned int col = Id - row * (ncols + 1);  // cheaper than Id % (ncols+1)

  // Neighbour indices into v-grid (nrows+1) x (ncols)
  const unsigned int ID_NE = Id - row;
  const unsigned int ID_NW = ID_NE - 1;
  const unsigned int ID_SE = ID_NE + ncols;
  const unsigned int ID_SW = ID_NW + ncols;

  // Velocity at this staggered node (u-grid)
  const double vel_x = __ldg(&u[Id]);
  const double vel_y = (__ldg(&v[ID_NE]) + __ldg(&v[ID_NW]) +
                        __ldg(&v[ID_SE]) + __ldg(&v[ID_SW])) * 0.25;

  // Back-trace position
  double x = (double)col - vel_x * scale;
  double y = (double)row - vel_y * scale;

  double x_1 = floor(x), y_1 = floor(y);
  double x_2 = x_1 + 1., y_2 = y_1 + 1.;

  // Same edge nudging as the CPU reference: when the departure point lands
  // exactly on (or within one cell of) the last valid grid line, shift the
  // stencil in by one cell instead of reading out of bounds.
  if (x_1 == (double)ncols)          { x_2 -= 1.; x_1 -= 1.; x -= 1.; }
  if (x_2 == 0.)                     { x_2 += 1.; x_1 += 1.; x += 1.; }
  if (y_1 == (double)(nrows - 1))    { y_2 -= 1.; y_1 -= 1.; y -= 1.; }
  if (y_2 == 0.)                     { y_2 += 1.; y_1 += 1.; y += 1.; }

  // Departure point more than one cell beyond the domain edge: no upstream
  // information is available, so treat it as zero (matches CPU reference).
  if (x_2 < 0. || x_1 > (double)ncols || y_2 < 0. || y_1 > (double)(nrows - 1)) {
    u_star[Id] = 0.;
    return;
  }

  // Flat indices into u-grid, stride is (ncols+1)
  const int stride = (int)(ncols + 1);
  const int Id_11  = (int)x_1 + (int)y_1 * stride;
  const int Id_12  = (int)x_1 + (int)y_2 * stride;
  const int Id_21  = (int)x_2 + (int)y_1 * stride;
  const int Id_22  = (int)x_2 + (int)y_2 * stride;

  // Bilinear weights
  const double w_x2 = x_2 - x;
  const double w_x1 = x - x_1;
  const double w_y2 = y_2 - y;
  const double w_y1 = y - y_1;

  // Interpolate with FMA
  const double a = fma(__ldg(&u[Id_11]), w_x2, __ldg(&u[Id_21]) * w_x1);
  const double b = fma(__ldg(&u[Id_12]), w_x2, __ldg(&u[Id_22]) * w_x1);

  u_star[Id] = fma(a, w_y2, b * w_y1);
}


// Merged West+East boundary back-trace (tag 1=West, 2=East). Launched over the
// full concatenated list as a second pass after bilinearInterpolationHorizontal;
// each branch is the verbatim body of the kernel it replaces.
__global__ void bilinearInterpolationHorizontalBoundary(
		const unsigned int* ids,
                const unsigned int* tag,
                const double* u,
		const double* v,
                double* u_star,
		const double scale,
		const unsigned int N_rows, const unsigned int N_cols,
		unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  const unsigned int t = tag[i];
  if (t == 0u) return;

  const auto Id = ids[i];
  if (t == 1u) {                       // West
    const auto Idd = Id + 1;
    u_star[Id] = u_star[Idd] * (u[Idd] < 0.);
  } else {                             // East
    const auto Idd = Id - 1;
    u_star[Id] = u_star[Idd] * (u[Idd] > 0.);
  }
}


void bilinearInterpolationHorizontal_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
                const thrust::device_vector<unsigned int>& tag,
                const thrust::device_vector<double>& u,
		const thrust::device_vector<double>& v,
		thrust::device_vector<double>& u_star,
		const double scale,
		const unsigned int N_rows, const unsigned int N_cols, cudaStream_t stream) {

    launch_kernel(
      bilinearInterpolationHorizontal,
      idAll.size(),
      stream,
      thrust::raw_pointer_cast(idAll.data()),
      thrust::raw_pointer_cast(tag.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(u_star.data()),
      scale, N_rows, N_cols
    );

    launch_kernel(
      bilinearInterpolationHorizontalBoundary,
      idAll.size(),
      stream,
      thrust::raw_pointer_cast(idAll.data()),
      thrust::raw_pointer_cast(tag.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(u_star.data()),
      scale, N_rows, N_cols
    );
}


// Internal-face back-trace; see bilinearInterpolationHorizontal for why this is
// launched over the whole internal+North+South list but only handles tag==0.
__global__ void bilinearInterpolationVertical(
    const unsigned int* __restrict__ ids,
    const unsigned int* __restrict__ tag,
    const double* __restrict__ u,
    const double* __restrict__ v,
    double* __restrict__ v_star,
    const double scale,
    const unsigned int nrows,
    const unsigned int ncols,
    unsigned int n) {

  const unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= n) return;
  if (__ldg(&tag[tid]) != 0u) return;

  const unsigned int Id   = __ldg(&ids[tid]);

  const unsigned int row = Id / ncols;
  const unsigned int col = Id - row * ncols; //Id % ncols;

  const unsigned int ID_SW = Id + row;
  const unsigned int ID_SE = ID_SW + 1;
  const unsigned int ID_NW = ID_SW - (ncols + 1);
  const unsigned int ID_NE = ID_NW + 1;

  // Back-trace position
  const double vel_x = (__ldg(&u[ID_SW]) + __ldg(&u[ID_SE]) +
                        __ldg(&u[ID_NW]) + __ldg(&u[ID_NE])) * 0.25;
  const double vel_y = __ldg(&v[Id]);

  double x = (double)col - vel_x * scale;
  double y = (double)row - vel_y * scale;

  double x_1 = floor(x), y_1 = floor(y);
  double x_2 = x_1 + 1., y_2 = y_1 + 1.;

  // Same edge nudging as the CPU reference: when the departure point lands
  // exactly on (or within one cell of) the last valid grid line, shift the
  // stencil in by one cell instead of reading out of bounds.
  if (x_1 == (double)(ncols - 1))    { x_2 -= 1.; x_1 -= 1.; x -= 1.; }
  if (x_2 == 0.)                     { x_2 += 1.; x_1 += 1.; x += 1.; }
  if (y_1 == (double)nrows)          { y_2 -= 1.; y_1 -= 1.; y -= 1.; }
  if (y_2 == 0.)                     { y_2 += 1.; y_1 += 1.; y += 1.; }

  // Departure point more than one cell beyond the domain edge: no upstream
  // information is available, so treat it as zero (matches CPU reference).
  if (x_2 < 0. || x_1 > (double)(ncols - 1) || y_2 < 0. || y_1 > (double)nrows) {
    v_star[Id] = 0.;
    return;
  }

  // Flat indices into v-grid
  const int Id_11 = (int)x_1 + (int)y_1 * (int)ncols;
  const int Id_12 = (int)x_1 + (int)y_2 * (int)ncols;
  const int Id_21 = (int)x_2 + (int)y_1 * (int)ncols;
  const int Id_22 = (int)x_2 + (int)y_2 * (int)ncols;

  // Bilinear weights
  const double w_x2 = x_2 - x;
  const double w_x1 = x - x_1;
  const double w_y2 = y_2 - y;
  const double w_y1 = y - y_1;

  // Interpolate with FMA
  const double a = fma(__ldg(&v[Id_11]), w_x2, __ldg(&v[Id_21]) * w_x1);
  const double b = fma(__ldg(&v[Id_12]), w_x2, __ldg(&v[Id_22]) * w_x1);

  v_star[Id] = fma(a, w_y2, b * w_y1);
}


// Merged North+South boundary back-trace (tag 1=North, 2=South); second pass
// after bilinearInterpolationVertical, each branch verbatim from the kernel
// it replaces.
__global__ void bilinearInterpolationVerticalBoundary(
		const unsigned int* ids,
                const unsigned int* tag,
                const double* u,
		const double* v,
                double* v_star,
		const double scale,
		const unsigned int nrows, const unsigned int ncols,
		unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  const unsigned int t = tag[i];
  if (t == 0u) return;

  const auto Id = ids[i];
  if (t == 1u) {                       // North
    const auto Idd = Id + ncols;
    v_star[Id] = v_star[Idd] * (v[Idd] < 0.);
  } else {                             // South
    const auto Idd = Id - ncols;
    v_star[Id] = v_star[Idd] * (v[Idd] > 0.);
  }
}



void bilinearInterpolationVertical_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
                const thrust::device_vector<unsigned int>& tag,
                const thrust::device_vector<double>& u,
		const thrust::device_vector<double>& v,
		thrust::device_vector<double>& v_star,
		const double scale,
		const unsigned int N_rows, const unsigned int N_cols, cudaStream_t stream) {

    launch_kernel(
      bilinearInterpolationVertical,
      idAll.size(),
      stream,
      thrust::raw_pointer_cast(idAll.data()),
      thrust::raw_pointer_cast(tag.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(v_star.data()),
      scale, N_rows, N_cols
    );

    launch_kernel(
      bilinearInterpolationVerticalBoundary,
      idAll.size(),
      stream,
      thrust::raw_pointer_cast(idAll.data()),
      thrust::raw_pointer_cast(tag.data()),
      thrust::raw_pointer_cast(u.data()),
      thrust::raw_pointer_cast(v.data()),
      thrust::raw_pointer_cast(v_star.data()),
      scale, N_rows, N_cols
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

// Sediment residuals, one direction. computeKernel_sediment's formula depends
// only on already-fully-populated per-face arrays (S_dir_mod, vel) and so is
// identical for internal and boundary faces -- a single launch over the
// concatenated id list is equivalent to the three internal/West-or-North/
// East-or-South launches it replaces (same rationale as compute_friction_wrapper).
void computeResidualsTruncated_wrapper(
		const thrust::device_vector<unsigned int>& idAll,
		thrust::device_vector<double>& Gamma_dir_1, thrust::device_vector<double>& Gamma_dir_2,
		const thrust::device_vector<double>& S_dir_mod,
		const thrust::device_vector<double>& vel,
                const double c1, cudaStream_t stream) {

    launch_kernel(
      computeKernel_sediment,
      idAll.size(),
      stream,
      thrust::raw_pointer_cast(idAll.data()),
      thrust::raw_pointer_cast(Gamma_dir_1.data()),
      thrust::raw_pointer_cast(Gamma_dir_2.data()),
      thrust::raw_pointer_cast(S_dir_mod.data()),
      thrust::raw_pointer_cast(vel.data()),
      c1
    );
}

//==============================================================================
// Gravitational-layer sediment production W_Gav (EPM). Elementwise over basin
// cells. The CPU reference adds a per-pour-point additional_source_term coming
// from the static excluded-subbasin approximation; that approximation is not
// carried on the GPU, so this kernel computes only the base EPM term.
__global__ void computeKernel_WGav(
    const unsigned int* ids,
    const double* Z_Gav,
    const double* T_raster,
    const double* melt_mask,
    const double* DP_total,
    double* W_Gav,
    const double pi_1e3,   // 1.e-3 * M_PI, folded on the host to match the CPU
    const double dt_sed,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto k = ids[i];

  W_Gav[k] = pi_1e3 * Z_Gav[k] *
             sqrt(fabs((.1 + .1 * T_raster[k]) * melt_mask[k])) *
             DP_total[k] * dt_sed;
}

void computeWGav_wrapper(
    const thrust::device_vector<unsigned int>& idBasinVect,
    const thrust::device_vector<double>& Z_Gav,
    const thrust::device_vector<double>& T_raster,
    const thrust::device_vector<double>& melt_mask,
    const thrust::device_vector<double>& DP_total,
          thrust::device_vector<double>& W_Gav,
    const double pi_1e3, const double dt_sed, cudaStream_t stream) {

    launch_kernel(
      computeKernel_WGav,
      idBasinVect.size(),
      stream,
      thrust::raw_pointer_cast(idBasinVect.data()),
      thrust::raw_pointer_cast(Z_Gav.data()),
      thrust::raw_pointer_cast(T_raster.data()),
      thrust::raw_pointer_cast(melt_mask.data()),
      thrust::raw_pointer_cast(DP_total.data()),
      thrust::raw_pointer_cast(W_Gav.data()),
      pi_1e3, dt_sed
    );
}

//==============================================================================
// Sediment transport sub-stepping (EPM). Per sub-step: assemble the h_sd
// advection fluxes h_interface_{x,y} from the Gamma coefficients and the
// current h_sd, then update h_sd. The CPU reference also routes fluxes from
// excluded (static-subbasin) cells to their pour point via additional_source
// _term; that approximation is not carried on the GPU, so it is dropped here
// (== 0 everywhere on the GPU).
//
// Face-flux kernels: same tag layout as the gravitational merge
// (0=internal, 1=West/North, 2=East/South) and the same index stencils; only
// the closing formula differs (Gamma_1*h_right + Gamma_2*h_left, with
// Gamma_1/Gamma_2 the right/left coefficients from computeResidualsTruncated).

__global__ void computeKernel_sedFluxHorizontal(
    const unsigned int* ids,
    const unsigned int* tag,
    const double* Gamma_x_1,   // right coeff
    const double* Gamma_x_2,   // left coeff
    const double* h_sd,
    double* h_interface_x,
    const unsigned int N_cols,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
  const unsigned int t  = tag[i];
  const unsigned int ii = Id / (N_cols + 1);

  double h_left, h_right;
  if (t == 0) {              // internal
    h_left  = h_sd[Id - ii - 1];
    h_right = h_sd[Id - ii];
  } else if (t == 1) {       // West
    h_left  = 0.0;
    h_right = h_sd[Id - ii];
  } else {                   // East
    h_left  = h_sd[Id - ii - 1];
    h_right = 0.0;
  }

  h_interface_x[Id] = Gamma_x_1[Id] * h_right + Gamma_x_2[Id] * h_left;
}

__global__ void computeKernel_sedFluxVertical(
    const unsigned int* ids,
    const unsigned int* tag,
    const double* Gamma_y_1,
    const double* Gamma_y_2,
    const double* h_sd,
    double* h_interface_y,
    const unsigned int N_cols,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
  const unsigned int t = tag[i];

  double h_left, h_right;
  if (t == 0) {              // internal: IDnorth = Id - N_cols, IDsouth = Id
    h_left  = h_sd[Id - N_cols];
    h_right = h_sd[Id];
  } else if (t == 1) {       // North
    h_left  = 0.0;
    h_right = h_sd[Id];
  } else {                   // South
    h_left  = h_sd[Id - N_cols];
    h_right = 0.0;
  }

  h_interface_y[Id] = Gamma_y_1[Id] * h_right + Gamma_y_2[Id] * h_left;
}

__global__ void computeKernel_sedUpdate(
    const unsigned int* ids,
    const double* h_interface_x,
    const double* h_interface_y,
    const double* W_Gav,
    double* h_sd,
    const unsigned int N_cols,
    unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
  const unsigned int ii = Id / N_cols;

  const double Res_x_cell = h_interface_x[Id + 1 + ii] - h_interface_x[Id + ii];
  const double Res_y_cell = h_interface_y[Id + N_cols] - h_interface_y[Id];

  h_sd[Id] += -(Res_x_cell + Res_y_cell) + W_Gav[Id];
}

void computeSedimentTransport_wrapper(
    const unsigned int numberOfSteps,
    const thrust::device_vector<unsigned int>& idHorizontalAll,
    const thrust::device_vector<unsigned int>& horizontalTag,
    const thrust::device_vector<unsigned int>& idVerticalAll,
    const thrust::device_vector<unsigned int>& verticalTag,
    const thrust::device_vector<unsigned int>& idBasinVect,
    const thrust::device_vector<double>& Gamma_x_1,
    const thrust::device_vector<double>& Gamma_x_2,
    const thrust::device_vector<double>& Gamma_y_1,
    const thrust::device_vector<double>& Gamma_y_2,
          thrust::device_vector<double>& h_sd,
          thrust::device_vector<double>& h_interface_x,
          thrust::device_vector<double>& h_interface_y,
    const thrust::device_vector<double>& W_Gav,
    const unsigned int N_cols, cudaStream_t stream) {

  const double* gx1 = thrust::raw_pointer_cast(Gamma_x_1.data());
  const double* gx2 = thrust::raw_pointer_cast(Gamma_x_2.data());
  const double* gy1 = thrust::raw_pointer_cast(Gamma_y_1.data());
  const double* gy2 = thrust::raw_pointer_cast(Gamma_y_2.data());
  double* hsd  = thrust::raw_pointer_cast(h_sd.data());
  double* hix  = thrust::raw_pointer_cast(h_interface_x.data());
  double* hiy  = thrust::raw_pointer_cast(h_interface_y.data());
  const double* wg = thrust::raw_pointer_cast(W_Gav.data());

  // All three launches per sub-step share `stream`, so on that stream they
  // run in order: flux_x, flux_y, then the h_sd update that reads both, then
  // the next sub-step's flux kernels that read the updated h_sd.
  for (unsigned int kk = 0; kk < numberOfSteps; ++kk) {
    launch_kernel(computeKernel_sedFluxHorizontal, idHorizontalAll.size(), stream,
      thrust::raw_pointer_cast(idHorizontalAll.data()),
      thrust::raw_pointer_cast(horizontalTag.data()),
      gx1, gx2, hsd, hix, N_cols);

    launch_kernel(computeKernel_sedFluxVertical, idVerticalAll.size(), stream,
      thrust::raw_pointer_cast(idVerticalAll.data()),
      thrust::raw_pointer_cast(verticalTag.data()),
      gy1, gy2, hsd, hiy, N_cols);

    launch_kernel(computeKernel_sedUpdate, idBasinVect.size(), stream,
      thrust::raw_pointer_cast(idBasinVect.data()),
      hix, hiy, wg, hsd, N_cols);
  }
}

//==============================================================================

__device__ int findPosition(const int* d_A_rows, const int* d_A_columns,
                             int row, int col) {

    // linear search method: \mathcal{O}(k), replace with binary search method for wider stencil
    for (int k = d_A_rows[row]; k < d_A_rows[row + 1]; k++)
        if (d_A_columns[k] == col)
            return k;
    return -1;
}

// A miss means the precomputed CSR sparsity pattern does not contain a
// (row,col) slot that a buildMatrix_* kernel needs. Writing to d_A_values[-1]
// would silently corrupt device memory instead of failing, so trap instead:
// this aborts the kernel (and the CUDA context) loudly rather than letting
// an out-of-bounds write surface later as an unrelated cuSPARSE failure.
__device__ __forceinline__ void checkPosition(int pos, int row, int col) {
    if (pos < 0) {
        printf("buildMatrix: CSR position not found for (row=%d, col=%d) "
               "-- sparsity pattern mismatch, aborting.\n", row, col);
        __trap();
    }
}

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
	      unsigned int n) {
 
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const auto Id = ids[i];
 
  const auto IDreIndex = idBasinVectReIndex[Id];

  // Bring the first appetizer to the rhs
  rhs[IDreIndex] = H[Id] + precipitation[Id] * dt_DSV;

  // Initialize now the diagonal of the system matrix
  int pos_diag = findPosition(d_A_rows, d_A_columns, IDreIndex, IDreIndex);
  checkPosition(pos_diag, IDreIndex, IDreIndex);
  d_A_values[pos_diag] = 1.0;
}


// Cell-centered gather replacing the old face-centered scatter
// (buildMatrix_horizontal_internal + buildMatrix_vertical_internal). One
// thread per basin cell owns exactly one CSR row and one rhs slot, so every
// write here is a plain (non-atomic) read-modify-write: no other thread ever
// touches this row's diagonal/off-diagonal entries or this cell's rhs slot.
// Boundary faces (this cell's neighbor is outside the basin) are left to the
// unchanged buildMatrix_{horizontal,vertical}_{West,East,North,South}
// kernels, exactly as before -- this kernel only ever handles internal
// (basin-basin) neighbor pairs.
__global__ void buildMatrix_cell_gather(
		const unsigned int* ids,               // idBasinVect, size n=m
                const double* H_int_x,
                const double* H_int_y,
		const double* alfa_x,
		const double* alfa_y,
		const double* orography,
		const double* u_star,
		const double* v_star,
		const unsigned int* basin_mask,        // size N (raw grid)
                const unsigned int* idBasinVectReIndex,
		const unsigned int N_cols,
		const unsigned int N,                  // basin_mask.size()
		const double c1,
		const double c3,
		const double H_min,
		double* rhs,
		const int* d_A_rows,
		const int* d_A_columns,
		double* d_A_values,
		unsigned int n) {

  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= n) return;

  const unsigned int k  = ids[idx];
  const unsigned int ii = k / N_cols;
  const unsigned int jj = k - ii * N_cols;
  const unsigned int R  = idBasinVectReIndex[k];

  const int pos_diag = findPosition(d_A_rows, d_A_columns, (int)R, (int)R);
  checkPosition(pos_diag, R, R);

  double diag_acc = 0.0;
  double rhs_acc  = 0.0;

  // ---- West ----
  if (jj > 0 && basin_mask[k - 1] == 1) {
    const unsigned int face   = k + ii;      // my west face (u-grid)
    const unsigned int neigh  = k - 1;
    const unsigned int neighR = idBasinVectReIndex[neigh];
    const double H_interface  = H_int_x[face];
    const double coeff_m      = H_interface * alfa_x[face];
    const bool   is_nonDry    = H_interface > H_min;

    const int pos_off = findPosition(d_A_rows, d_A_columns, (int)R, (int)neighR);
    checkPosition(pos_off, R, neighR);

    if (is_nonDry) {
      diag_acc += c3 * coeff_m;
      rhs_acc  += -c1 * (-coeff_m * u_star[face]) -
                  (orography[k] - orography[neigh]) * c3 * coeff_m;
      d_A_values[pos_off] = -c3 * coeff_m;
    } else {
      diag_acc += 1e-6;
      d_A_values[pos_off] = 0.0;
    }
  }

  // ---- East ----
  if (jj < N_cols - 1 && basin_mask[k + 1] == 1) {
    const unsigned int face   = k + ii + 1;  // my east face (u-grid)
    const unsigned int neigh  = k + 1;
    const unsigned int neighR = idBasinVectReIndex[neigh];
    const double H_interface  = H_int_x[face];
    const double coeff_m      = H_interface * alfa_x[face];
    const bool   is_nonDry    = H_interface > H_min;

    const int pos_off = findPosition(d_A_rows, d_A_columns, (int)R, (int)neighR);
    checkPosition(pos_off, R, neighR);

    if (is_nonDry) {
      diag_acc += c3 * coeff_m;
      rhs_acc  += -c1 * (+coeff_m * u_star[face]) -
                  (orography[k] - orography[neigh]) * c3 * coeff_m;
      d_A_values[pos_off] = -c3 * coeff_m;
    } else {
      diag_acc += 1e-6;
      d_A_values[pos_off] = 0.0;
    }
  }

  // ---- North ----
  if (k >= N_cols && basin_mask[k - N_cols] == 1) {
    const unsigned int face   = k;           // my north face (v-grid)
    const unsigned int neigh  = k - N_cols;
    const unsigned int neighR = idBasinVectReIndex[neigh];
    const double H_interface  = H_int_y[face];
    const double coeff_m      = H_interface * alfa_y[face];
    const bool   is_nonDry    = H_interface > H_min;

    const int pos_off = findPosition(d_A_rows, d_A_columns, (int)R, (int)neighR);
    checkPosition(pos_off, R, neighR);

    if (is_nonDry) {
      diag_acc += c3 * coeff_m;
      rhs_acc  += -c1 * (-coeff_m * v_star[face]) -
                  (orography[k] - orography[neigh]) * c3 * coeff_m;
      d_A_values[pos_off] = -c3 * coeff_m;
    } else {
      diag_acc += 1e-6;
      d_A_values[pos_off] = 0.0;
    }
  }

  // ---- South ----
  if (k + N_cols < N && basin_mask[k + N_cols] == 1) {
    const unsigned int face   = k + N_cols;  // my south face (v-grid)
    const unsigned int neigh  = k + N_cols;
    const unsigned int neighR = idBasinVectReIndex[neigh];
    const double H_interface  = H_int_y[face];
    const double coeff_m      = H_interface * alfa_y[face];
    const bool   is_nonDry    = H_interface > H_min;

    const int pos_off = findPosition(d_A_rows, d_A_columns, (int)R, (int)neighR);
    checkPosition(pos_off, R, neighR);

    if (is_nonDry) {
      diag_acc += c3 * coeff_m;
      rhs_acc  += -c1 * (+coeff_m * v_star[face]) -
                  (orography[k] - orography[neigh]) * c3 * coeff_m;
      d_A_values[pos_off] = -c3 * coeff_m;
    } else {
      diag_acc += 1e-6;
      d_A_values[pos_off] = 0.0;
    }
  }

  // Only this thread ever touches row R / rhs[R] beyond buildMatrix_cell_center
  // (which ran earlier on the same stream), so a plain += is safe here.
  d_A_values[pos_diag] += diag_acc;
  rhs[R] += rhs_acc;
}

// Merged boundary rhs contribution, one direction. Launched twice: once over
// the horizontal (internal+West+East) list with H_int_x/alfa_x/u_star and
// isHorizontal=1, once over the vertical (internal+North+South) list with the
// _y arrays and isHorizontal=0. tag==0 (internal) entries are skipped -- they
// are handled by buildMatrix_cell_gather. This replaces the four separate
// West/East/North/South kernels; each branch below is the verbatim formula of
// the kernel it replaces (West==North sign, East==South sign; the reindex
// stencil differs per direction/side).
//
// rhs is updated with atomicAdd because, unlike the old one-kernel-per-side
// layout, West and East (resp. North and South) now run in the same launch and
// a 1-cell-wide basin neck would let both touch the same row. Contention is at
// most 2 and only on perimeter cells, so this is not the O(N) scatter the
// internal assembly was deliberately de-atomized away from.
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
		unsigned int n) {

  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;

  const unsigned int t = tag[i];
  if (t == 0u) return;                       // internal: see buildMatrix_cell_gather

  const auto Id = ids[i];

  unsigned int IDReIndex;
  if (isHorizontal) {
    const unsigned int ii = Id / (N_cols + 1);
    IDReIndex = idBasinVectReIndex[(t == 1u) ? (Id - ii) : (Id - ii - 1)]; // West : East
  } else {
    IDReIndex = idBasinVectReIndex[(t == 1u) ? Id : (Id - N_cols)];        // North : South
  }

  const auto H_interface = H_int[Id];
  const double coeff_m = H_interface * alfa[Id];

  if (H_interface > H_min) {
    // West/North: -coeff_m ; East/South: +coeff_m  (negation is exact in IEEE754)
    const double signed_coeff_m = (t == 1u) ? -coeff_m : coeff_m;
    atomicAdd(&rhs[IDReIndex],
              isNonReflectingBC * (-c1 * (signed_coeff_m * vel_star[Id])));
  }
}




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
		Vec& rhs, cudaStream_t stream) {

    // zero BEFORE kernels
    cudaMemsetAsync(d_A_values, 0, nnz * sizeof(double), stream);

    //--------------------------------------------------------------------------
    // Cell-centered contributions
    launch_kernel(
	buildMatrix_cell_center,
	idBasinVect.size(),
	stream,
	thrust::raw_pointer_cast(idBasinVect.data()),
        thrust::raw_pointer_cast(H.data()),
	thrust::raw_pointer_cast(precipitation.data()),
	thrust::raw_pointer_cast(idBasinVectReIndex.data()),
	dt_DSV, rhs.ptr, d_A_rows, d_A_columns, d_A_values
    );

    //--------------------------------------------------------------------------
    // Internal (basin-basin) horizontal + vertical contributions, gathered
    // per basin cell -- no atomics (see buildMatrix_cell_gather).
    launch_kernel(
	buildMatrix_cell_gather,
	idBasinVect.size(),
	stream,
	thrust::raw_pointer_cast(idBasinVect.data()),
        thrust::raw_pointer_cast(H_int_x.data()),
        thrust::raw_pointer_cast(H_int_y.data()),
        thrust::raw_pointer_cast(alfa_x.data()),
        thrust::raw_pointer_cast(alfa_y.data()),
	thrust::raw_pointer_cast(orography.data()),
	thrust::raw_pointer_cast(u_star.data()),
	thrust::raw_pointer_cast(v_star.data()),
	thrust::raw_pointer_cast(basin_mask.data()),
	thrust::raw_pointer_cast(idBasinVectReIndex.data()),
	N_cols, static_cast<unsigned int>(basin_mask.size()),
	c1, c3, H_min,
	rhs.ptr, d_A_rows, d_A_columns, d_A_values
    );

    //--------------------------------------------------------------------------
    // Boundary contributions (rhs only, matrix untouched). Two launches -- one
    // over the horizontal (internal+West+East) list, one over the vertical
    // (internal+North+South) list -- replace the four West/East/North/South
    // kernels; tag==0 entries are skipped inside the kernel.
    launch_kernel(
	buildMatrix_boundary_merged,
        idHorizontalAll.size(),
	stream,
        thrust::raw_pointer_cast(idHorizontalAll.data()),
        thrust::raw_pointer_cast(horizontalTag.data()),
	thrust::raw_pointer_cast(H_int_x.data()),
        thrust::raw_pointer_cast(alfa_x.data()),
	thrust::raw_pointer_cast(u_star.data()),
	thrust::raw_pointer_cast(idBasinVectReIndex.data()),
	N_cols, c1, H_min, isNonReflectingBC, /*isHorizontal=*/1,
	rhs.ptr
    );

    launch_kernel(
	buildMatrix_boundary_merged,
        idVerticalAll.size(),
	stream,
        thrust::raw_pointer_cast(idVerticalAll.data()),
        thrust::raw_pointer_cast(verticalTag.data()),
	thrust::raw_pointer_cast(H_int_y.data()),
        thrust::raw_pointer_cast(alfa_y.data()),
	thrust::raw_pointer_cast(v_star.data()),
	thrust::raw_pointer_cast(idBasinVectReIndex.data()),
	N_cols, c1, H_min, isNonReflectingBC, /*isHorizontal=*/0,
	rhs.ptr
    );

}
 
//==============================================================================
// https://github.com/NVIDIA/cuda-samples
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
	   const bool           use_preconditioner = true) {

    const double zero      = 0.0;
    const double one       = 1.0;
    const double minus_one = -1.0;

    // ### 1 ### R0 = b - A * X0
    CHECK_CUDA( cudaMemcpy(d_R.ptr, d_B.ptr, m * sizeof(double),
                           cudaMemcpyDeviceToDevice) )
    CHECK_CUSPARSE( cusparseSpMV(cusparseHandle,
                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &minus_one, matA, d_X.vec, &one, d_R.vec,
                                 CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT,
                                 d_bufferMV) )

    // ### 2 ### R_aux_0 = L^-1 L^-T R_0
    if (use_preconditioner) {
        // IC solve: L^-1 L^-T R. No memset of d_tmp/d_R_aux needed first:
        // cusparseSpSV_solve is a full triangular solve, every element of
        // the output is written, none are accumulated onto prior content.
	CHECK_CUSPARSE( cusparseSpSV_solve(cusparseHandle,
				CUSPARSE_OPERATION_NON_TRANSPOSE,
				&one, matL, d_R.vec, d_tmp.vec,
				CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT,
				spsvDescrL) )
	CHECK_CUSPARSE( cusparseSpSV_solve(cusparseHandle,
				CUSPARSE_OPERATION_TRANSPOSE,
				&one, matL, d_tmp.vec, d_R_aux.vec,
				CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT,
				spsvDescrLT) )
    } else {
        // no preconditioner: R_aux = R
        CHECK_CUDA( cudaMemcpy(d_R_aux.ptr, d_R.ptr, m * sizeof(double),
                               cudaMemcpyDeviceToDevice) )
    }

    // ### 3 ### P0 = R_aux_0
    CHECK_CUDA( cudaMemcpy(d_P.ptr, d_R_aux.ptr, m * sizeof(double),
                           cudaMemcpyDeviceToDevice) )

    double nrm_R;
    CHECK_CUBLAS( cublasDnrm2(cublasHandle, m, d_R.ptr, 1, &nrm_R) )
    double threshold = tolerance * nrm_R;

    double delta;
    CHECK_CUBLAS( cublasDdot(cublasHandle, m, d_R.ptr, 1, d_R_aux.ptr, 1, &delta) )

    // ### 4 ### CG iterations
    int iterations_used = maxIterations;
    for (int i = 0; i < maxIterations; i++) {

        // ### 5 ### alpha = delta / (A*P, P)
        CHECK_CUSPARSE( cusparseSpMV(cusparseHandle,
                                     CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                                     matA, d_P.vec, &zero, d_T.vec,
                                     CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT,
                                     d_bufferMV) )
        double denominator;
        CHECK_CUBLAS( cublasDdot(cublasHandle, m, d_T.ptr, 1, d_P.ptr, 1,
                                 &denominator) )
        double alpha = delta / denominator;

        // ### 6 ### X_i+1 = X_i + alpha * P
        CHECK_CUBLAS( cublasDaxpy(cublasHandle, m, &alpha, d_P.ptr, 1,
                                  d_X.ptr, 1) )

        // ### 7 ### R_i+1 = R_i - alpha * T
        double minus_alpha = -alpha;
        CHECK_CUBLAS( cublasDaxpy(cublasHandle, m, &minus_alpha, d_T.ptr, 1,
                                  d_R.ptr, 1) )

        // ### 8 ### convergence check
        CHECK_CUBLAS( cublasDnrm2(cublasHandle, m, d_R.ptr, 1, &nrm_R) )
        if (nrm_R < threshold) {
            iterations_used = i + 1;
            break;
        }

        // ### 9 ### R_aux_i+1 = L^-1 L^-T R_i+1 (no memset needed, see above)
	if (use_preconditioner) {
	    CHECK_CUSPARSE( cusparseSpSV_solve(cusparseHandle,
				CUSPARSE_OPERATION_NON_TRANSPOSE,
				&one, matL, d_R.vec, d_tmp.vec,
				CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT,
				spsvDescrL) )
	    CHECK_CUSPARSE( cusparseSpSV_solve(cusparseHandle,
				CUSPARSE_OPERATION_TRANSPOSE,
				&one, matL, d_tmp.vec, d_R_aux.vec,
				CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT,
				spsvDescrLT) )
	}
	else {
	    // no preconditioner: R_aux = R
            CHECK_CUDA( cudaMemcpy(d_R_aux.ptr, d_R.ptr, m * sizeof(double),
                                   cudaMemcpyDeviceToDevice) )
	}

        // ### 10 ### beta = delta_new / delta
        double delta_new;
        CHECK_CUBLAS( cublasDdot(cublasHandle, m, d_R.ptr, 1, d_R_aux.ptr, 1,
                                 &delta_new) )
        double beta = delta_new / delta;
        delta       = delta_new;

        // ### 11 ### P_i+1 = R_aux_i+1 + beta * P_i
        CHECK_CUBLAS( cublasDscal(cublasHandle, m, &beta, d_P.ptr, 1) )
        CHECK_CUBLAS( cublasDaxpy(cublasHandle, m, &one, d_R_aux.ptr, 1,
                                  d_P.ptr, 1) )
    }

    printf("GPU CG: iterations performed: %d, residual achieved: %e\n",
           iterations_used, nrm_R);

    return EXIT_SUCCESS;
}

//==============================================================================

__global__ void updateH_kernel(
    const unsigned int* ids,
    const unsigned int* idBasinVectReIndex,
    double* H,
    double* eta,
    const double* orography,
    const double* d_X,
    unsigned int n)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    const unsigned int Id        = ids[i];
    const unsigned int IDreIndex = idBasinVectReIndex[Id];

    H[Id]   = fabs(d_X[IDreIndex]);
    eta[Id] = H[Id] + orography[Id];
}

void updateH_wrapper(
    const thrust::device_vector<unsigned int>& idBasinVect,
    const thrust::device_vector<unsigned int>& idBasinVectReIndex,
    thrust::device_vector<double>& H,
    thrust::device_vector<double>& eta,
    const thrust::device_vector<double>& orography,
    const double* d_X,
    cudaStream_t stream)
{
    launch_kernel(updateH_kernel, idBasinVect.size(), stream,
        thrust::raw_pointer_cast(idBasinVect.data()),
        thrust::raw_pointer_cast(idBasinVectReIndex.data()),
        thrust::raw_pointer_cast(H.data()),
        thrust::raw_pointer_cast(eta.data()),
        thrust::raw_pointer_cast(orography.data()),
        d_X);
}

//==============================================================================
// updateVel kernels

// Merged internal+West+East horizontal update: one launch instead of three.
// `tag[i]` (0=internal, 1=West boundary, 2=East boundary) says which of the
// three original id lists ids[i] came from -- set once at setup by directly
// copying list membership, not re-derived here, so this carries no risk of
// misclassifying a face (see the tag-construction comment in main_final_H.cpp
// for why: the boundary lists mix index conventions in a way that isn't safe
// to re-infer per-face from basin_mask alone). Each branch below is the
// verbatim formula of the kernel it replaces.
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
    unsigned int n)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    const unsigned int Id = ids[i];
    const unsigned int ii = Id / (N_cols + 1);
    const unsigned int t  = tag[i];

    if (t == 0) {
        // ---- internal ----
        const unsigned int IDeast = Id - ii;
        const unsigned int IDwest = Id - ii - 1;
        const double H_e = H[IDeast], H_w = H[IDwest];
        const double eta_diff = -eta[IDeast] + eta[IDwest];
        const double s = (eta_diff > 0.0) - (eta_diff < 0.0);
        const double H_interface = (H_e + H_w) * 0.5 + s * (-H_e + H_w) * 0.5;
        u[Id] = (H_interface > H_min)
                    ? alfa_x[Id] * (u_star[Id] - c2 * (eta[IDeast] - eta[IDwest]))
                    : 0.0;
    } else if (t == 1) {
        // ---- West boundary ----
        const unsigned int IDeast = Id - ii + 1;
        const unsigned int IDwest = Id - ii;
        const double H_w = H[IDwest];
        const double eta_diff = -eta[IDeast] + eta[IDwest];
        const double s = (eta_diff > 0.0) - (eta_diff < 0.0);
        const double H_interface = H_w * 0.5 + s * (-H_w) * 0.5;
        u[Id] = (H_interface > H_min) ? isNonReflectingBC * alfa_x[Id] * u_star[Id] : 0.0;
    } else {
        // ---- East boundary ----
        const unsigned int IDeast = Id - ii - 1;
        const unsigned int IDwest = Id - ii - 2;
        const double H_e = H[IDeast];
        const double eta_diff = -eta[IDeast] + eta[IDwest];
        const double s = (eta_diff > 0.0) - (eta_diff < 0.0);
        const double H_interface = H_e * 0.5 + s * H_e * 0.5;
        u[Id] = (H_interface > H_min) ? isNonReflectingBC * alfa_x[Id] * u_star[Id] : 0.0;
    }
}

// Merged internal+North+South vertical update, same rationale as above.
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
    unsigned int n)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    const unsigned int Id = ids[i];
    const unsigned int t  = tag[i];

    if (t == 0) {
        // ---- internal ----
        const unsigned int IDsouth = Id;
        const unsigned int IDnorth = Id - N_cols;
        const double H_s = H[IDsouth], H_n = H[IDnorth];
        const double eta_diff = -eta[IDsouth] + eta[IDnorth];
        const double s = (eta_diff > 0.0) - (eta_diff < 0.0);
        const double H_interface = (H_s + H_n) * 0.5 + s * (-H_s + H_n) * 0.5;
        v[Id] = (H_interface > H_min)
                    ? alfa_y[Id] * (v_star[Id] - c2 * (eta[IDsouth] - eta[IDnorth]))
                    : 0.0;
    } else if (t == 1) {
        // ---- North boundary ----
        const unsigned int IDnorth = Id;
        const unsigned int IDsouth = Id + N_cols;
        const double H_n = H[IDnorth];
        const double eta_diff = -eta[IDsouth] + eta[IDnorth];
        const double s = (eta_diff > 0.0) - (eta_diff < 0.0);
        const double H_interface = H_n * 0.5 + s * (-H_n) * 0.5;
        v[Id] = (H_interface > H_min) ? isNonReflectingBC * alfa_y[Id] * v_star[Id] : 0.0;
    } else {
        // ---- South boundary ----
        const unsigned int IDsouth = Id - N_cols;
        const unsigned int IDnorth = Id - 2 * N_cols;
        const double H_s = H[IDsouth];
        const double eta_diff = -eta[IDsouth] + eta[IDnorth];
        const double s = (eta_diff > 0.0) - (eta_diff < 0.0);
        const double H_interface = H_s * 0.5 + s * H_s * 0.5;
        v[Id] = (H_interface > H_min) ? isNonReflectingBC * alfa_y[Id] * v_star[Id] : 0.0;
    }
}

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
    cudaStream_t stream)
{
    const double* H_ptr   = thrust::raw_pointer_cast(H.data());
    const double* eta_ptr = thrust::raw_pointer_cast(eta.data());
    double* u_ptr         = thrust::raw_pointer_cast(u.data());
    double* v_ptr         = thrust::raw_pointer_cast(v.data());
    const double* us_ptr  = thrust::raw_pointer_cast(u_star.data());
    const double* vs_ptr  = thrust::raw_pointer_cast(v_star.data());
    const double* ax_ptr  = thrust::raw_pointer_cast(alfa_x.data());
    const double* ay_ptr  = thrust::raw_pointer_cast(alfa_y.data());

    launch_kernel(updateVelVerticalMerged,
        idVerticalAll.size(), stream,
        thrust::raw_pointer_cast(idVerticalAll.data()),
        thrust::raw_pointer_cast(verticalTag.data()),
        v_ptr, vs_ptr, ay_ptr, H_ptr, eta_ptr, c2, H_min, isNonReflectingBC, N_cols);

    launch_kernel(updateVelHorizontalMerged,
        idHorizontalAll.size(), stream,
        thrust::raw_pointer_cast(idHorizontalAll.data()),
        thrust::raw_pointer_cast(horizontalTag.data()),
        u_ptr, us_ptr, ax_ptr, H_ptr, eta_ptr, c2, H_min, isNonReflectingBC, N_cols);
}

//==============================================================================

#if 0 == 1
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
