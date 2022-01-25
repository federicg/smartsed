#ifndef UTILS_H_H
#  define UTILS_H_H

#include "typedefs_H.h"
#include <array>

#include <omp.h>

// you have to declare both to make it work
#include <mpi.h>
#include <lis.h>


class
Vector2D
{
  
public:

  //! Empty constructor (all components are set to zero)
  Vector2D () {}

  Vector2D (std::array<Real,2> const& indices) // Note the use of Real indices
    : M_coords (indices) {}

  //! Copy constructor
  Vector2D
  (Vector2D const& vector)
  { *this = vector; }

  ~Vector2D() = default;

  //! Operator +=
  Vector2D&
  operator+= (Vector2D const& vector)
  {
    for ( UInt i = 0; i < 2; i++ )
      M_coords[ i ] += vector.M_coords[ i ];
    return *this;
  }

  //! Assignment operator
  Vector2D&
  operator= (Vector2D const& vector)
  {
    for (UInt i = 0; i < 2; i++)
      M_coords[ i ] = vector.M_coords[ i ];
    return *this;
  }

  Vector2D
  operator+ (Vector2D const& vector) const
  {
    Vector2D tmp( *this );
    return tmp += vector;
  }

  //! Operator -=
  Vector2D&
  operator-= (Vector2D const& vector)
  {
    for (UInt i = 0; i < 2; i++)
      M_coords[ i ] -= vector.M_coords[ i ];
    return *this;
  }

  //! Operator -
  Vector2D
  operator- (Vector2D const& vector) const
  {
    Vector2D tmp( *this );
    return tmp -= vector;
  }

  //! Operator *= (multiplication by scalar)
  Vector2D&
  operator*= ( Real const& factor )
  {
    for ( UInt i = 0; i < 2; i++ )
      M_coords[ i ] *= factor;
    return *this;
  }

  //! Operator /= (division by scalar)
  Vector2D&
  operator/= ( Real const& factor )
  {
    *this *= 1. / factor;
    return *this;
  }

  //! Operator / (division by scalar)
  Vector2D
  operator/ ( Real const& factor ) const
  {
    Vector2D tmp (*this);
    return tmp /= factor;
  }


  Vector2D
  operator* (Real const& factor)
  {
    Vector2D tmp( *this );
    return tmp *= factor;
  }


  Real
  dot ( Vector2D const& vector ) const
  {
    Real scalarProduct = 0.;
    for ( UInt i = 0; i < 2; i++ )
      scalarProduct += M_coords[ i ] * vector.M_coords[ i ];
    return scalarProduct;
  }

  Real norm () const
  { return std::sqrt (this->dot (*this)); }


  //! Operator ()b
  Real const& operator() ( UInt const& i ) const
  { return M_coords[ i ]; }

  //! Operator ()
  Real& operator() ( UInt const& i )
  { return M_coords[ i ]; }

private:

  std::array<Real,2> M_coords;

};

//! Operator * (multiplication by scalar on the right)
Vector2D
operator* (Vector2D const& vector, Real const& factor);

//! Operator * (multiplication by scalar on the left)
Vector2D
operator* (Real const& factor, Vector2D const& vector);

std::map<Int, std::array<Real,2>>
createCN_map_Gav ();

std::map<std::array<Int,2>, Int>
createCN_map ();


class
Raster
{

public:

  Raster (const std::string& file, const int& rank);

  ~Raster() = default;

  UInt ncols, nrows;

  Real xllcorner, yllcorner,
    cellsize, NODATA_value;

  SpMat Coords;  // forse mettere una matrice densa

};


Real
signum (const Real& x);

class Rain // Previous interpolation not Linear
{

public:

  Rain( const std::string& infiltrationModel,
        const UInt&        N,
        const bool&        isInitialLoss,
        const Real&        perc_initialLoss );

  Rain() = delete;
  ~Rain() = default;

  void
  constant_precipitation (const std::string&       file,
                          const UInt&              ndata,
                          const bool&              is_precipitation,
                          const Real&              time_spacing);


  void
  IDW_precipitation (const std::vector<std::string>& file_vect,
                     const std::vector<UInt>&        ndata_vec,
                     const std::vector<Real>&        time_spacing_vect,
                     const std::vector<Real>&        X,
                     const std::vector<Real>&        Y,
                     const Real&                     xllcorner,
                     const Real&                     yllcorner,
                     const Real&                     pixel_size,
                     const UInt&                     N_rows,
                     const UInt&                     N_cols,
                     const std::vector<UInt>&        idBasinVect);

  void
  computePrecipitation (const Real&              time,
                        const std::vector<Real>& S,
                        const std::vector<Real>& melt_mask,
                        const std::vector<Real>& h_G,
                        const std::vector<Real>& H,
                        const UInt&              N_rows,
                        const UInt&              N_cols,
                        const std::vector<UInt>& idBasinVect);


  std::vector<Real> DP_total,  DP_cumulative,  DP_infiltrated;

private:

  std::vector<std::vector<Real> > Hyetograph, // # station times ndata
    IDW_weights;

  std::vector<Real> M_time_spacing_vect;
  bool M_isInitialLoss;
  Real rainfall_intensity = 0;
  Real c;

};


class
Temperature
{

public:

  Temperature (const std::string&       file,
               const UInt&              N,
               const UInt&              max_Days,
               const Real&              T_crit,
               const std::vector<Real>& orography,
               const UInt&              ndata,
               const UInt&              steps_per_hour,
               const Real&              time_spacing,
               const Real&              height_thermometer,
               const std::string        format_temp );

  Temperature() = delete;
  ~Temperature() = default;


  void
  computeTemperature (const UInt&              i,
                      const std::vector<Real>& orography,
                      const std::vector<UInt>& idBasinVect);

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

class
evapoTranspiration
{

public:

  evapoTranspiration (const std::string&       ET_model,
                      const UInt&              N,
                      const std::vector<Real>& orography,
                      const std::vector<Real>& J,
                      const UInt&              max_Days,
                      const Real&              phi_rad,
                      const Real&              height_thermometer);

  evapoTranspiration() = delete;
  ~evapoTranspiration() = default;

  void ET (const std::vector<Real>& T_mean, // length nstep: vector of temperature in °C
           const std::vector<Real>& T_min,  // length nstep
           const std::vector<Real>& T_max,  // length nstep
           const Int&               i,       
           const std::vector<UInt>& idBasinVect,
           const std::vector<Real>& orography);


  std::vector<Real> ET_vec;

private:

  std::vector<Real> Ra;
  UInt M_evapoTranspiration_model;
  static constexpr Real M_Gsc = .082; // Solar constant
  const Real height_th;
  static constexpr Real Temp_diff = -6.5e-3;

};


class frictionClass
{

public:

  frictionClass (const std::vector<Real>& H_interface_horizontal, 
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
                 const std::vector<Real>& S_y);

  frictionClass() = delete;
  ~frictionClass() = default;

  void
  f_x ();

  void
  f_y ();

  std::vector<Real> alfa_x,
    alfa_y;

private:

  UInt M_frictionModel;

  const Real& M_n_manning;
  const Real& M_dt_DSV;

  Real M_coeff;
  std::function<Real(Real const&, Real const&)> M_gamma_dt_DSV = [](Real const& dt, Real const& cc){return dt * cc;};
  
  

  Real M_H_min;

  const std::vector<Real>& u; 
  const std::vector<Real>& v;
  const std::vector<Real>& H_interface_horizontal;
  const std::vector<Real>& H_interface_vertical;

  const std::vector<UInt>& idStaggeredInternalVectHorizontal;
  const std::vector<UInt>& idStaggeredBoundaryVectWest;
  const std::vector<UInt>& idStaggeredBoundaryVectEast;
  const std::vector<UInt>& idStaggeredInternalVectVertical;
  const std::vector<UInt>& idStaggeredBoundaryVectNorth;
  const std::vector<UInt>& idStaggeredBoundaryVectSouth;

  const UInt& N_rows;
  const UInt& N_cols;

  static constexpr Real M_expo    = 4./3.;
  static constexpr Real M_expo_r1 = .11 * 2;
  static constexpr Real M_expo_r2 = .03 * 2;
  static constexpr Real M_g       = 9.81;
  static constexpr Real M_toll    = 1.e-4;

  std::vector<Real> M_expo_r_x_vect,
    M_expo_r_y_vect,
    M_gamma_dt_DSV_x_,
    M_gamma_dt_DSV_y_;

};

class
upwind
{

public:

  upwind (const std::vector<Real>& H,
          const std::vector<Real>& u,
          const std::vector<Real>& v,

          const std::vector<UInt>& idStaggeredInternalVectHorizontal,
          const std::vector<UInt>& idStaggeredBoundaryVectWest,
          const std::vector<UInt>& idStaggeredBoundaryVectEast,
          const std::vector<UInt>& idStaggeredInternalVectVertical,
          const std::vector<UInt>& idStaggeredBoundaryVectNorth,
          const std::vector<UInt>& idStaggeredBoundaryVectSouth,

          const UInt& N_rows, const UInt& N_cols)
    : H(H), u(u), v(v), idStaggeredInternalVectHorizontal(idStaggeredInternalVectHorizontal),
                        idStaggeredBoundaryVectWest(idStaggeredBoundaryVectWest),
                        idStaggeredBoundaryVectEast(idStaggeredBoundaryVectEast),
                        idStaggeredInternalVectVertical(idStaggeredInternalVectVertical),
                        idStaggeredBoundaryVectNorth(idStaggeredBoundaryVectNorth),
                        idStaggeredBoundaryVectSouth(idStaggeredBoundaryVectSouth), N_cols( N_cols ), N_rows( N_rows )
  {
    horizontal.resize( u.size() );
    vertical.  resize( v.size() );
  }

  upwind() = delete;
  ~upwind() = default;

  void
  computeHorizontal ();

  void
  computeVertical ();

  std::vector<Real> horizontal,
    vertical;

private:

  const std::vector<Real>& H;
  const std::vector<Real>& u;
  const std::vector<Real>& v;

  const std::vector<UInt>& idStaggeredInternalVectHorizontal;
  const std::vector<UInt>& idStaggeredBoundaryVectWest;
  const std::vector<UInt>& idStaggeredBoundaryVectEast;
  const std::vector<UInt>& idStaggeredInternalVectVertical;
  const std::vector<UInt>& idStaggeredBoundaryVectNorth;
  const std::vector<UInt>& idStaggeredBoundaryVectSouth;

  const UInt& N_cols;
  const UInt& N_rows;

};


bool is_file_exist(const char *fileName);

void
bilinearInterpolation (const std::vector<Real>& u,
                       const std::vector<Real>& v,
                       std::vector<Real>& u_star,
                       std::vector<Real>& v_star,
                       const UInt&              nrows,
                       const UInt&              ncols,
                       const Real&              dt,
                       const Real&              pixel_size);


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
                       const std::vector<UInt>& idStaggeredBoundaryVectSouth);


int
computePourCell(const int& IDcell,
                const UInt& N_cols,
                const std::vector<Real>& oro,
                const std::set<UInt>& idBasinVect,
                const std::set<UInt>& idStaggeredBoundaryVectSouth,
                const std::set<UInt>& idStaggeredBoundaryVectNorth,
                const std::set<UInt>& idStaggeredBoundaryVectWest,
                const std::set<UInt>& idStaggeredBoundaryVectEast); 





void
computeAdjacencies (const std::vector<Real>& basin_mask_Vec,
                    const std::vector<Real>& basin_mask_Vec_mpi,

                    std::vector<UInt>& idStaggeredBoundaryVectSouth_mpi,
                    std::vector<UInt>& idStaggeredBoundaryVectNorth_mpi,
                    std::vector<UInt>& idStaggeredBoundaryVectWest_mpi,
                    std::vector<UInt>& idStaggeredBoundaryVectEast_mpi,

                    std::vector<UInt>& idStaggeredInternalVectHorizontal_mpi,
                    std::vector<UInt>& idStaggeredInternalVectVertical_mpi,

                    std::vector<UInt>& idStaggeredBoundaryVectHorizontal_among_ranks,
                    std::vector<UInt>& idStaggeredBoundaryVectVertical_among_ranks,

                    std::vector<UInt>& idStaggeredInternalVectHorizontal_in_rank,
                    std::vector<UInt>& idStaggeredInternalVectVertical_in_rank,

                    std::vector<UInt>& idBasinVect_mpi,
                    std::vector<UInt>& idBasinVectReIndex_mpi,

                    const UInt&              N_rows,
                    const UInt&              N_cols);


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
                    const UInt&              N_cols);


void
buildMatrix (const std::vector<Real>& H_int_x,
             const std::vector<Real>& H_int_y,
             const std::vector<Real>& orography,
             const std::vector<Real>& u_star,
             const std::vector<Real>& v_star,
             const std::vector<Real>& u,
             const std::vector<Real>& v,
             const std::vector<Real>& H,
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

             const std::vector<Real>& basin_mask_Vec_mpi,

             const std::vector<UInt>& idBasinVectReIndex,
             const bool&              isNonReflectingBC,
             const bool&              isH,
             
             LIS_MATRIX&    A,
             LIS_VECTOR&                    rhs);



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
           const std::vector<Real>& H,
           const std::vector<Real>& orography,
           const std::vector<UInt>& idStaggeredInternalVectHorizontal,
           const std::vector<UInt>& idStaggeredInternalVectVertical,
           const std::vector<UInt>& idStaggeredBoundaryVectWest,
           const std::vector<UInt>& idStaggeredBoundaryVectEast,
           const std::vector<UInt>& idStaggeredBoundaryVectNorth,
           const std::vector<UInt>& idStaggeredBoundaryVectSouth,
           const bool&              isNonReflectingBC);








Real
maxdt (const Real& u_abs_max,
       const Real& v_abs_max,
       const Real& gravity,
       const Real& Hmax,
       const Real& pixel_size);

Real
maxCourant (const std::vector<Real>& u,
            const std::vector<Real>& v,
            const Real&              c1);

Real
maxCourant (const std::vector<Real>& H,
            const Real&            gravity,
            const Real&            c1);

Real
compute_dt_sediment (const Real&              alpha,
                     const Real&              beta,
                     const std::vector<Real>& S_x,
                     const std::vector<Real>& S_y,
                     const std::vector<Real>& u,
                     const std::vector<Real>& v,
                     const Real&              pixel_size,
                     const Real&              dt_DSV,
                     const std::vector<UInt>& idStaggeredInternalVectHorizontal,
                     const std::vector<UInt>& idStaggeredInternalVectVertical,
                     const std::vector<UInt>& idStaggeredBoundaryVectWest,
                     const std::vector<UInt>& idStaggeredBoundaryVectEast,
                     const std::vector<UInt>& idStaggeredBoundaryVectNorth,
                     const std::vector<UInt>& idStaggeredBoundaryVectSouth);

int
current_start_chunk(const int& rank, const std::vector<int>& chunk_length_vec);

void
saveVector (const std::vector<Real>& b,
            const std::string& Name);

void
saveMatrix (const SpMat& A,
            const std::string& Name);


void
saveSolution(const std::string& preName,
             const std::string& flag,
             const UInt& N_rows,
             const UInt& N_cols,
             const Real& xllcorner,
             const Real& yllcorner,
             const Real& cellsize,
             const Real& NODATA_value,
             const std::vector<Real>& H); // it is H or orography

void
saveSolution (const std::string& preName,
              const std::string& flag,
              const UInt& N_rows,
              const UInt& N_cols,
              const Real& xllcorner,
              const Real& yllcorner,
              const Real& cellsize,
              const Real& NODATA_value,
              const std::vector<Int>& H); // it is H or orography


void
saveSolution (const std::string& preName,
              const UInt& N_rows,
              const UInt& N_cols,
              const Real& xllcorner,
              const Real& yllcorner,
              const Real& cellsize,
              const Real& NODATA_value,
              const std::vector<std::tuple<bool,int>> excluded_ids); // excluded regions, high slopes I hope



void
saveSolution(const std::string& preName,
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
             const std::vector<Real>& H); // it is H or orography

void
saveTemporalSequence (const Vector2D&    X_gauges,
                      const Real&        time,
                      const std::string& preName,
                      const Real&        H);

void
saveTemporalSequence (const Real&        time,
                      const std::string& preName,
                      const Real&        H);

void
comunicationStencil(std::vector<MPI_Request>& requests, 
               const std::vector<int>& corresponding_rank_given_id, 
                     std::vector<Real>& h,
               const std::vector<UInt>& idStaggeredBoundaryVectHorizontal_among_ranks,
               const std::vector<UInt>& idStaggeredBoundaryVectVertical_among_ranks,
               const int& rank,
               const UInt& N_cols);

// For gravitational layer
void
computeResiduals (const std::vector<Real>& n_x,
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
                              std::vector<Real>& Res_y);

// For sediment transport
void
computeResidualsTruncated (const std::vector<Real>&             u,
                           const std::vector<Real>&             v,
                           const UInt&                          N_cols,
                           const UInt&                          N_rows,
                           const UInt&                          N,
                           const Real&                          c1,
                           const std::vector<Real>&             S_x,
                           const std::vector<Real>&             S_y,
                           const Real&                          alpha,
                           const Real&                          beta,
                           const Real&                          gamma,
                           const std::vector<UInt>&             idStaggeredInternalVectHorizontal,
                           const std::vector<UInt>&             idStaggeredInternalVectVertical,
                           const std::vector<UInt>&             idStaggeredBoundaryVectWest,
                           const std::vector<UInt>&             idStaggeredBoundaryVectEast,
                           const std::vector<UInt>&             idStaggeredBoundaryVectNorth,
                           const std::vector<UInt>&             idStaggeredBoundaryVectSouth,
                           std::vector<std::array<Real, 2>>&    Gamma_x,
                           std::vector<std::array<Real, 2>>&    Gamma_y);



std::vector<Real>
compute_d_perc (const std::vector<Real>& clay, const std::vector<Real>& sand, const Real& perc);


#endif
