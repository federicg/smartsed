#include "code_init.h"

//==============================================================================

int computePourCell(const int &IDcell, const unsigned int &N_cols,
                    const std::vector<double> &oro,
                    const std::set<unsigned int> &idBasinVect,
                    const std::set<unsigned int> &idStaggeredBoundaryVectSouth,
                    const std::set<unsigned int> &idStaggeredBoundaryVectNorth,
                    const std::set<unsigned int> &idStaggeredBoundaryVectWest,
                    const std::set<unsigned int> &idStaggeredBoundaryVectEast) {

  int candidate_id = -1;

  // +-----------------------------------------------+
  // |               "Cartesian" cells               |
  // +-----------------------------------------------+

  const int i = IDcell / N_cols;

  const int IDsouth = IDcell + N_cols, IDnorth = IDcell, IDwest = IDcell + i,
            IDeast = IDcell + i + 1;

  const auto iterator_south = idStaggeredBoundaryVectSouth.find(IDsouth),
             iterator_north = idStaggeredBoundaryVectNorth.find(IDnorth),
             iterator_west = idStaggeredBoundaryVectWest.find(IDwest),
             iterator_east = idStaggeredBoundaryVectEast.find(IDeast);

  std::vector<unsigned int> candidates;
  candidates.reserve(8);

  bool is_north = false, is_south = false, is_west = false, is_east = false;

  if (iterator_south == idStaggeredBoundaryVectSouth.end()) {
    const int IDcell_south = IDcell + N_cols;
    candidates.push_back(IDcell_south);

    is_south = true;
  }

  if (iterator_north == idStaggeredBoundaryVectNorth.end()) {
    const int IDcell_north = IDcell - N_cols;
    candidates.push_back(IDcell_north);

    is_north = true;
  }

  if (iterator_west == idStaggeredBoundaryVectWest.end()) {
    const int IDcell_west = IDcell - 1;
    candidates.push_back(IDcell_west);

    is_west = true;
  }

  if (iterator_east == idStaggeredBoundaryVectEast.end()) {
    const int IDcell_east = IDcell + 1;
    candidates.push_back(IDcell_east);

    is_east = true;
  }

  // +-----------------------------------------------+
  // |               Diagonal cells                  |
  // +-----------------------------------------------+

  if (is_north && is_west) {
    // exist north-west cell for sure
    const int IDcel_nw = IDcell - 1 - N_cols;

    if (idBasinVect.find(IDcel_nw) != idBasinVect.end()) {
      candidates.push_back(IDcel_nw);
    }
  }

  if (is_north && is_east) {
    // exist north-east cell for sure
    const int IDcel_ne = IDcell + 1 - N_cols;

    if (idBasinVect.find(IDcel_ne) != idBasinVect.end()) {
      candidates.push_back(IDcel_ne);
    }
  }

  if (is_south && is_west) {
    // exist south-west cell for sure
    const int IDcel_sw = IDcell - 1 + N_cols;

    if (idBasinVect.find(IDcel_sw) != idBasinVect.end()) {
      candidates.push_back(IDcel_sw);
    }
  }

  if (is_south && is_east) {
    // exist south-east cell for sure
    const int IDcel_se = IDcell + 1 + N_cols;

    if (idBasinVect.find(IDcel_se) != idBasinVect.end()) {
      candidates.push_back(IDcel_se);
    }
  }

  // +-----------------------------------------------+
  // |       Compute min oro cell (pour point)       |
  // +-----------------------------------------------+

  double current_minimum_oro = oro[IDcell];
  for (const auto &it : candidates) {
    const auto &candidate_minimum_oro = oro[it];
    if (candidate_minimum_oro < current_minimum_oro) {
      current_minimum_oro = candidate_minimum_oro;
      candidate_id = it;
    }
  }

  return candidate_id;
}

//==============================================================================

void computeAdjacencies(
    const std::vector<double> &basin_mask_Vec,

    std::vector<unsigned int> &idStaggeredBoundaryVectSouth,
    std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    std::vector<unsigned int> &idStaggeredBoundaryVectEast,

    std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    std::vector<unsigned int> &idStaggeredInternalVectVertical,

    std::vector<unsigned int> &idBasinVect,
    std::vector<unsigned int> &idBasinVectReIndex,

    const unsigned int &N_rows, const unsigned int &N_cols) {

  // +-----------------------------------------------+
  // |                 Basin H IDs                   |
  // +-----------------------------------------------+

  unsigned int h = 0;
  for (unsigned int i = 0; i < N_rows; i++) {
    for (unsigned int j = 0; j < N_cols; j++) {
      const unsigned int k = j + i * N_cols;
      idBasinVectReIndex.push_back(h);

      if (basin_mask_Vec[k] == 1) {
        idBasinVect.push_back(k);
        h++;
      }
    }
  }

  // +-----------------------------------------------+
  // |         Vertical Vel. Staggered IDs           |
  // +-----------------------------------------------+

  // cycle on centered cells
  for (unsigned int i = 0; i < N_rows; i++) {
    for (unsigned int j = 0; j < N_cols; j++) {
      const unsigned int
          IDcell = j + i * N_cols,
          IDcell_south = IDcell + N_cols,
          IDvel = IDcell + N_cols; // interface between IDcell and IDcell_south

      if (i != (N_rows - 1)) {
        if ((basin_mask_Vec[IDcell] + basin_mask_Vec[IDcell_south]) ==
            1) // interface cell
        {
          if (basin_mask_Vec[IDcell] == 0) {
            idStaggeredBoundaryVectNorth.push_back(IDvel);
          } else {
            idStaggeredBoundaryVectSouth.push_back(IDvel);
          }
        }

        if (i == 0 && basin_mask_Vec[IDcell] == 1) {
          const auto IDvel_north = IDcell;

          idStaggeredBoundaryVectNorth.push_back(IDvel_north); // it is ok
        }

        if ((basin_mask_Vec[IDcell] + basin_mask_Vec[IDcell_south]) == 2) {
          idStaggeredInternalVectVertical.push_back(
              IDvel); // it is ok no repetition
        }

      } else {
        if (basin_mask_Vec[IDcell] == 1) {
          const auto IDvel_south = IDvel;

          idStaggeredBoundaryVectSouth.push_back(IDvel_south);
        }
      }
    }
  }

  // +-----------------------------------------------+
  // |         Horizontal Vel. Staggered IDs         |
  // +-----------------------------------------------+

  // cycle on centered cells
  for (unsigned int i = 0; i < N_rows; i++) {
    for (unsigned int j = 0; j < N_cols; j++) {
      const unsigned int IDcell = j + i * N_cols, IDcell_east = IDcell + 1,
                         IDvel = IDcell + i +
                                 1; // interface between IDcell and IDcell_east

      if (j != (N_cols - 1)) {
        if ((basin_mask_Vec[IDcell] + basin_mask_Vec[IDcell_east]) == 1) {
          if (basin_mask_Vec[IDcell] == 0) {
            idStaggeredBoundaryVectWest.push_back(IDvel);
          } else {
            idStaggeredBoundaryVectEast.push_back(IDvel);
          }
        }

        if (j == 0 && basin_mask_Vec[IDcell] == 1) {
          const auto IDvel_west = IDcell + i;

          idStaggeredBoundaryVectWest.push_back(IDvel_west);
        }

        if ((basin_mask_Vec[IDcell] + basin_mask_Vec[IDcell_east]) == 2) {
          idStaggeredInternalVectHorizontal.push_back(IDvel);
        }
      } else {
        if (basin_mask_Vec[IDcell] == 1) {
          const auto IDvel_east = IDvel;
          idStaggeredBoundaryVectEast.push_back(IDvel_east);
        }
      }
    }
  }
}

//==============================================================================

void computeAdjacencies(
    const std::vector<double> &basin_mask_Vec_input,
    const std::vector<std::tuple<bool, int>> &excluded_ids,

    std::vector<unsigned int> &idStaggeredBoundaryVectSouth,
    std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    std::vector<unsigned int> &idStaggeredBoundaryVectEast,

    std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    std::vector<unsigned int> &idStaggeredInternalVectVertical,

    std::vector<unsigned int> &idBasinVect,
    std::vector<unsigned int> &idBasinVectReIndex,

    const unsigned int &N_rows, const unsigned int &N_cols) {

  std::vector<double> basin_mask_Vec = basin_mask_Vec_input;

  for (unsigned int i = 0; i < basin_mask_Vec.size(); i++) {
    if (std::get<0>(excluded_ids[i])) {
      basin_mask_Vec[i] = 0.;
    }
  }

  // +-----------------------------------------------+
  // |                 Basin H IDs                   |
  // +-----------------------------------------------+

  unsigned int h = 0;
  for (unsigned int i = 0; i < N_rows; i++) {
    for (unsigned int j = 0; j < N_cols; j++) {
      const unsigned int k = j + i * N_cols;

      idBasinVectReIndex.push_back(h);

      if (basin_mask_Vec[k] == 1) {
        idBasinVect.push_back(k);
        h++;
      }
    }
  }

  // +-----------------------------------------------+
  // |         Vertical Vel. Staggered IDs           |
  // +-----------------------------------------------+

  // cycle on centered cells
  for (unsigned int i = 0; i < N_rows; i++) {
    for (unsigned int j = 0; j < N_cols; j++) {
      const unsigned int
          IDcell = j + i * N_cols,
          IDcell_south = IDcell + N_cols,
          IDvel = IDcell + N_cols; // interface between IDcell and IDcell_south

      if (i != (N_rows - 1)) {
        if ((basin_mask_Vec[IDcell] + basin_mask_Vec[IDcell_south]) ==
            1) // interface cell
        {
          if (basin_mask_Vec[IDcell] == 0) {
            idStaggeredBoundaryVectNorth.push_back(IDvel);
          } else {
            idStaggeredBoundaryVectSouth.push_back(IDvel);
          }
        }

        if (i == 0 && basin_mask_Vec[IDcell] == 1) {
          const auto IDvel_north = IDcell;

          idStaggeredBoundaryVectNorth.push_back(IDvel_north); // it is ok
        }

        if ((basin_mask_Vec[IDcell] + basin_mask_Vec[IDcell_south]) == 2) {
          idStaggeredInternalVectVertical.push_back(
              IDvel); // it is ok no repetition
        }

      } else {

        if (basin_mask_Vec[IDcell] == 1) {
          const auto IDvel_south = IDvel;

          idStaggeredBoundaryVectSouth.push_back(IDvel_south);
        }
      }
    }
  }

  // +-----------------------------------------------+
  // |         Horizontal Vel. Staggered IDs         |
  // +-----------------------------------------------+

  // cicle on centered cells
  for (unsigned int i = 0; i < N_rows; i++) {
    for (unsigned int j = 0; j < N_cols; j++) {
      const unsigned int IDcell = j + i * N_cols, IDcell_east = IDcell + 1,
                         IDvel = IDcell + i +
                                 1; // interface between IDcell and IDcell_east

      if (j != (N_cols - 1)) {
        if ((basin_mask_Vec[IDcell] + basin_mask_Vec[IDcell_east]) == 1) {
          if (basin_mask_Vec[IDcell] == 0) {
            idStaggeredBoundaryVectWest.push_back(IDvel);
          } else {
            idStaggeredBoundaryVectEast.push_back(IDvel);
          }
        }

        if (j == 0 && basin_mask_Vec[IDcell] == 1) {
          const auto IDvel_west = IDcell + i;

          idStaggeredBoundaryVectWest.push_back(IDvel_west);
        }

        if ((basin_mask_Vec[IDcell] + basin_mask_Vec[IDcell_east]) == 2) {
          idStaggeredInternalVectHorizontal.push_back(IDvel);
        }
      } else {
        if (basin_mask_Vec[IDcell] == 1) {
          const auto IDvel_east = IDvel;

          idStaggeredBoundaryVectEast.push_back(IDvel_east);
        }
      }
    }
  }
}

//==============================================================================

void saveVector(const Eigen::VectorXd &b, const std::string &Name) {
  std::ofstream ff(Name);
  for (unsigned int k = 0; k < b.size(); k++) {
    ff << b(k) << " ";
    ff << std::endl;
  }
  ff.close();
}

//==============================================================================

void saveMatrix(const Eigen::SparseMatrix<double> &A, const std::string &Name) {
  std::ofstream ff(Name);
  for (unsigned int k = 0; k < A.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
      ff << it.row() + 1 << " " << it.col() + 1 << " " << it.value()
         << std::endl; // row index
    }
  }
  ff.close();
}

//==============================================================================

void saveTemporalSequence(const Vector2D &X_gauges, const double &time,
                          const std::string &preName, const double &H) {
  std::ofstream ff(preName + ".txt", std::ofstream::out | std::ofstream::app);
  ff << H << " " << time << std::endl;
  ff.close();
}

//==============================================================================

void saveTemporalSequence(const double &time, const std::string &preName,
                          const double &H) {
  std::ofstream ff(preName + ".txt", std::ofstream::out | std::ofstream::app);
  ff << H << " " << time << std::endl;
  ff.close();
}

//==============================================================================

void putDry_excludedNodes(
    const std::vector<unsigned int> &idStaggeredInternalVectHorizontal,
    const std::vector<unsigned int> &idStaggeredInternalVectVertical,
    const std::vector<unsigned int> &idStaggeredBoundaryVectWest,
    const std::vector<unsigned int> &idStaggeredBoundaryVectEast,
    const std::vector<unsigned int> &idStaggeredBoundaryVectNorth,
    const std::vector<unsigned int> &idStaggeredBoundaryVectSouth,
    const std::vector<unsigned int> &idBasinVect, const unsigned int &N_cols,
    const std::vector<std::tuple<bool, int>> &excluded_ids,

    Eigen::VectorXd &H, Eigen::VectorXd &eta,
    const std::vector<double> &orography, std::vector<double> &u,
    std::vector<double> &v) {

  for (const auto &Id : idBasinVect) {
    const auto &current_tuple = excluded_ids[Id];
    if (std::get<0>(current_tuple)) {
      H(Id) = 0.;
      eta(Id) = orography[Id];
    }
  }

  for (const auto &Id : idStaggeredInternalVectHorizontal) {
    const unsigned int i = Id / (N_cols + 1),

                       IDleft = Id - i - 1, // H
        IDright = Id - i;                   // H

    if (std::get<0>(excluded_ids[IDleft]) &&
        std::get<0>(excluded_ids[IDright])) // se entrambe le celle sono escluse
    {
      H(IDleft) = 0.;
      H(IDright) = 0.;
      u[Id] = 0.;
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectWest) {
    const unsigned int i = Id / (N_cols + 1), IDright = Id - i;

    if (std::get<0>(excluded_ids[IDright])) {
      H(IDright) = 0.;
      u[Id] = 0.;
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectEast) {
    const unsigned int i = Id / (N_cols + 1), IDleft = Id - i - 1;

    if (std::get<0>(excluded_ids[IDleft])) {
      H(IDleft) = 0.;
      u[Id] = 0.;
    }
  }

  for (const auto &Id : idStaggeredInternalVectVertical) {
    const unsigned int IDleft = Id - N_cols, // H
        IDright = Id;                        // H

    if (std::get<0>(excluded_ids[IDleft]) &&
        std::get<0>(excluded_ids[IDright])) // se entrambe le celle sono escluse
    {
      H(IDleft) = 0.;
      H(IDright) = 0.;
      v[Id] = 0.;
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectNorth) {
    const unsigned int IDright = Id;

    if (std::get<0>(excluded_ids[IDright])) // se entrambe le celle sono escluse
    {
      H(IDright) = 0.;
      v[Id] = 0.;
    }
  }

  for (const auto &Id : idStaggeredBoundaryVectSouth) {
    const unsigned int IDleft = Id - N_cols;

    if (std::get<0>(excluded_ids[IDleft])) // se entrambe le celle sono escluse
    {
      H(IDleft) = 0.;
      v[Id] = 0.;
    }
  }
}

//==============================================================================

Raster::Raster(const std::string &file) {
  std::cout << "Reading file, " << file << std::endl;

  std::vector<Eigen::Triplet<double>> cc;

  std::ifstream ff(file);
  if (ff.is_open()) {
    std::string str;

    ff >> str;
    if (str != "ncols" && str != "NCOLS") {
      std::cout << "Wrong file format in the Raster .txt files" << std::endl;
      exit(-1.);
    }
    ff >> ncols;

    ff >> str;
    if (str != "nrows" && str != "NROWS") {
      std::cout << "Wrong file format in the Raster .txt files" << std::endl;
      exit(-1.);
    }
    ff >> nrows;

    ff >> str;
    if (str != "xllcorner" && str != "XLLCORNER") {
      std::cout << "Wrong file format in the Raster .txt files" << std::endl;
      exit(-1.);
    }
    ff >> xllcorner;

    ff >> str;
    if (str != "yllcorner" && str != "YLLCORNER") {
      std::cout << "Wrong file format in the Raster .txt files" << std::endl;
      exit(-1.);
    }
    ff >> yllcorner;

    ff >> str;
    if (str != "cellsize" && str != "CELLSIZE") {
      std::cout << "Wrong file format in the Raster .txt files" << std::endl;
      exit(-1.);
    }
    ff >> cellsize;

    ff >> str;
    if (str != "NODATA_value") {
      std::cout << "Wrong file format in the Raster .txt files" << std::endl;
      exit(-1.);
    }
    ff >> NODATA_value;

    cc.reserve(nrows * ncols);

    double value;
    for (unsigned int i = 0; i < nrows; i++) {
      for (unsigned int j = 0; j < ncols; j++) {

        ff >> value;

        cc.push_back(Eigen::Triplet<double>(i, j, value));
      }
    }
  } else {
    std::cout << "Unable to open the file, check " << file
              << " in the SMARTSED_input file" << std::endl;
    exit(-1.);
  }

  Coords.resize(nrows, ncols);
  Coords.setFromTriplets(cc.begin(), cc.end());
};

//==============================================================================

Vector2D operator*(Vector2D const &vector, double const &factor) {
  Vector2D tmp(vector);
  return tmp *= factor;
}

//==============================================================================

Vector2D operator*(double const &factor, Vector2D const &vector) {
  Vector2D tmp(vector);
  return tmp *= factor;
}

//==============================================================================

std::map<int, std::array<double, 2>> createCN_map_Gav(const std::string &file) {
  std::cout << "Reading Gavrilovic coefficients, ... " << std::endl;

  std::map<int, std::array<double, 2>> CN;

  std::vector<std::array<double, 2>> vector_Gav(44);

  std::ifstream ff(file);
  if (ff.is_open()) {
    double X, Y;
    for (unsigned int i = 0; i < 44; i++) {
      ff >> X;
      ff >> Y;
      vector_Gav[i] = std::array<double, 2>{{X, Y}};
    }

  } else {
    std::cout << "Unable to open the file, check " << file
              << " in the SMARTSED_input file" << std::endl;
    exit(-1.);
  }

  CN[111] = std::array<double, 2>{{vector_Gav[0][0], vector_Gav[0][1]}};
  CN[112] = std::array<double, 2>{{vector_Gav[1][0], vector_Gav[1][1]}};
  CN[121] = std::array<double, 2>{{vector_Gav[2][0], vector_Gav[2][1]}};
  CN[122] = std::array<double, 2>{{vector_Gav[3][0], vector_Gav[3][1]}};
  CN[123] = std::array<double, 2>{{vector_Gav[4][0], vector_Gav[4][1]}};
  CN[124] = std::array<double, 2>{{vector_Gav[5][0], vector_Gav[5][1]}};
  CN[131] = std::array<double, 2>{{vector_Gav[6][0], vector_Gav[6][1]}};
  CN[132] = std::array<double, 2>{{vector_Gav[7][0], vector_Gav[7][1]}};
  CN[133] = std::array<double, 2>{{vector_Gav[8][0], vector_Gav[8][1]}};
  CN[141] = std::array<double, 2>{{vector_Gav[9][0], vector_Gav[9][1]}};
  CN[142] = std::array<double, 2>{{vector_Gav[10][0], vector_Gav[10][1]}};
  CN[211] = std::array<double, 2>{{vector_Gav[11][0], vector_Gav[11][1]}};
  CN[212] = std::array<double, 2>{{vector_Gav[12][0], vector_Gav[12][1]}};
  CN[213] = std::array<double, 2>{{vector_Gav[13][0], vector_Gav[13][1]}};
  CN[221] = std::array<double, 2>{{vector_Gav[14][0], vector_Gav[14][1]}};
  CN[222] = std::array<double, 2>{{vector_Gav[15][0], vector_Gav[15][1]}};
  CN[223] = std::array<double, 2>{{vector_Gav[16][0], vector_Gav[16][1]}};
  CN[231] = std::array<double, 2>{{vector_Gav[17][0], vector_Gav[17][1]}};
  CN[241] = std::array<double, 2>{{vector_Gav[18][0], vector_Gav[18][1]}};
  CN[242] = std::array<double, 2>{{vector_Gav[19][0], vector_Gav[19][1]}};
  CN[243] = std::array<double, 2>{{vector_Gav[20][0], vector_Gav[20][1]}};
  CN[244] = std::array<double, 2>{{vector_Gav[21][0], vector_Gav[21][1]}};
  CN[311] = std::array<double, 2>{{vector_Gav[22][0], vector_Gav[22][1]}};
  CN[312] = std::array<double, 2>{{vector_Gav[23][0], vector_Gav[23][1]}};
  CN[313] = std::array<double, 2>{{vector_Gav[24][0], vector_Gav[24][1]}};
  CN[321] = std::array<double, 2>{{vector_Gav[25][0], vector_Gav[25][1]}};
  CN[322] = std::array<double, 2>{{vector_Gav[26][0], vector_Gav[26][1]}};
  CN[323] = std::array<double, 2>{{vector_Gav[27][0], vector_Gav[27][1]}};
  CN[324] = std::array<double, 2>{{vector_Gav[28][0], vector_Gav[28][1]}};
  CN[331] = std::array<double, 2>{{vector_Gav[29][0], vector_Gav[29][1]}};
  CN[332] = std::array<double, 2>{{vector_Gav[30][0], vector_Gav[30][1]}};
  CN[333] = std::array<double, 2>{{vector_Gav[31][0], vector_Gav[31][1]}};
  CN[334] = std::array<double, 2>{{vector_Gav[32][0], vector_Gav[32][1]}};
  CN[335] = std::array<double, 2>{{vector_Gav[33][0], vector_Gav[33][1]}};
  CN[411] = std::array<double, 2>{{vector_Gav[34][0], vector_Gav[34][1]}};
  CN[412] = std::array<double, 2>{{vector_Gav[35][0], vector_Gav[35][1]}};
  CN[421] = std::array<double, 2>{{vector_Gav[36][0], vector_Gav[36][1]}};
  CN[422] = std::array<double, 2>{{vector_Gav[37][0], vector_Gav[37][1]}};
  CN[423] = std::array<double, 2>{{vector_Gav[38][0], vector_Gav[38][1]}};
  CN[511] = std::array<double, 2>{{vector_Gav[39][0], vector_Gav[39][1]}};
  CN[512] = std::array<double, 2>{{vector_Gav[40][0], vector_Gav[40][1]}};
  CN[521] = std::array<double, 2>{{vector_Gav[41][0], vector_Gav[41][1]}};
  CN[522] = std::array<double, 2>{{vector_Gav[42][0], vector_Gav[42][1]}};
  CN[523] = std::array<double, 2>{{vector_Gav[43][0], vector_Gav[43][1]}};

  return CN;
}

//==============================================================================

std::map<std::array<int, 2>, int> createCN_map() {
  std::cout << "Reading CN coefficients, ... " << std::endl;

  std::map<std::array<int, 2>, int> CN;

  CN[std::array<int, 2>{{111, 0}}] = 89;
  CN[std::array<int, 2>{{111, 1}}] = 92;
  CN[std::array<int, 2>{{111, 2}}] = 94;
  CN[std::array<int, 2>{{111, 3}}] = 95;

  CN[std::array<int, 2>{{112, 0}}] = 77;
  CN[std::array<int, 2>{{112, 1}}] = 85;
  CN[std::array<int, 2>{{112, 2}}] = 90;
  CN[std::array<int, 2>{{112, 3}}] = 92;

  CN[std::array<int, 2>{{121, 0}}] = 81;
  CN[std::array<int, 2>{{121, 1}}] = 88;
  CN[std::array<int, 2>{{121, 2}}] = 91;
  CN[std::array<int, 2>{{121, 3}}] = 93;

  CN[std::array<int, 2>{{122, 0}}] = 83;
  CN[std::array<int, 2>{{122, 1}}] = 89;
  CN[std::array<int, 2>{{122, 2}}] = 92;
  CN[std::array<int, 2>{{122, 3}}] = 93;

  CN[std::array<int, 2>{{123, 0}}] = 83;
  CN[std::array<int, 2>{{123, 1}}] = 89;
  CN[std::array<int, 2>{{123, 2}}] = 92;
  CN[std::array<int, 2>{{123, 3}}] = 93;

  CN[std::array<int, 2>{{124, 0}}] = 83;
  CN[std::array<int, 2>{{124, 1}}] = 89;
  CN[std::array<int, 2>{{124, 2}}] = 92;
  CN[std::array<int, 2>{{124, 3}}] = 93;

  CN[std::array<int, 2>{{131, 0}}] = 81;
  CN[std::array<int, 2>{{131, 1}}] = 88;
  CN[std::array<int, 2>{{131, 2}}] = 91;
  CN[std::array<int, 2>{{131, 3}}] = 93;

  CN[std::array<int, 2>{{132, 0}}] = 81;
  CN[std::array<int, 2>{{132, 1}}] = 88;
  CN[std::array<int, 2>{{132, 2}}] = 91;
  CN[std::array<int, 2>{{132, 3}}] = 93;

  CN[std::array<int, 2>{{133, 0}}] = 77;
  CN[std::array<int, 2>{{133, 1}}] = 86;
  CN[std::array<int, 2>{{133, 2}}] = 91;
  CN[std::array<int, 2>{{133, 3}}] = 94;

  CN[std::array<int, 2>{{141, 0}}] = 49;
  CN[std::array<int, 2>{{141, 1}}] = 69;
  CN[std::array<int, 2>{{141, 2}}] = 79;
  CN[std::array<int, 2>{{141, 3}}] = 84;

  CN[std::array<int, 2>{{142, 0}}] = 49;
  CN[std::array<int, 2>{{142, 1}}] = 69;
  CN[std::array<int, 2>{{142, 2}}] = 79;
  CN[std::array<int, 2>{{142, 3}}] = 84;

  CN[std::array<int, 2>{{211, 0}}] = 60;
  CN[std::array<int, 2>{{211, 1}}] = 72;
  CN[std::array<int, 2>{{211, 2}}] = 80;
  CN[std::array<int, 2>{{211, 3}}] = 84;

  CN[std::array<int, 2>{{212, 0}}] = 60;
  CN[std::array<int, 2>{{212, 1}}] = 72;
  CN[std::array<int, 2>{{212, 2}}] = 80;
  CN[std::array<int, 2>{{212, 3}}] = 84;

  CN[std::array<int, 2>{{213, 0}}] = 62;
  CN[std::array<int, 2>{{213, 1}}] = 71;
  CN[std::array<int, 2>{{213, 2}}] = 78;
  CN[std::array<int, 2>{{213, 3}}] = 81;

  CN[std::array<int, 2>{{221, 0}}] = 43;
  CN[std::array<int, 2>{{221, 1}}] = 65;
  CN[std::array<int, 2>{{221, 2}}] = 76;
  CN[std::array<int, 2>{{221, 3}}] = 82;

  CN[std::array<int, 2>{{222, 0}}] = 43;
  CN[std::array<int, 2>{{222, 1}}] = 65;
  CN[std::array<int, 2>{{222, 2}}] = 76;
  CN[std::array<int, 2>{{222, 3}}] = 82;

  CN[std::array<int, 2>{{223, 0}}] = 43;
  CN[std::array<int, 2>{{223, 1}}] = 65;
  CN[std::array<int, 2>{{223, 2}}] = 76;
  CN[std::array<int, 2>{{223, 3}}] = 82;

  CN[std::array<int, 2>{{231, 0}}] = 30;
  CN[std::array<int, 2>{{231, 1}}] = 58;
  CN[std::array<int, 2>{{231, 2}}] = 71;
  CN[std::array<int, 2>{{231, 3}}] = 78;

  CN[std::array<int, 2>{{241, 0}}] = 58;
  CN[std::array<int, 2>{{241, 1}}] = 72;
  CN[std::array<int, 2>{{241, 2}}] = 81;
  CN[std::array<int, 2>{{241, 3}}] = 85;

  CN[std::array<int, 2>{{242, 0}}] = 59;
  CN[std::array<int, 2>{{242, 1}}] = 74;
  CN[std::array<int, 2>{{242, 2}}] = 82;
  CN[std::array<int, 2>{{242, 3}}] = 86;

  CN[std::array<int, 2>{{243, 0}}] = 59;
  CN[std::array<int, 2>{{243, 1}}] = 74;
  CN[std::array<int, 2>{{243, 2}}] = 82;
  CN[std::array<int, 2>{{243, 3}}] = 86;

  CN[std::array<int, 2>{{244, 0}}] = 43;
  CN[std::array<int, 2>{{244, 1}}] = 65;
  CN[std::array<int, 2>{{244, 2}}] = 76;
  CN[std::array<int, 2>{{244, 3}}] = 82;

  CN[std::array<int, 2>{{311, 0}}] = 36;
  CN[std::array<int, 2>{{311, 1}}] = 60;
  CN[std::array<int, 2>{{311, 2}}] = 73;
  CN[std::array<int, 2>{{311, 3}}] = 79;

  CN[std::array<int, 2>{{312, 0}}] = 36;
  CN[std::array<int, 2>{{312, 1}}] = 60;
  CN[std::array<int, 2>{{312, 2}}] = 73;
  CN[std::array<int, 2>{{312, 3}}] = 79;

  CN[std::array<int, 2>{{313, 0}}] = 36;
  CN[std::array<int, 2>{{313, 1}}] = 60;
  CN[std::array<int, 2>{{313, 2}}] = 73;
  CN[std::array<int, 2>{{313, 3}}] = 79;

  CN[std::array<int, 2>{{321, 0}}] = 39;
  CN[std::array<int, 2>{{321, 1}}] = 61;
  CN[std::array<int, 2>{{321, 2}}] = 74;
  CN[std::array<int, 2>{{321, 3}}] = 80;

  CN[std::array<int, 2>{{322, 0}}] = 77;
  CN[std::array<int, 2>{{322, 1}}] = 86;
  CN[std::array<int, 2>{{322, 2}}] = 91;
  CN[std::array<int, 2>{{322, 3}}] = 94;

  CN[std::array<int, 2>{{323, 0}}] = 77;
  CN[std::array<int, 2>{{323, 1}}] = 86;
  CN[std::array<int, 2>{{323, 2}}] = 91;
  CN[std::array<int, 2>{{323, 3}}] = 94;

  CN[std::array<int, 2>{{324, 0}}] = 35;
  CN[std::array<int, 2>{{324, 1}}] = 56;
  CN[std::array<int, 2>{{324, 2}}] = 70;
  CN[std::array<int, 2>{{324, 3}}] = 77;

  CN[std::array<int, 2>{{331, 0}}] = 55;
  CN[std::array<int, 2>{{331, 1}}] = 72;
  CN[std::array<int, 2>{{331, 2}}] = 81;
  CN[std::array<int, 2>{{331, 3}}] = 86;

  CN[std::array<int, 2>{{332, 0}}] = 77;
  CN[std::array<int, 2>{{332, 1}}] = 86;
  CN[std::array<int, 2>{{332, 2}}] = 91;
  CN[std::array<int, 2>{{332, 3}}] = 94;

  CN[std::array<int, 2>{{333, 0}}] = 74;
  CN[std::array<int, 2>{{333, 1}}] = 83;
  CN[std::array<int, 2>{{333, 2}}] = 88;
  CN[std::array<int, 2>{{333, 3}}] = 90;

  CN[std::array<int, 2>{{334, 0}}] = 77;
  CN[std::array<int, 2>{{334, 1}}] = 86;
  CN[std::array<int, 2>{{334, 2}}] = 91;
  CN[std::array<int, 2>{{334, 3}}] = 94;

  CN[std::array<int, 2>{{335, 0}}] = 77;
  CN[std::array<int, 2>{{335, 1}}] = 86;
  CN[std::array<int, 2>{{335, 2}}] = 91;
  CN[std::array<int, 2>{{335, 3}}] = 94;

  CN[std::array<int, 2>{{411, 0}}] = 98;
  CN[std::array<int, 2>{{411, 1}}] = 98;
  CN[std::array<int, 2>{{411, 2}}] = 98;
  CN[std::array<int, 2>{{411, 3}}] = 98;

  CN[std::array<int, 2>{{412, 0}}] = 98;
  CN[std::array<int, 2>{{412, 1}}] = 98;
  CN[std::array<int, 2>{{412, 2}}] = 98;
  CN[std::array<int, 2>{{412, 3}}] = 98;

  CN[std::array<int, 2>{{421, 0}}] = 98;
  CN[std::array<int, 2>{{421, 1}}] = 98;
  CN[std::array<int, 2>{{421, 2}}] = 98;
  CN[std::array<int, 2>{{421, 3}}] = 98;

  CN[std::array<int, 2>{{422, 0}}] = 98;
  CN[std::array<int, 2>{{422, 1}}] = 98;
  CN[std::array<int, 2>{{422, 2}}] = 98;
  CN[std::array<int, 2>{{422, 3}}] = 98;

  CN[std::array<int, 2>{{423, 0}}] = 98;
  CN[std::array<int, 2>{{423, 1}}] = 98;
  CN[std::array<int, 2>{{423, 2}}] = 98;
  CN[std::array<int, 2>{{423, 3}}] = 98;

  CN[std::array<int, 2>{{511, 0}}] = 100;
  CN[std::array<int, 2>{{511, 1}}] = 100;
  CN[std::array<int, 2>{{511, 2}}] = 100;
  CN[std::array<int, 2>{{511, 3}}] = 100;

  CN[std::array<int, 2>{{512, 0}}] = 100;
  CN[std::array<int, 2>{{512, 1}}] = 100;
  CN[std::array<int, 2>{{512, 2}}] = 100;
  CN[std::array<int, 2>{{512, 3}}] = 100;

  CN[std::array<int, 2>{{521, 0}}] = 100;
  CN[std::array<int, 2>{{521, 1}}] = 100;
  CN[std::array<int, 2>{{521, 2}}] = 100;
  CN[std::array<int, 2>{{521, 3}}] = 100;

  CN[std::array<int, 2>{{522, 0}}] = 100;
  CN[std::array<int, 2>{{522, 1}}] = 100;
  CN[std::array<int, 2>{{522, 2}}] = 100;
  CN[std::array<int, 2>{{522, 3}}] = 100;

  CN[std::array<int, 2>{{523, 0}}] = 100;
  CN[std::array<int, 2>{{523, 1}}] = 100;
  CN[std::array<int, 2>{{523, 2}}] = 100;
  CN[std::array<int, 2>{{523, 3}}] = 100;

  return CN;
}

//==============================================================================

bool is_file_exist(const char *fileName) {
  std::ifstream infile(fileName);
  return infile.good();
}

//==============================================================================

std::vector<double> compute_d_perc(const std::vector<double> &clay,
                                   const std::vector<double> &sand,
                                   const double &perc) {
  // linear interpolation in log10 x-scale
  std::vector<double> d_perc(clay.size());

  for (unsigned int i = 0; i < clay.size(); i++) {
    auto &d_perc_cell = d_perc[i];

    const auto &sand_cell = sand[i];
    const auto &clay_cell = clay[i];

    const auto Y_0 = (1 - sand_cell) * 100;
    const auto Y_1 = clay_cell * 100; // ( 1 - sand_cell - silt_cell ) * 100
    const auto Y_2 = 100;

    if (perc <= Y_0) {
      const double angular_coeff =
          (Y_0 - Y_1) / (std::log10(25.)); // 50 \mu m / 2 \mu m
      const double DY = perc - Y_0;

      d_perc_cell = 50.e-6 * std::pow(10, DY / angular_coeff);
    } else {
      const double angular_coeff =
          (Y_2 - Y_0) / (std::log10(40)); // 2 mm / 50 \mu m
      const double DY = perc - Y_0;

      d_perc_cell = 50.e-6 * std::pow(10, DY / angular_coeff);
    }
  }
  return d_perc;
}

void saveSolution(const std::string &preName, const std::string &flag,
                  const unsigned int &N_rows, const unsigned int &N_cols,
                  const double &xllcorner, const double &yllcorner,
                  const double &cellsize, const double &NODATA_value,
                  const Eigen::VectorXd &H) // it is H or orography
{

  std::ofstream ff(preName + ".asc");

  if (flag == "u") {
    ff << "ncols ";
    ff << N_cols + 1;
    ff << std::endl;

    ff << "nrows ";
    ff << N_rows;
    ff << std::endl;
  } else if (flag == "v") {
    ff << "ncols ";
    ff << N_cols;
    ff << std::endl;

    ff << "nrows ";
    ff << N_rows + 1;
    ff << std::endl;
  } else {
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

  if (flag == "u") {

    for (unsigned int i = 0; i < N_rows; i++) {

      for (unsigned int j = 0; j <= N_cols; j++) {

        const auto Id = j + i * (N_cols + 1);

        ff << H[Id] << " ";
      }

      ff << std::endl;
    }

  } else if (flag == "v") {

    for (unsigned int i = 0; i <= N_rows; i++) {

      for (unsigned int j = 0; j < N_cols; j++) {

        const auto Id = j + i * N_cols;

        ff << H[Id] << " ";
      }

      ff << std::endl;
    }

  } else // H or orography
  {

    for (unsigned int i = 0; i < N_rows; i++) {

      for (unsigned int j = 0; j < N_cols; j++) {

        const auto k = j + i * N_cols; // H

        ff << H[k] << " ";
      }

      ff << std::endl;
    }
  }

  ff.close();
}

//==============================================================================

void saveSolution(const std::string &preName, const unsigned int &N_rows,
                  const unsigned int &N_cols, const double &xllcorner,
                  const double &yllcorner, const double &cellsize,
                  const double &NODATA_value,
                  const std::vector<std::tuple<bool, int>>
                      excluded_ids) // excluded regions, high slopes I hope
{

  {
    std::ofstream ff(preName + "_bool" + ".asc");

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

    for (unsigned int i = 0; i < N_rows; i++) {

      for (unsigned int j = 0; j < N_cols; j++) {

        const auto k = j + i * N_cols; // H

        ff << std::get<0>(excluded_ids[k]) << " ";
      }

      ff << std::endl;
    }

    ff.close();
  }

  std::ofstream ff(preName + "_pour" + ".asc");

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

  for (unsigned int i = 0; i < N_rows; i++) {

    for (unsigned int j = 0; j < N_cols; j++) {

      const auto k = j + i * N_cols; // H

      ff << std::get<1>(excluded_ids[k]) << " ";
    }

    ff << std::endl;
  }

  ff.close();
}

//==============================================================================

void saveSolution(const std::string &preName, const std::string &flag,
                  const unsigned int &N_rows, const unsigned int &N_cols,
                  const double &xllcorner, const double &yllcorner,
                  const double &cellsize, const double &NODATA_value,
                  const std::vector<double> &H) // it is H or orography
{

  std::ofstream ff(preName + ".asc");

  if (flag == "u") {
    ff << "ncols ";
    ff << N_cols + 1;
    ff << std::endl;

    ff << "nrows ";
    ff << N_rows;
    ff << std::endl;
  } else if (flag == "v") {
    ff << "ncols ";
    ff << N_cols;
    ff << std::endl;

    ff << "nrows ";
    ff << N_rows + 1;
    ff << std::endl;
  } else {
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

  if (flag == "u") {

    for (unsigned int i = 0; i < N_rows; i++) {

      for (unsigned int j = 0; j <= N_cols; j++) {

        const auto Id = j + i * (N_cols + 1);

        ff << H[Id] << " ";
      }

      ff << std::endl;
    }

  } else if (flag == "v") {

    for (unsigned int i = 0; i <= N_rows; i++) {

      for (unsigned int j = 0; j < N_cols; j++) {

        const auto Id = j + i * N_cols;

        ff << H[Id] << " ";
      }

      ff << std::endl;
    }

  } else // H or orography
  {

    for (unsigned int i = 0; i < N_rows; i++) {

      for (unsigned int j = 0; j < N_cols; j++) {

        const auto k = j + i * N_cols; // H

        ff << H[k] << " ";
      }

      ff << std::endl;
    }
  }

  ff.close();
}

//==============================================================================

void saveSolution(const std::string &preName, const std::string &flag,
                  const unsigned int &N_rows, const unsigned int &N_cols,
                  const double &xllcorner, const double &yllcorner,
                  const double &cellsize, const double &NODATA_value,
                  const std::vector<int> &H) // it is H or orography
{

  std::ofstream ff(preName + ".asc");

  if (flag == "u") {
    ff << "ncols ";
    ff << N_cols + 1;
    ff << std::endl;

    ff << "nrows ";
    ff << N_rows;
    ff << std::endl;
  } else if (flag == "v") {
    ff << "ncols ";
    ff << N_cols;
    ff << std::endl;

    ff << "nrows ";
    ff << N_rows + 1;
    ff << std::endl;
  } else {
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

  if (flag == "u") {

    for (unsigned int i = 0; i < N_rows; i++) {

      for (unsigned int j = 0; j <= N_cols; j++) {

        const auto Id = j + i * (N_cols + 1);

        ff << H[Id] << " ";
      }

      ff << std::endl;
    }

  } else if (flag == "v") {

    for (unsigned int i = 0; i <= N_rows; i++) {

      for (unsigned int j = 0; j < N_cols; j++) {

        const auto Id = j + i * N_cols;

        ff << H[Id] << " ";
      }

      ff << std::endl;
    }

  } else // H or orography
  {

    for (unsigned int i = 0; i < N_rows; i++) {

      for (unsigned int j = 0; j < N_cols; j++) {

        const auto k = j + i * N_cols; // H

        ff << H[k] << " ";
      }

      ff << std::endl;
    }
  }

  ff.close();
}

//==============================================================================

void saveSolution(const std::string &preName, const std::string &flag,
                  const unsigned int &N_rows, const unsigned int &N_cols,
                  const double &xllcorner, const double &yllcorner,
                  const double &cellsize, const double &NODATA_value,
                  const unsigned int &n, const std::vector<double> &u,
                  const std::vector<double> &v,
                  const Eigen::VectorXd &H) // it is H or orography
{

  std::ofstream ff(preName + std::to_string(n) + ".asc");

  if (flag == "u") {
    ff << "ncols ";
    ff << N_cols + 1;
    ff << std::endl;

    ff << "nrows ";
    ff << N_rows;
    ff << std::endl;
  } else if (flag == "v") {
    ff << "ncols ";
    ff << N_cols;
    ff << std::endl;

    ff << "nrows ";
    ff << N_rows + 1;
    ff << std::endl;
  } else {
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

  if (flag == "u") {

    for (unsigned int i = 0; i < N_rows; i++) {

      for (unsigned int j = 0; j <= N_cols; j++) {

        const auto Id = j + i * (N_cols + 1);

        ff << u[Id] << " ";
      }

      ff << std::endl;
    }

  } else if (flag == "v") {

    for (unsigned int i = 0; i <= N_rows; i++) {

      for (unsigned int j = 0; j < N_cols; j++) {

        const auto Id = j + i * N_cols;

        ff << v[Id] << " ";
      }

      ff << std::endl;
    }

  } else // H or orography
  {

    for (unsigned int i = 0; i < N_rows; i++) {

      for (unsigned int j = 0; j < N_cols; j++) {

        const auto k = j + i * N_cols; // H

        ff << H(k) << " ";
      }

      ff << std::endl;
    }
  }

  ff.close();
}

//==============================================================================

void saveSolution(const std::string &preName, const std::string &flag,
                  const unsigned int &N_rows, const unsigned int &N_cols,
                  const double &xllcorner, const double &yllcorner,
                  const double &cellsize, const double &NODATA_value,
                  const unsigned int &n, const std::vector<double> &u,
                  const std::vector<double> &v,
                  const std::vector<double> &H) // it is H or orography
{

  std::ofstream ff(preName + std::to_string(n) + ".asc");

  if (flag == "u") {
    ff << "ncols ";
    ff << N_cols + 1;
    ff << std::endl;

    ff << "nrows ";
    ff << N_rows;
    ff << std::endl;
  } else if (flag == "v") {
    ff << "ncols ";
    ff << N_cols;
    ff << std::endl;

    ff << "nrows ";
    ff << N_rows + 1;
    ff << std::endl;
  } else {
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

  if (flag == "u") {

    for (unsigned int i = 0; i < N_rows; i++) {

      for (unsigned int j = 0; j <= N_cols; j++) {

        const auto Id = j + i * (N_cols + 1);

        ff << u[Id] << " ";
      }

      ff << std::endl;
    }

  } else if (flag == "v") {

    for (unsigned int i = 0; i <= N_rows; i++) {

      for (unsigned int j = 0; j < N_cols; j++) {

        const auto Id = j + i * N_cols;

        ff << v[Id] << " ";
      }

      ff << std::endl;
    }

  } else // H or orography
  {

    for (unsigned int i = 0; i < N_rows; i++) {

      for (unsigned int j = 0; j < N_cols; j++) {

        const auto k = j + i * N_cols; // H

        ff << H[k] << " ";
      }

      ff << std::endl;
    }
  }

  ff.close();
}

//==============================================================================
