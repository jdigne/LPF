/**
 * This file is part of lpf
 *
 * @file LPF.h
 * @author Julie Digne
 *
 * Copyright (c) 2015-2020 Julie Digne
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */             
#ifndef LPF_H
#define LPF_H

#include <iostream>
#include <fstream>
#include<vector>

#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;

class LPF
{

public :

  LPF();
  
  LPF(unsigned int nbins);
  
  LPF(std::vector<Vector3d> &qpoints, std::vector<Vector3d> &locpoints);
  
  ~LPF();
  
  void optimize_efficiency(std::vector<Vector3d> &locpoints,
                           const MatrixXd &genpat,
                           Vector3d &seed,
                           Vector3d &t1,
                           Vector3d &t2,
                           Vector3d &n,
                           double radius);
  

  unsigned int getSizeDesc() const;//get the number of double to store a vectorized descriptors
  
  bool vectorize(size_t i, double &val) const;
  
  bool write(std::ofstream &f);
  
  bool writeBin(std::ofstream &f);
  
  bool read(unsigned int nbins, std::ifstream &f);

  bool readBin(unsigned int nbins, std::ifstream &f);
  
  bool readMatrix(const VectorXd &X);
  
  
  bool reconstruct(const Vector3d &seed,
                   const Vector3d &t1,
                   const Vector3d &t2,
                   const Vector3d &n,
                   const double radius,
                   const std::vector<Vector3d> &qpoints,
                   std::vector<Vector3d> &recpoints);
  
  
  bool save_target_zone(const Vector3d &seed,
                   const Vector3d &t1,
                   const Vector3d &t2,
                   const Vector3d &n,
                   const double radius,
                   std::vector<Vector3d> &recpoints);
  
  void update(const VectorXd &val);
  
  void tovector(VectorXd &val) const;
  
  void recompute(const VectorXd &qpoints);
  
  void recompute(const MatrixXd &qpoints);

  double vote(Vector3d &point,
              const Vector3d seed,
              const Vector3d &t1,
              const Vector3d &t2,
              const Vector3d &normal,
              const std::vector<Vector3d> &qpoints,
              const double radius);
  
  void update_target_zone(Eigen::Matrix3d &R, Eigen::Vector3d &t);
  
  void update_target_zone_after_ICP(Eigen::Matrix3d &R, Eigen::Vector3d &t);
  
private :

  unsigned int m_nbins;

  bool m_irreg;
  
  Eigen::MatrixXd m_v;

  double parallelism_measure();
  
  std::vector<Vector3d> m_target;
};

#endif
