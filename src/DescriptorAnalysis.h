/**
 * This file is part of lpf
 *
 * @file
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
#ifndef DESCRIPTOR_ANALYSIS_H
#define DESCRIPTOR_ANALYSIS_H

#include "LPF.h"
#include "Eigen/Dense"
#include <vector>

using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::VectorXd;

class DescriptorAnalysis
{
public :
  
  DescriptorAnalysis();
  
  template<class T>
  DescriptorAnalysis(const std::vector<T*> &descs)
  {
    unsigned int sizedesc = descs[0]->getSizeDesc();
    unsigned int ndesc = descs.size();
    m_X.resize(sizedesc, ndesc);

    typename std::vector<T*>::const_iterator vi;
    size_t j = 0;
    for(vi = descs.begin(); vi != descs.end(); ++vi, ++j)
    {
      for(size_t i = 0 ; i < sizedesc ; ++i)
      {
        double val = 0;
        (*vi)->vectorize(i, val);
        m_X(i,j) = val;
      }
    }
  }
  
  ~DescriptorAnalysis();
  
  bool learnLPFdictionaryAndFramesICP(std::vector<LPF*> &descs, 
                       MatrixXd &D,
                       MatrixXd &alpha,
                       MatrixXd &t,
                       std::vector<Vector3d> &t1s,
                       std::vector<Vector3d> &t2s,
                       std::vector<Vector3d> &ns,
                       unsigned int natoms,
                       double lambda,
                       unsigned int dic_update_niter,
                       unsigned int niter);
  
  
  bool setGenericPattern(std::vector<Vector3d> &genpat);
  
  MatrixXd & getX();
  
private :
  
  MatrixXd m_X;
  
  VectorXd m_genpat; //generic pattern
  
  
  bool recompute();
  
  
  
  bool optimizeLPFbyICP(MatrixXd& V,
                        MatrixXd& T,
                        MatrixXd& tloc,
                        MatrixXd &R,
                        std::vector<Vector3d> &t1s,
                        std::vector<Vector3d> &t2s,
                        std::vector<Vector3d> &ns);
  
  
  bool tovector(const Matrix3d &X, VectorXd &v) const;
  
  bool tovector(const MatrixXd &X, VectorXd &v) const;
  
  bool tomatrix(const VectorXd &v, Matrix3d &X) const;
  
  bool tomatrix(const VectorXd &v, MatrixXd &X) const;
 
  template<class T>
  bool recompute(std::vector<T*> &descs)
  {
#pragma omp parallel for
    for(int i = 0; i < descs.size(); ++i)
    {
      //descs[i]->update(m_X.col(i));
      descs[i]->recompute(m_genpat);
      VectorXd v;
      descs[i]->tovector(v);
      m_X.col(i) = v;
    }
    return true;
  }
  
  template<class T>
  bool update_target_zones(std::vector<T*> &descs, const MatrixXd &R, const MatrixXd &Tr)
  {
#pragma omp parallel for
    for(int i = 0; i < descs.size(); ++i)
    {
      Matrix3d Ri;
      Vector3d Ti;
      tomatrix(R.row(i), Ri);
      Ti = Tr.row(i);
      descs[i]->update_target_zone_after_ICP(Ri, Ti);
    }
    return true;
  }
};

#endif
