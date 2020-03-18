/**
 * This file is part of lpf
 *
 * @file DescriptorAnalysis.cpp
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
#include "DescriptorAnalysis.h"
#include <fstream>
#include <Eigen/SVD> 
#include <cassert>
#include "utils.h"
#include "Optimization.h"

using Eigen::ArrayXd;
using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::JacobiSVD;


bool checkNAN(const MatrixXd &M)
{
  for(size_t i = 0; i < (unsigned int) M.rows(); ++i)
    for(size_t j = 0; j < (unsigned int)M.cols(); ++j)
      if(std::isnan(M(i,j)))
        return true;
  return false;    
}

DescriptorAnalysis::DescriptorAnalysis()
{

}


DescriptorAnalysis::~DescriptorAnalysis()
{

}

MatrixXd& DescriptorAnalysis::getX()
{
  return m_X;
}



bool DescriptorAnalysis::learnLPFdictionaryAndFramesICP(std::vector<LPF*> &descs, 
                                                        MatrixXd &D,
                                                        MatrixXd &alpha,
                                                        MatrixXd &t,
                                                        std::vector<Vector3d> &t1s,
                                                        std::vector<Vector3d> &t2s,
                                                        std::vector<Vector3d> &ns,
                                                        unsigned int natoms,
                                                        double lambda,
                                                        unsigned int dic_update_niter,
                                                        unsigned int niter)
{
  std::cout<<"Optimizing the dictionary, positions and orientations of the LPFs with recomputation step"<<std::endl;
  t.resize(m_X.cols(),3);
  t.setZero();
  
  MatrixXd tloc(m_X.cols(),3);
  tloc.setZero();
  
  MatrixXd R(m_X.cols(),9);
  R.setZero();
  
  Optimization optim;
  
  if(D.cols() == 0)
  {
    D.resize(m_X.rows(), natoms);
    optim.initializeRandomDictionary(m_X, D);
  }
  
  optim.omp(m_X, D, alpha, lambda);
  int batchsize = m_X.cols()/20;
  
  MatrixXd V(m_X.rows(), m_X.cols());
  
  std::ofstream f;
  f.open("error.txt");

  double error = (0.5 * (m_X - D * alpha).squaredNorm() + lambda * alpha.lpNorm<1>())/m_X.cols();
  f << error << std::endl;
  
  for(size_t iter = 0; iter < niter; ++iter)
  {
    std::cout << "iteration " << iter << std::endl;
    optim.mod_l1(m_X, D, alpha, lambda, dic_update_niter);
    
    double l = lambda * alpha.lpNorm<1>();
    
    V = D * alpha;
    error = (0.5 * (m_X - V).squaredNorm() + l) / m_X.cols();
    std::cout << "Dico optim error: " << error << std::endl;
    
    f << error << std::endl;
 
    optimizeLPFbyICP(V, t, tloc, R, t1s, t2s, ns);
    
    error = (0.5 * (m_X - V).squaredNorm() + l) / m_X.cols();
    std::cout << "ICP optim error: " << error << std::endl;
    f << error << std::endl;
    
    //debug
    update_target_zones(descs, R, tloc);
    recompute(descs);
    
    error = (0.5 * (m_X - V).squaredNorm() + l) / m_X.cols();
    std::cout << "Recomputation error: " << error << std::endl;
    f << error << std::endl;
  }
  f.close();
  
  return true;
}

bool DescriptorAnalysis::tomatrix(const VectorXd& v, Matrix3d& X) const
{
   if(X.rows()*X.cols() != v.rows())
   {
     std::cerr<<"bad vector to matrix conversion"<<std::endl; 
     return false;
   }
   
   size_t ncols = X.cols();
   size_t nrows = X.rows();
   
   for(size_t i = 0 ; i < ncols ; ++i)
       X.col(i) = v.segment(i*nrows, nrows);
   
   return true;
}

bool DescriptorAnalysis::tomatrix(const VectorXd& v, MatrixXd& X) const
{
   if(X.rows()*X.cols() != v.rows())
   {
     std::cerr<<"bad vector to matrix conversion"<<std::endl; 
     return false;
   }
   size_t ncols = X.cols();
   size_t nrows = X.rows();
   
   for(size_t i = 0 ; i < ncols ; ++i)
       X.col(i) = v.segment(i*nrows, nrows);
   
   return true;
}

bool DescriptorAnalysis::tovector(const Matrix3d& X, VectorXd& v) const
{
    if(X.rows()*X.cols() != v.rows())
    {
      std::cerr<<"bad matrix to vector conversion"<<std::endl; 
      return false;
    }
    size_t ncols = X.cols();
    size_t nrows = X.rows();
    
    for(size_t i = 0 ; i < ncols ; ++i)
        v.segment(i*nrows, nrows) = X.col(i);

    return true;
}

bool DescriptorAnalysis::tovector(const MatrixXd& X, VectorXd& v) const
{
    if(X.rows()*X.cols() != v.rows())
    {
      std::cerr<<"bad matrix to vector conversion"<<std::endl; 
      return false;
    }
    size_t ncols = X.cols();
    size_t nrows = X.rows();
    
    for(size_t i = 0 ; i < ncols ; ++i)
        v.segment(i*nrows, nrows) = X.col(i);

    return true;
}



bool DescriptorAnalysis::setGenericPattern(std::vector< Vector3d > &genpat)
{
  if(genpat.size() == 0)
    return false;
  
  size_t size = genpat.size();
  
  m_genpat.resize(3*size);
  
  for(size_t i = 0; i < size; ++i)
  {
    Vector3d &g = genpat[i];
     m_genpat(i) = g(0);
     m_genpat(i + size) = g(1);
     m_genpat(i + 2*size) = g(2);
  }
  return true;
}


bool DescriptorAnalysis::optimizeLPFbyICP(MatrixXd& V,
                                          MatrixXd& T,
                                          MatrixXd& tloc,
                                          MatrixXd& R,
                                          std::vector<Vector3d> &t1s,
                                          std::vector<Vector3d> &t2s,
                                          std::vector<Vector3d> &ns   
                                          )
{
  int size = m_X.rows()/3;
  double w = 1./((double)size);
  R.resize(m_X.cols(),9);
  tloc.resize(m_X.cols(),3);
  

#pragma omp parallel for
  for(int i = 0; i < m_X.cols(); ++i)
  {
      VectorXd v = V.col(i) + m_genpat;
      MatrixXd y(size, 3);
      y.col(0) =v.segment(0,size);
      y.col(1) =v.segment(size, size);
      y.col(2) =v.segment(2*size, size);

      MatrixXd x(size,3);
      v = m_X.col(i) + m_genpat;
      x.col(0) =v.segment(0,size);
      x.col(1) =v.segment(size, size);
      x.col(2) =v.segment(2*size, size);
      
      
      Vector3d mux = w * x.colwise().sum();
      Vector3d muy = w * y.colwise().sum();
      MatrixXd cx,cy;
      cx = x.rowwise() - mux.transpose();
      cy = y.rowwise() - muy.transpose();
      
      
      Matrix3d S = w * cx.transpose() * cy;
      
      JacobiSVD<Matrix3d> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
      
      Matrix3d temp = Matrix3d::Identity();
      Matrix3d U = svd.matrixU();
      Matrix3d V = svd.matrixV();
      temp(2,2) = (V*U.transpose()).determinant();
      Matrix3d Rloc = V * temp * U.transpose();
      
      Vector3d t = (muy - Rloc*mux);
      
      tloc.row(i) = t.transpose();
      
      //updating m_X
      x = x * Rloc.transpose();
      x.rowwise() += t.transpose();
      tovector(x, v);
      m_X.col(i) = v - m_genpat;
      
      //compute new parameterization;
      Vector3d t1n = Rloc(0,0) * t1s[i] + Rloc(0,1) * t2s[i] + Rloc(0,2) * ns[i];
      Vector3d t2n = Rloc(1,0) * t1s[i] + Rloc(1,1) * t2s[i] + Rloc(1,2) * ns[i];
      Vector3d nn  = Rloc(2,0) * t1s[i] + Rloc(2,1) * t2s[i] + Rloc(2,2) * ns[i];
      
      VectorXd r(9);
      tovector(Rloc,r);
      
      R.row(i) = r;
      
      normalize(t1n);
      normalize(t2n);
      normalize(nn);
      
      t1s[i]=t1n;
      t2s[i]=t2n;
      ns[i]=nn;
      
      T.row(i) -= t(0)*t1s[i] + t(1)*t2s[i] + t(2)*ns[i];
  }
  return true;
}
