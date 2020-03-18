/**
 * This file is part of lpf
 *
 * @file Optimization.cpp
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
#include "Optimization.h"
#include<iostream>
#include<fstream>
#include <algorithm>

using Eigen::JacobiSVD;

Optimization::Optimization()
{

}


Optimization::~Optimization()
{

}

bool Optimization::initializeRandomDictionary(const MatrixXd& X, MatrixXd& D)
{
  int *indices = new int[X.cols()];
  rand_permutation(X.cols(), indices);

  for(size_t i = 0; i < D.cols(); ++i)
  {
    double t = X.col(indices[i]).norm();
    D.col(i) = X.col(indices[i])/t;
  }

  delete indices;
  return true;
}


//OMP-Cholesky
void Optimization::omp(const MatrixXd& X, const MatrixXd& D, MatrixXd& coeffs, int sparsity)
{
  coeffs.resize(D.cols(),X.cols());
  coeffs.setZero();
  
#pragma omp parallel for
  for(int i = 0; i < X.cols(); ++i)
  {
    VectorXd x = X.col(i);
    VectorXd r = x;
    VectorXd c(coeffs.rows());
    VectorXd alpha = D.transpose()*x;
    vector<int> indices;
    indices.reserve(sparsity);
    MatrixXd L(1,1);
    MatrixXd Dt;
    L(0,0) = 1.0;

    for(int j = 0; j< sparsity ; ++j)
    {
      int index = 0;
      double max = -1;
      for(int k = 0; k < D.cols(); ++k)
      {
        if(std::find(indices.begin(), indices.end(), k) == indices.end())
        {
           double t = (D.col(k).transpose() * r).array().abs()(0);
           if(t > max)
           {
             max = t;
             index = k;
           }
        }
      }
      
      if(j > 0)
      {
        VectorXd w = Dt.transpose()*D.col(index);
        L.triangularView<Eigen::Lower>().solveInPlace(w);
        
        MatrixXd Ltemp = L;
        int nrows = L.rows();
        int ncols = L.cols();
        L.resize(nrows+1, ncols+1);
        L.setZero();
        L.block(0,0, nrows, ncols) = Ltemp;
        L.block(nrows, 0 , 1, ncols) = w.transpose();
        double lastcoeff =  1.0 - w.transpose()*w;
       
        assert( lastcoeff > -1e-16);
        lastcoeff = 0 > lastcoeff ? 0 : lastcoeff;
        
        lastcoeff = std::sqrt(lastcoeff);
        L(nrows, ncols) = lastcoeff;
      }
      indices.push_back(index);
      VectorXd alphaI;
      subvector(alpha, indices, alphaI);
      
      //solve Ly = alpha_I
      VectorXd y = L.triangularView<Eigen::Lower>().solve(alphaI);
      //solve L^Tx = y
      c = L.transpose().triangularView<Eigen::Upper>().solve(y);

      submatrix(D, indices, Dt);
      r = x - Dt*c;
      if(r.squaredNorm()<1e-12)
      {
        break;//residual already minimized
      }
    }
    
    for(int k = 0; k < indices.size();++k)
    {
      coeffs(indices[k],i) = c(k);
    }
  }
}




void Optimization::homotopy(const MatrixXd& X, const MatrixXd& D, MatrixXd& coeffs, double lambdaf)
{
  
  assert(checkNonZero(X));
  assert(checkNonZero(D));
  
  coeffs.resize(D.cols(),X.cols());
  coeffs.setZero();
  MatrixXd DtD = D.transpose()*D;

  #pragma omp parallel for
  for(int i = 0; i < X.cols(); ++i)
  {
    VectorXd x = X.col(i);
    VectorXd dir(coeffs.rows());
    VectorXd alpha(coeffs.rows());
    alpha.setZero();
    VectorXd cI, sgncI, v, dirt;
    MatrixXd DI,M;
    VectorXd c =  D.transpose()*x;
    
    vector<int> indices;
    
    int index = 0;
    VectorXd t = c.array().abs();
    double lambda = t.maxCoeff(&index);
    indices.push_back(index);
      
    M.resize(1,1);
    M(0,0) = 1.0 / DtD(index, index);
    
    while(lambda > lambdaf)
    {
      subvector(c, indices, cI);
      sign(cI, sgncI);

      dirt = M * sgncI;

      dir.setZero();
      for(int k = 0; k < indices.size();++k)
      {
        dir(indices[k]) = dirt(k);
      }

      submatrix(D, indices, DI);
      v = DI*dirt;

      double minp, minm;
      int indp, indm;

      minm = minp = lambda+1;
      for(int k = 0; k < alpha.rows();++k)
      {
          if(std::find(indices.begin(), indices.end(), k) != indices.end())
          {
            double v3 = - alpha(k)/dir(k);
            if(v3 > 0 && v3 < minm)
            {
              minm = v3;
              indm = k;
            }
          }
          else
          {
            double v1 = (lambda - c(k))/(1 - D.col(k).transpose()*v);
            double v2 = (lambda + c(k))/(1 + D.col(k).transpose()*v);
            if(v1 > 0 && v1 < minp)
            {
              minp = v1;
              indp = k;
            }
            if(v2 > 0 && v2 < minp)
            {
              minp = v2;
              indp = k;
            }
          }
      }
      
      if((minm > lambda+0.5) && (minp > lambda+0.5))
        break;//min value is reached
      
      
      double min = 0;
      bool remove = false;
      if(minm < minp )
      {
         vector<int>::iterator it = std::find(indices.begin(), indices.end(), indm);
         int index_in_vector = it - indices.begin();
         if(index_in_vector >= M.cols())
           std::cout<<"Index in vector"<<index_in_vector<<std::endl;
         updateInverseIndexRemoved(M, index_in_vector);
         min = minm;
         indices[index_in_vector] = indices.back();
         indices.pop_back();
         remove = true;
      }
      else
      {
        updateInverseIndexAdded(M, DI, D.col(indp));
        indices.push_back(indp);
        min = minp;
      }
      
      if (lambda - min < lambdaf)
      {
        alpha = alpha + (lambda - lambdaf)*dir;
        lambda = lambdaf;
      }
      else
      {
        alpha = alpha + min * dir;
        lambda = lambda - min;
      }
      
      if(remove)
        alpha(indm)=0;
      
      c = D.transpose()*x - DtD*alpha;
      
      if(c.squaredNorm()<1e-12)//min error is reached
        break;
    }
    coeffs.col(i) = alpha;
  }
}


//dictionary update

//D should be initialized
void Optimization::mod_l1(const MatrixXd& X, MatrixXd& D, MatrixXd& coeffs, double lambda, int niter)
{
  MatrixXd B, C;
  
  for(int n = 0; n < niter ; ++n)
  {
    //compute sparse codes
    homotopy(X, D, coeffs, lambda);
    
    assert(!checkNaN(X));
    assert(!checkNaN(D));
    assert(!checkNaN(coeffs));
    B = X * coeffs.transpose();
    C = coeffs * coeffs.transpose();  
    //update dictionary D //not exactly mod_l1 as described by Lee 2007, but should fit
    
#pragma omp parallel for
    for(size_t j = 0; j < (unsigned int)D.cols(); ++j)
    {
      VectorXd dj(D.rows());
      if(abs(C(j,j))<1e-16)
        continue;

      dj = 1./C(j,j) * (B.col(j) - D * C.col(j)) + D.col(j);
      double norm = dj.norm();
      norm = norm > 1.0 ? norm : 1.0;
      dj = 1.0/norm *dj;
      D.col(j)=dj;
    }
    
  }
  homotopy(X, D, coeffs, lambda);
}





//private methods

void Optimization::rand_permutation(unsigned int n, int* numbers)
{
  for(size_t i = 0; i< n; ++i)
  {
    int r = rand() % (i+1);
    numbers[i] = numbers[r];
    numbers[r] = i;
  }
}



void Optimization::submatrix(const MatrixXd& M, vector< int >& indices, MatrixXd &subM)
{
  subM.resize(M.rows(), indices.size());
  for(int i = 0; i < indices.size(); ++i)
  {
    subM.col(i) = M.col(indices[i]);
  }
}




void Optimization::subvector(const VectorXd& V, vector< int >& indices, VectorXd &subV)
{
  subV.resize(indices.size());
  for(int i = 0; i < indices.size(); ++i)
  {
    subV(i) = V(indices[i]);
  }
}



void Optimization::sign(const VectorXd& X, VectorXd& sgnX)
{
  sgnX.resize(X.rows());
  for(int i = 0; i < X.rows(); ++i)
  {
    if(X(i) < 0)
      sgnX(i) = -1;
    else if(X(i) > 0)
      sgnX(i) = 1;
    else
      sgnX(i) = 0;
  }
}

bool Optimization::checkNaN(const MatrixXd &M)
{
  for(size_t i = 0; i < M.rows(); ++i)
    for(size_t j = 0; j < M.cols(); ++j)
      if(std::isnan(M(i,j)))
        return true;
  return false;    
}

bool Optimization::checkNaN(const VectorXd &v)
{
  for(size_t i = 0; i < v.rows(); ++i)
      if(std::isnan(v(i)))
        return true;
  return false;    
}

bool Optimization::checkNonZero(const VectorXd& v)
{
    if(v.squaredNorm() < 1e-16)
      return false;
    else return true;
}


bool Optimization::checkNonZero(const MatrixXd &M)
{
  for(int i = 0; i < M.cols() ; ++i)
  {
    if(! checkNonZero((VectorXd) M.col(i)))
      return false;
  }
  return true;
}




void Optimization::updateInverseIndexAdded(MatrixXd& M, const MatrixXd &DI, const VectorXd &d)
{
   int m = M.rows();
   MatrixXd Mnew(m + 1, m + 1);
   
   VectorXd b = DI.transpose() * d;
   VectorXd c = M*b;
   
   double s = (d.squaredNorm() - b.transpose()*c);
   s = 1.0/s;
   
   Mnew.block(0, 0, m, m) = M + s * c * c.transpose();
   Mnew.block(0, m, m, 1) = -s * c;
   Mnew.block(m, 0, 1, m) = -s * c.transpose();
   Mnew(m,m) = s;
   
   M = Mnew;
}


void Optimization::updateInverseIndexRemoved(MatrixXd& M, const int index)
{
   int mn = M.rows()-1;
   MatrixXd Mnew(mn, mn);
   
   permuteRowAndColLast(M, index);
   
   
   Mnew = M.block(0, 0, mn, mn);
   double s = M(mn, mn);
   VectorXd u3 = - M.block(0, mn, mn, 1);
   VectorXd u2 = u3/s;
   Mnew -= s*u2*u2.transpose();
   
   M = Mnew;
}

//The order of the indices will change because of this swap, to be handled if this function is used.
void Optimization::permuteRowAndColLast(MatrixXd& M, int i)
{
  int m = M.rows()-1;

  assert(i <= m);
  assert(i >= 0);
  M.row(i).swap(M.row(m));
  M.col(i).swap(M.col(m));
}

