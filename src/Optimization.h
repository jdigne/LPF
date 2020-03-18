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
#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <Eigen/Dense>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

class Optimization
{
public :
  Optimization();
  
  ~Optimization();
  
  bool initializeRandomDictionary(const MatrixXd& X, MatrixXd& D);
  

  void omp(const MatrixXd &X, const MatrixXd &D, MatrixXd &alpha, int sparsity);
  
  
  void homotopy(const MatrixXd &X, const MatrixXd &D, MatrixXd &alpha, double lambda);
  
  //dictionary learning
  void mod_l1(const MatrixXd& X, MatrixXd& D, MatrixXd& coeffs, double lambda, int niter);
  
  
private :
  
  void rand_permutation(unsigned int n, int* numbers);

  void submatrix(const MatrixXd &M, vector<int> &indices, MatrixXd &subM);
  
  void subvector(const VectorXd &V, vector<int> &indices, VectorXd &subV);
  
  void sign(const VectorXd &X, VectorXd &sgnX);
  
  bool checkNaN(const MatrixXd &M);
  
  bool checkNaN(const VectorXd &v);
  
  bool checkNonZero(const MatrixXd &M);
  
  bool checkNonZero(const VectorXd &v);
  
  void updateInverseIndexRemoved(MatrixXd& M, const int index);
  
  void updateInverseIndexAdded(MatrixXd& M, const MatrixXd &DI, const VectorXd& d);
  
  void permuteRowAndColLast(MatrixXd &M, int i);//permutes row i of M with last row of M and column i of M with last column of M.
  
};

#endif
