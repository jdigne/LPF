/**
 * This file is part of lpf
 *
 * @file utils.h
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
#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <Eigen/Dense>


using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::Matrix2d;
using Eigen::Vector2d;
using Eigen::VectorXd;
using Eigen::JacobiSVD;


static bool read_qpoints(const char* filename, std::vector<Vector3d> &queries)
{
  std::ifstream f;
  f.open(filename);
  if(!f)
  {
    std::cerr<<"cannot read query points from file "<<filename<<std::endl;
    return false;
  }
  
  double x,y,z;
  while(f>>x>>y>>z)
  {
    queries.push_back(Vector3d(x,y,z));
  }
  f.close();
  std::cout<<queries.size()<<" query points in the irregular grid."<<std::endl;
  return true;
}

static bool generate_1d_grid(unsigned int nbins, std::vector<Vector3d> &queries)
{
  //building the grid
  double s = 2./((double)nbins);
  double hs = 1.0;
  queries.clear();

  for(size_t i = 0; i < nbins; ++i)
    queries.push_back(Vector3d(i * s - hs,0,0));
  return true;
}

static bool generate_2d_grid(unsigned int nbins, std::vector<Vector3d> &queries)
{
  //building the grid
  double s = 2./((double)nbins);
  double hs= 1.0;
  queries.clear();

  for(size_t i = 0; i < nbins; ++i)
    for(size_t j = 0; j < nbins; ++j)
    {
      double x = i * s - hs; 
      double y = j * s - hs; 
      if(x*x + y*y < 1.0)
        queries.push_back(Vector3d(x, y, 0));
    }
  std::cout<<"Generic Pattern generated with "<<queries.size()<<" points"<<std::endl;
  return true;
}

static bool generate_3d_grid(unsigned int nbins, std::vector<Vector3d> &queries)
{
  //building the grid
  double s = 2./((double)nbins);
  double hs= 1.0;
  queries.clear();

  for(size_t i = 0; i < nbins; ++i)
    for(size_t j = 0; j < nbins; ++j)
      for(size_t k = 0; k < nbins; ++k)
      {
        double x = i * s - hs; 
        double y = j * s - hs; 
        double z = k * s - hs; 
        if(x*x + y*y + z*z < 1.0)
         queries.push_back(Vector3d(x, y, z));
      }
  return true;
}

static bool random_points(unsigned int nbins, std::vector<Vector3d> &queries)
{
  //building the grid
  srand(std::time(NULL));
  queries.clear();

  for(size_t i = 0; i < nbins; ++i)
  {
    double x = 2.0 * rand()/((double) RAND_MAX) - 1.0;
    double y = 2.0 * rand()/((double) RAND_MAX) - 1.0;
    double z = 2.0 * rand()/((double) RAND_MAX) - 1.0;
    queries.push_back(Vector3d(x,y,z));
  }
  return true;
}


static bool saveQueries(std::vector<Vector3d> &queries, const char* filename)
{
  std::ofstream f;
  f.open(filename);
  
  if(!f)
  {
    return false;
  }
  
  for(size_t i = 0; i < queries.size(); ++i)
  {
    f<<queries[i]<<std::endl;
  }
  
  f.close();
  return true;
}






static double project(Vector3d &t, double x, double y, double z)
{
  return t.transpose()*(Vector3d() << x, y, z).finished();
}

static void normalize(Vector3d &t)
{
   double n = t.norm();
   t = t/n;
}

inline double sqr(double t)
{
    return t*t;
}


inline void check_frames(std::vector<Vector3d> &t1s, std::vector<Vector3d> &t2s, std::vector<Vector3d> &ns)
{
#pragma omp parallel for
  for(int i = 0; i < t1s.size(); ++i)
  {
    assert(t1s[i].dot(t2s[i])<1e-10);
    assert(t1s[i].dot(ns[i])<1e-10);
    assert(t2s[i].dot(ns[i])<1e-10);
    assert(abs(t1s[i].norm()-1.0)<1e-10);
    assert(abs(t2s[i].norm()-1.0)<1e-10);
    assert(abs(ns[i].norm()-1.0)<1e-10);
    assert(abs((t1s[i].cross(t2s[i])).dot(ns[i]) - 1.0) <1e-10);
    assert(abs((t2s[i].cross(ns[i])).dot(t1s[i]) - 1.0) <1e-10);
    assert(abs((ns[i].cross(t1s[i])).dot(t2s[i]) - 1.0) <1e-10);
  }
}


inline void svd_transform_estimate_3(const MatrixXd &p, const MatrixXd &q, Matrix3d &R, Vector3d &t)
{
  assert(p.rows() == q.rows() && p.cols() == q.cols());
  assert(p.cols() == 3);
  
  Vector3d barp = p.colwise().mean();
  Vector3d barq = q.colwise().mean();
  
  Matrix3d M  = (p.rowwise() - barp.transpose()).transpose()*(q.rowwise() - barq.transpose());
  
  JacobiSVD<Matrix3d> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  
  //R= V *diag(1,1,det(V*U^T)))*U^T
  MatrixXd temp = Matrix3d::Identity();
  MatrixXd U = svd.matrixU();
  MatrixXd V = svd.matrixV();
  temp(2,2) = (V*U.transpose()).determinant();
  R = V * temp * U.transpose();
  
  t = (barq - R*barp);
}


inline void svd_transform_estimate_2(const MatrixXd &p, const MatrixXd &q, Matrix2d &R, Vector2d &t)
{
  assert(p.rows() == q.rows() && p.cols() == q.cols());
  assert(p.cols() == 2);
  
  Vector2d barp = p.colwise().mean();
  Vector2d barq = q.colwise().mean();
  
  
  Matrix2d M  = (p.rowwise() - barp.transpose()).transpose()*(q.rowwise() - barq.transpose());
  JacobiSVD<Matrix2d> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  
  //R= V *diag(1,1,det(V*U^T)))*U^T
  Matrix2d temp = Matrix2d::Identity();
  Matrix2d U = svd.matrixU();
  Matrix2d V = svd.matrixV();
  temp(1,1) = (V*U.transpose()).determinant();
  R = V * temp * U.transpose();
  
  t = (barq - R*barp);
}

#endif
