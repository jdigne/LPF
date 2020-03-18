/**
 * This file is part of lpf
 *
 * @file Parameterizer.cpp
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
#include "Parameterizer.h"
#include "Eigen/Dense"
#include "utils.h"
#include "OctreeIterator.h"

using Eigen::Matrix3d;
using Eigen::Vector3d;
typedef TOctreeIterator<Vector3d> OctreeIterator;


Parameterizer::Parameterizer()
{

}



Parameterizer::Parameterizer(double radius, Octree *octree)
{
  m_octree = octree;
  m_radius = radius;
}

Parameterizer::~Parameterizer()
{
  m_octree = NULL;
}


unsigned int Parameterizer::compute_tangent_params(const std::vector< Vector3d >& points,
                                                   std::vector< Vector3d >& t1s,
                                                   std::vector< Vector3d >& t2s,
                                                   std::vector< Vector3d >& normals)
{
  t1s.clear();
  t2s.clear();
  normals.clear();
  unsigned int npoints = points.size();
  unsigned int nneighbor = 0;
  unsigned int max_neighbor = 0;
  t1s.assign(npoints, Vector3d(1,0,0));
  t2s.assign(npoints, Vector3d(0,1,0));
  normals.assign(npoints, Vector3d(0,0,1));
  
  unsigned int ndegenerate = 0;
  double sigma = -0.5/(m_radius*m_radius);
  

  #pragma omp parallel for
  for(int pi = 0; pi < npoints; ++pi)
  {
    OctreeIterator iterator(m_octree);
    iterator.setDepth(1);
    
    std::vector<Vector3d*> neighbors;
    std::vector<double> distances;
    iterator.getNeighbors(points[pi], neighbors, distances);
    
    assert(neighbors.size()>0);
    
    if(neighbors.size()<6)
    {
#pragma omp critical
      ndegenerate = ndegenerate + 1;
      continue;
    }
    //build covariance matrix
    double sum_w = 0;
    unsigned int npts = neighbors.size();
    
    Matrix3d a;
    Vector3d m;
    a.setZero();
    m.setZero();
    for(int i = 0; i < distances.size(); ++i)
    {
      double d = distances[i];
      double w = exp(sigma*d);
      a += w * (*neighbors[i]) * (*neighbors[i]).transpose();
      m += w * (*neighbors[i]);
      sum_w += w;
    }
    
    m = m / sum_w;
    a = a / sum_w - m * m.transpose();
    
    Eigen::SelfAdjointEigenSolver<Matrix3d> eigensolver(a);
    if (eigensolver.info() != Eigen::Success) abort();
    Matrix3d evec = eigensolver.eigenvectors();
    normals[pi] << evec(0,0), evec(1,0), evec(2,0);
    t2s[pi] << evec(0,1), evec(1,1), evec(2,1);
    t1s[pi] << evec(0,2), evec(1,2), evec(2,2);
    
    Vector3d test;
    test = t1s[pi].cross(t2s[pi]);
    normalize(test);
    
    if(test.dot(normals[pi]) < 0)
    {
      normals[pi] *= -1;
    }
    
#pragma omp atomic
    nneighbor = nneighbor + npts;
    
    if(npts>max_neighbor)
    {
      #pragma omp critical
      max_neighbor = npts;
      
    }
  }
  nneighbor = nneighbor / ((double)npoints);
 
  std::cout<<ndegenerate<<" points with degenerate neighborhoods"<<std::endl;
  std::cout<<"Mean number of neighbors "<< nneighbor <<"."<<std::endl;
  std::cout<<"Max number of neighbors "<< max_neighbor <<"."<<std::endl;
  
  return nneighbor;
}



unsigned int Parameterizer::compute_random_params(const std::vector< Vector3d >& points,
                                                   std::vector< Vector3d >& t1s,
                                                   std::vector< Vector3d >& t2s,
                                                   std::vector< Vector3d >& normals)
{
  t1s.clear();
  t2s.clear();
  normals.clear();
  unsigned int npoints = points.size();
  unsigned int nneighbor = 0;
  unsigned int max_neighbor = 0;
  t1s.assign(npoints, Vector3d(1,0,0));
  t2s.assign(npoints, Vector3d(0,1,0));
  normals.assign(npoints, Vector3d(0,0,1));
  
  unsigned int ndegenerate = 0;
  double sigma = -0.5/(m_radius*m_radius);
  

  #pragma omp parallel for
  for(int pi = 0; pi < npoints; ++pi)
  {
      Vector3d u = Vector3d::Random();
      normalize(u);
      
      Vector3d v = Vector3d::Random();
      v = v - v.dot(u)*u;
      normalize(v);
      
      Vector3d w = u.cross(v);
      normalize(w);
      
      t1s[pi] = u;
      t2s[pi] = v;
      normals[pi] = w;
  }
  return 0;
}

