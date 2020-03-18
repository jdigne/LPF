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
#ifndef DESCRIPTOR_COMPUTATION_H
#define DESCRIPTOR_COMPUTATION_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include "Octree.h"
#include "OctreeIterator.h"

using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::MatrixXd;



class DescriptorComputation
{
  typedef TOctree<Vector3d> Octree;
  typedef TOctreeIterator<Vector3d> OctreeIterator;
  
public :
  
  DescriptorComputation();
  
  ~DescriptorComputation();
 
  //queries is the generic pattern
  DescriptorComputation(double radius, std::vector<Vector3d> &queries);
  
  //queries is the generic pattern
  DescriptorComputation(double radius, std::vector<Vector3d> &queries, Octree *octree);
  
  //to avoid memory consumption: saves directly in a binary file
  template<class T>
  bool computeAndSaveDescriptors(std::vector<Vector3d> &seeds,
                        std::vector<Vector3d> &points,
                        std::vector<Vector3d> &t1s,
                        std::vector<Vector3d> &t2s,
                        std::vector<Vector3d> &normals,
                        const char *filename)
  {
     double sqrad = m_radius*m_radius;
     std::ofstream f;
     f.open(filename, std::ios::out | std::ios::binary);
     if(!f)
     {
       std::cout<<"Could not open file "<<filename<<" for saving the descriptors. Exiting."<<std::endl;
       return false;
     }
     typename std::vector<T*>::const_iterator it;
     unsigned int size_desc = m_queries.size();
     unsigned int ndesc = seeds.size();
     
     f.write(reinterpret_cast<const char*>(&size_desc), sizeof(size_desc));
     f.write(reinterpret_cast<const char*>(&ndesc), sizeof(ndesc));
     
     
    for(size_t ind = 0; ind < ndesc; ++ind)
    {
      std::vector<Vector3d*> neighbors;
      std::vector<Vector3d> locpoints;
      
      OctreeIterator iterator(m_octree);
      //iterator.setDepth(0);
      iterator.setR(m_wide_radius);
     
      iterator.getNeighbors(seeds[ind], neighbors);
  
      express_locally(seeds[ind], t1s[ind], t2s[ind], normals[ind], neighbors, locpoints);
  
      T locdescriptor(m_queries, locpoints);
      locdescriptor.writeBin(f);
      
      if( (ind+1) % 100000 == 0 )
        std::cout << ind+1 << " processed points." << std::endl; 
    }
  
    f.close();
    return true;
  }
  
  
  
  template<class T>
  bool computeDescriptors(std::vector<Vector3d> &seeds,
                        std::vector<Vector3d> &points,
                        std::vector<Vector3d> &t1s,
                        std::vector<Vector3d> &t2s,
                        std::vector<Vector3d> &normals,
                        std::vector<T*> &descriptors)
  {
    double sqrad = 2.25 * m_radius*m_radius;
    int npoints = seeds.size();
    
    descriptors.assign(npoints, NULL);
  
    #pragma omp parallel for
    for(int ind = 0; ind < npoints; ++ind)
    {
      std::vector<Vector3d*> neighbors;
      std::vector<Vector3d> locpoints;
      
      OctreeIterator iterator(m_octree);
      //iterator.setDepth(0);
      iterator.setR(m_wide_radius);
      iterator.getNeighbors(seeds[ind], neighbors);
  
      express_locally(seeds[ind], t1s[ind], t2s[ind], normals[ind], neighbors, locpoints);
  
      descriptors[ind] = new T(m_queries, locpoints);
      
      if( (ind+1) % 100000 == 0 )
        std::cout << ind+1 << " processed points." << std::endl; 
    }
  
    return true;
  }
  
  
  
  
 
  template<class T>
  bool computeDescriptorsWithEfficiencyOptim(std::vector<Vector3d> &seeds,
                        std::vector<Vector3d> &points,
                        std::vector<Vector3d> &t1s,
                        std::vector<Vector3d> &t2s,
                        std::vector<Vector3d> &normals,
                        std::vector<T*> &descriptors)
{
  unsigned int nseeds = seeds.size();
  unsigned int nneighbor = 0;
  unsigned int max_neighbor = 0;
  descriptors.assign(nseeds, NULL);
  
  unsigned int ndegenerate = 0;
  double sigma = -0.5/(m_radius*m_radius);
  

#pragma omp parallel for
  for(int pi = 0; pi < nseeds; ++pi)
  {
    OctreeIterator iterator(m_octree);
    //iterator.setDepth(0);
    iterator.setR(m_wide_radius);
    
    std::vector<Vector3d*> neighbors;
    std::vector<double> distances;
    iterator.getNeighbors(seeds[pi], neighbors, distances);
    
    assert(neighbors.size()>0);
    
    std::vector<Vector3d> locpoints;
    express_locally(seeds[pi], t1s[pi], t2s[pi], normals[pi], neighbors, locpoints);
    descriptors[pi] = new T(m_queries, locpoints);
    
    descriptors[pi]->optimize_efficiency(locpoints, m_genpat, seeds[pi], t1s[pi], t2s[pi], normals[pi], m_radius);
    
    if( (pi+1) % 10000 == 0 )
    {
      #pragma omp critical
      std::cout << pi+1 << " processed points." << std::endl; 
    }
  }
  return true;
}


private :
  double m_radius;
  
  double m_wide_radius;
  
  unsigned int m_maxneighbor;
  
  Octree *m_octree;
  
  std::vector<Vector3d> m_queries;
  
  MatrixXd m_genpat;
  
  
  bool express_locally(const Vector3d &seed,
                       const Vector3d &t1,
                       const Vector3d &t2,
                       const Vector3d &n,
                       const std::vector<Vector3d*> &neighbors,
                       std::vector<Vector3d> &locpoints);
};

#endif
