/**
 * This file is part of lpf
 *
 * @file DescriptorReconstruction.h
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
#ifndef DESCRIPTOR_RECONSTRUCTION_H
#define DESCRIPTOR_RECONSTRUCTION_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Octree.h"
#include "OctreeIterator.h"
#include "IndexedPoint.h"


#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::Vector3d;


class DescriptorReconstruction
{
public :

typedef TOctree<IndexedPoint> Octree;
typedef TOctreeIterator<IndexedPoint> OctreeIterator;
  
  DescriptorReconstruction(const double radius);
  
  DescriptorReconstruction(const double radius, Vector3d &origin, double size, unsigned int depth);
  
  ~DescriptorReconstruction();
  
  void buildOctree(const Vector3d &origin, double size, const double precision);

  template<class T>
  bool reconstructAndMerge(const std::vector< Vector3d >& seeds,
                           const double precision,
                           const std::vector< Vector3d >& t1s,
                           const std::vector< Vector3d >& t2s,
                           const std::vector< Vector3d >& normals,
                           const std::vector< Vector3d >& qpoints,
                           const std::vector< T* >& descs,
                           std::vector< Vector3d >& final_points)
  {
    int nseeds = seeds.size();
    int nptspseed = qpoints.size();
    
    std::vector< IndexedPoint > temppts;
    
  #pragma omp parallel for
    for(int ind = 0; ind < nseeds; ++ind)
    {
      std::vector<Vector3d> locpoints;
      descs[ind]->reconstruct(seeds[ind], t1s[ind], t2s[ind], normals[ind], m_radius, qpoints, locpoints);
  #pragma omp critical
      {
        for(int i = 0; i < locpoints.size() ; ++i)
        {
          temppts.push_back(IndexedPoint(locpoints[i], i, 1));
        }
      }
    }
    
    
    m_octree->addInitialPoints(temppts.begin(), temppts.end());
  
    double sqr = precision * precision;
    
    std::cout<<"Merging precision:"<<precision<<std::endl;
    int N = temppts.size();
    
    int nresampled = 0;
    for(int ind = 0; ind < temppts.size(); ++ind)
    {
       IndexedPoint &p = temppts[ind];
       OctreeIterator iterator(m_octree);
       iterator.setR(precision);
       std::vector<IndexedPoint*> neighbors;
       std::vector<double> sqdists;
       iterator.getNeighbors(p.coordinates, neighbors, sqdists);
       
       int n = neighbors.size();
       Vector3d merged;
       merged.setZero();
       int t = 0;
       for(int j = 0 ; j < n ; ++j)
       {
         merged += neighbors[j]->coordinates;
         neighbors[j]->valid = -1;
         t ++;
       }
       
       if(t > 0)
       {
         merged = merged*1.0/((double)t);
         final_points.push_back(merged);
         nresampled++;
         
         if(nresampled % 100000 == 0)
           std::cout << nresampled << " points resampled (" << floor(100*(double)ind/(double)N) << "%)." << std::endl;
       }
    }
    return true;
  }  
  

  
//WARNING: Will produce huge point sets with duplicated points
  template<class T>
  bool reconstructAndSave(const std::vector<Vector3d> &seeds,
                   const std::vector<Vector3d> &t1s,
                   const std::vector<Vector3d> &t2s,
                   const std::vector<Vector3d> &normals,
                   const std::vector< Vector3d >& qpoints,
                   const std::vector<T*> &descs,
                   const char *filename)
  {
    int nseeds = seeds.size();
    std::cout<<"Saving only points given by u+v"<<std::endl;
    
    std::ofstream out;
    out.open(filename);
    if(! out)
    {
      std::cerr<<"cannot open file to save reconstructed points. Exiting."<<std::endl;
      return false;
    }
    
    for(int ind = 0; ind < nseeds; ++ind)
    {
      std::vector<Vector3d> locpoints;
  
      descs[ind]->reconstruct(seeds[ind], t1s[ind], t2s[ind], normals[ind], m_radius, qpoints, locpoints);

      while(! locpoints.empty())
      {
        out << locpoints.back()[0]<<"\t" <<locpoints.back()[1]<<"\t" <<locpoints.back()[2]<< "\n";
        locpoints.pop_back();
      }

      if( (ind+1) % 50000 == 0 )
        std::cout << ind+1 << " processed points." << std::endl; 
    }
    out.close();
    return true;
  }
  
  
//WARNING: Will produce huge point sets with duplicated points
  template<class T>
  bool reconstructAndSaveTargetZones(const std::vector<Vector3d> &seeds,
                   const std::vector<Vector3d> &t1s,
                   const std::vector<Vector3d> &t2s,
                   const std::vector<Vector3d> &normals,
                   const std::vector<T*> &descs,
                   const char *filename)
  {
    int nseeds = seeds.size();
    std::cout<<"Saving whole target zones"<<std::endl;
    
    std::ofstream out;
    out.open(filename);
    if(! out)
    {
      std::cerr<<"cannot open file to save reconstructed points. Exiting."<<std::endl;
      return false;
    }
    
    for(int ind = 0; ind < nseeds; ++ind)
    {
      std::vector<Vector3d> locpoints;
  
      descs[ind]->save_target_zone(seeds[ind], t1s[ind], t2s[ind], normals[ind], m_radius, locpoints);

      while(! locpoints.empty())
      {
        out << locpoints.back()[0]<<"\t" <<locpoints.back()[1]<<"\t" <<locpoints.back()[2]<< "\n";
        locpoints.pop_back();
      }

      if( (ind+1) % 50000 == 0 )
        std::cout << ind+1 << " processed points." << std::endl; 
    }
    out.close();
    return true;
  }
  
  template<class T>
  bool buildDescriptorsFromMatrix(const MatrixXd &X,
                                  std::vector<T*> &descs)
  {
    unsigned int n = X.rows()/3;
    descs.resize(X.cols());
    
    for(size_t i = 0; i < X.cols(); ++i)
      descs[i] = new T(n);
    
#pragma omp parallel for
    for(size_t i = 0; i < X.cols(); ++i)
    {
      descs[i]->readMatrix(X.col(i));
    }
    return true;
  }
  
  
  void adaptSeeds(std::vector<Vector3d> &seeds,
                  MatrixXd &t);
  
private :
  
  double m_radius;
  
  Octree * m_octree;
  
};

#endif
