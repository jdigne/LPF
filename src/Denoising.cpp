/**
 * This file is part of lpf
 *
 * @file Denoising.cpp
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
#include "Denoising.h"


Denoising::Denoising(Vector3d& origin, double size, unsigned int depth)
{
  origin = origin.array() - 0.1*size;
  m_octree = new Octree(origin, 1.2*size, depth);
}


Denoising::~Denoising()
{
  delete m_octree;
}




bool Denoising::denoise(std::vector< Vector3d >& points,
                        const double radius,
                        const std::vector< Vector3d >& seeds,
                        const std::vector< Vector3d >& t1s,
                        const std::vector< Vector3d >& t2s,
                        const std::vector< Vector3d >& normals,
                        const std::vector< Vector3d >& qpoints,
                        const std::vector< LPF* >& descs, 
                        const double lambda)
{
  for(int i = 0; i < seeds.size() ; ++i)
  {
    IndexedPoint p(seeds[i], i, 1);
    m_octree->addInitialPoint(p);
  }
  
  Vector3d origin = m_octree->getOrigin();
  double min = origin.x() < origin.y() ? origin.x() : origin.y();
  min = min < origin.z() ? min : origin.z();
  min = min - 2 * m_octree->getSize();
  
  double sqrad = radius*radius;
  double normalizer = 1.0/(1.0+lambda);
#pragma omp parallel for
  for(int ind = 0; ind < points.size(); ++ind)
  {
     std::vector<IndexedPoint*> neighbors;
     OctreeIterator iterator(m_octree);
     iterator.setR(radius);
     
     iterator.getNeighbors(points[ind], neighbors);
    
     Vector3d denoised;
     denoised.setZero();
     double sum_w = 0.0;
     
     if(neighbors.size()==0)
        continue;
     
     for(int i = 0; i < neighbors.size() ; ++i)
     {
       int index = neighbors[i]->index;
       
       LPF *lpf = descs[index];
       Vector3d point = points[ind];
       
       double w = lpf->vote(point, seeds[index],
                            t1s[index], t2s[index], normals[index],
                            qpoints, radius);
       if(w <1e-16)
         continue;
       
       denoised += w * point; 
       sum_w += w;
     }
     
     if(sum_w < 1e-16)
     {
       points[ind] = Vector3d(min,min,min);
       continue;
     }
     
     double t = 1./sum_w;
    points[ind] = t * denoised;
  }
  std::cout<<std::endl;
  min = min + m_octree->getSize();
  
  std::vector<Vector3d>::iterator it = points.begin();
  while(it != points.end())
  {
      if(it->x() < min )
        it = points.erase(it);
      else
        ++it;
  }
  return true;
}
