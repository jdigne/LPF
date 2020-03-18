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
#ifndef DENOISING_H
#define DENOISING_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "LPF.h"
#include "Octree.h"
#include "OctreeIterator.h"
#include "IndexedPoint.h"


#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::Vector3d;


class Denoising
{
  typedef TOctree<IndexedPoint> Octree;
  typedef TOctreeIterator<IndexedPoint> OctreeIterator;
  
public :

  Denoising(Vector3d& origin, double size, unsigned int depth);
  
  ~Denoising();

  
  bool denoise(std::vector< Vector3d >& points,
               const double radius,
               const std::vector< Vector3d >& seeds,
               const std::vector< Vector3d >& t1s,
               const std::vector< Vector3d >& t2s,
               const std::vector< Vector3d >& normals,
               const std::vector< Vector3d >& qpoints,
               const std::vector< LPF* >& descs,
               const double lambda);
  
private :

  Octree *m_octree;
  
  Vector3d m_origin;
  
  double size;
};

#endif
