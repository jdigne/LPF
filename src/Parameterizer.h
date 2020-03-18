/**
 * This file is part of lpf
 *
 * @file Parameterizer.h
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
#ifndef PARAMETERIZER_H
#define PARAMETERIZER_H

#include <Eigen/Dense>
#include "Octree.h"

typedef TOctree<Vector3d> Octree;
using Eigen::Vector3d;

class Parameterizer
{
  public :
  
  Parameterizer();
  
  Parameterizer(double radius, Octree *octree);
  
  ~Parameterizer();
  
  
  unsigned int compute_tangent_params(const std::vector<Vector3d> &qpoints,
                                      std::vector<Vector3d> &t1,
                                      std::vector<Vector3d> &t2,
                                      std::vector<Vector3d> &normals);
  
  unsigned int compute_random_params(const std::vector<Vector3d> &qpoints,
                                      std::vector<Vector3d> &t1,
                                      std::vector<Vector3d> &t2,
                                      std::vector<Vector3d> &normals);
  
  private :
  
  double m_radius;
  
  Octree *m_octree;
  
};

#endif
