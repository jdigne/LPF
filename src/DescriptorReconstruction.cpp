/**
 * This file is part of lpf
 *
 * @file DescriptorReconstruction.cpp
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
#include "DescriptorReconstruction.h"


using Eigen::Matrix3d;
using Eigen::Vector3d;



DescriptorReconstruction::DescriptorReconstruction(const double radius)
{
  m_radius = radius;
  m_octree = NULL;
}

DescriptorReconstruction::DescriptorReconstruction(const double radius, Vector3d& origin, double size, unsigned int depth)
{
  m_radius = radius;
  Vector3d morigin = origin.array() - 0.05*size;
  m_octree = new Octree(morigin, 1.1*size, depth);
}

void DescriptorReconstruction::buildOctree(const Vector3d& origin, double size, const double precision)
{
  Vector3d morigin = origin.array() - 0.05*size;
  double msize = 1.1 * size;
  int mdepth = floor( log2( msize / (2.0*precision) ));
  m_octree = new Octree(morigin, msize, mdepth);
}



DescriptorReconstruction::~DescriptorReconstruction()
{
  if(m_octree != NULL)
    delete m_octree;
}




void DescriptorReconstruction::adaptSeeds(std::vector<Vector3d> &seeds,
                                          MatrixXd &t)
{
  assert(t.rows() == seeds.size());
  
  t = t * m_radius;
  
#pragma omp parallel for
  for(int i = 0; i < seeds.size(); ++i)
  {
    Vector3d T = t.row(i);
    seeds[i] += T;
  }
}


