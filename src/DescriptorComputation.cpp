/**
 * This file is part of lpf
 *
 * @file DescriptorComputation.cpp
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
#include "DescriptorComputation.h"
#include "utils.h"

using Eigen::MatrixXd;

DescriptorComputation::DescriptorComputation()
{

}

DescriptorComputation::DescriptorComputation(double radius, std::vector<Vector3d>& queries)
{
  m_radius = radius;
  m_wide_radius = 1.2 * m_radius;
  m_maxneighbor = 200;
  m_queries.clear();
  m_genpat.resize(queries.size(),3);
  for(size_t i = 0; i < queries.size(); ++i)
  {
    Vector3d &v= queries[i];
    m_queries.push_back(v);
    m_genpat.row(i) = v;
  }
}

DescriptorComputation::DescriptorComputation(double radius, std::vector< Vector3d >& queries, Octree* octree)
{
  m_octree = octree;
  m_radius = radius;
  m_wide_radius = 1.0 * m_radius;
  m_queries.clear();
  m_genpat.resize(queries.size(),3);
  for(size_t i = 0; i < queries.size(); ++i)
  {
    Vector3d &v= queries[i];
    m_queries.push_back(v);
    m_genpat.row(i) = v;
  }
}


DescriptorComputation::~DescriptorComputation()
{
  m_radius = 0;
  m_wide_radius =  0;
  m_queries.clear();
}

bool DescriptorComputation::express_locally(const Vector3d &seed,
                                            const Vector3d& t1,
                                            const Vector3d& t2,
                                            const Vector3d& n,
                                            const std::vector<Vector3d*> &neighbors,
                                            std::vector< Vector3d >& locpoints)
{
  std::vector<Vector3d*>::const_iterator vi;
  for(vi = neighbors.begin(); vi != neighbors.end() ; ++vi)
  {
    Vector3d v = **vi;
    v -= seed;
    double xloc = v.dot(t1)/m_radius;
    double yloc = v.dot(t2)/m_radius;
    double zloc = v.dot(n)/m_radius;
    locpoints.push_back(Vector3d(xloc, yloc, zloc));
  }
  return true;
}
