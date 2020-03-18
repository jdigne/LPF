/**
 * This file is part of lpf
 *
 * @file IndexedPoint.h
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
#ifndef POINT_H
#define POINT_H

#include <Eigen/Dense>

using Eigen::Vector3d;

class IndexedPoint
{
public :
  IndexedPoint()
  {
    coordinates.setZero();
    valid = -1;
    index = -1;
  }
  
  
  IndexedPoint(const Vector3d &v, int _index, int _valid)
  {
    coordinates = v;
    valid = _valid;
    index = _index;
  }
  
  
  Vector3d coordinates;
  
  int index;
  
  int valid;
  
  
  double x(){return coordinates.x();}
  double y(){return coordinates.y();}
  double z(){return coordinates.z();}
};

inline  Vector3d operator- (const Vector3d &u, const IndexedPoint &v)
{
  return u - v.coordinates;
}

inline  Vector3d operator- (const IndexedPoint &v, const Vector3d &u)
{
  return v.coordinates - u;
}


#endif
