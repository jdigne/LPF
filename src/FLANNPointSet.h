/**
 * This file is part of lpf
 *
 * @file FLANNPointSet.h
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
#ifndef FLANN_POINT_SET
#define FLANN_POINT_SET

#include <Eigen/Dense>
using Eigen::Vector3d;

typedef struct Pt Pt;

class FLANNPointSet
{
public :
  std::vector<Vector3d>  *pts;
  
  FLANNPointSet()
  {
    xmin= 0;
    xmax = -10.0;
  }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return pts->size(); }

  // Returns the distance between the vector "p1[0:size-1]" and the data point 
  //with index "idx_p2" stored in the class:
  inline double kdtree_distance(const double *p1, const size_t idx_p2,size_t 
size) const
  {
    const double d0=p1[0]-(*pts)[idx_p2].x();
    const double d1=p1[1]-(*pts)[idx_p2].y();
    //const double d2=p1[2]-(*pts)[idx_p2].z();
    //return d0*d0+d1*d1+d2*d2;
    return d0*d0+d1*d1;
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate 
  //value, the  "if/else's" are actually solved at compile time.
  inline double kdtree_get_pt(const size_t idx, int dim) const
  {
    if (dim==0) return (*pts)[idx].x();
    else if (dim==1) return (*pts)[idx].y();
    else return (*pts)[idx].z();
  }

  // Optional bounding-box computation: return false to default to a standard 
  //bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned 
  //in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 
  //for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX &bb) const { return false; }

 
  double xmin,ymin,zmin;
  double xmax,ymax,zmax;
  
  void addPoint(Vector3d &mp)
  {
      pts->push_back(mp);

      if(xmax - xmin < -1.0)
      {
        xmin = xmax = mp.x();
        ymin = ymax = mp.y();
        zmin = zmax = mp.z();
      }
      else
      {
        xmax = xmax > mp.x() ? xmax : mp.x();
        ymax = ymax > mp.y() ? ymax : mp.y();
        zmax = zmax > mp.z() ? zmax : mp.z();
        xmin = xmin < mp.x() ? xmin : mp.x();
        ymin = ymin < mp.y() ? ymin : mp.y();
        zmin = zmin < mp.z() ? zmin : mp.z();
      }
  }
  
  void setPoints(std::vector<Vector3d> *locpoints)
  {
      pts = locpoints;
  }
  
  void addPoints(MatrixXd &v)
  {
      for(int i = 0; i <v.rows(); ++i)
      {
        pts->push_back(v.row(i));
      }
  }
  //void addSample(IndexedPoint &mp)
  //{
  //    Pt p;
  //    p.x = mp.x();
  //    p.y = mp.y();
  //    p.z = mp.z();
  //    p.nx= mp.nx();
  //    p.ny= mp.ny();
  //    p.nz= mp.nz();
  //    p.t1x = mp.t1x();
  //    p.t1y = mp.t1y();
  //    p.t1z = mp.t1z();
  //    p.t2x = mp.t2x();
  //    p.t2y = mp.t2y();
  //    p.t2z = mp.t2z();
  //    p.lambda0 = mp.lambda0();
  //    p.lambda1 = mp.lambda1();
  //
  //    if(xmax - xmin < -1.0)
  //    {
  //      xmin = xmax = p.x;
  //      ymin = ymax = p.y;
  //      zmin = zmax = p.z;
  //    }
  //    else
  //    {
  //      xmax = xmax > p.x ? xmax : p.x;
  //      ymax = ymax > p.y ? ymax : p.y;
  //      zmax = zmax > p.z ? zmax : p.z;
  //      xmin = xmin < p.x ? xmin : p.x;
  //      ymin = ymin < p.y ? ymin : p.y;
  //      zmin = zmin < p.z ? zmin : p.z;
  //    }
  //    pts.push_back(p);
  //}
};

#endif
