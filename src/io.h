/**
 * This file is part of lpf
 *
 * @file io.h
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
#ifndef IO_H
#define IO_H
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

using Eigen::Vector3d;
using Eigen::MatrixXd;

class IO
{
public :
  IO();
  
  ~IO();
  
  
  
  static bool readPoints(const char *filename, const double min_radius, std::vector<Vector3d> &data, Vector3d &octree_origin, double &size, unsigned int &depth);
  
  
  static bool readPoints(const char *filename, std::vector<Vector3d> &data);
  
  static bool readSeedParams(const char *filename, std::vector<Vector3d> &t1s, std::vector<Vector3d> &t2s, std::vector<Vector3d> &normals);
  
  static bool readMatrix(const char*filename, MatrixXd &M);
  
  static bool saveMatrix(const char *filename, const Eigen::MatrixXd &M);
  
  static bool savePoints(const char *filename, const std::vector<Vector3d> &data);
  
  static bool readParameterizedPoints(const char* filename, const double min_radius, std::vector<Vector3d> &data, std::vector<Vector3d> &normals, std::vector<Vector3d> &t1s, std::vector<Vector3d> &t2s, Vector3d &octree_origin, double &size, unsigned int &depth);
  
  template<class T>
  static bool saveDesc(const char *filename, std::vector<T*> &data)
  {
    std::ofstream f;
    f.open(filename, std::ios::out);
    if(!f)
    {
      std::cout<<"Could not open file "<<filename<<" for saving the descriptors. Exiting."<<std::endl;
      return false;
    }
    f.precision( std::numeric_limits<double>::digits10 + 1);
    
    typename std::vector<T*>::const_iterator it;
    for(it = data.begin(); it !=data.end(); ++it)
    {
      (*it)->write(f);
    }
    
    f.close();
    return true;
  }
  
  template<class T>
  static bool saveDescBin(const char* filename, std::vector<T*> &data)
  {
    //ofstream out;
    //out.open(filename);
    //out.precision( numeric_limits<double>::digits10 + 1);

   std::ofstream f;
   f.open(filename, std::ios::out | std::ios::binary);
   if(!f)
   {
     std::cout<<"Could not open file "<<filename<<" for saving the descriptors. Exiting."<<std::endl;
     return false;
   }
   typename std::vector<T*>::const_iterator it;
   unsigned int size_desc = data[0]->size_desc;
   unsigned int ndesc = data.size();

   f.write(reinterpret_cast<const char*>(&size_desc), sizeof(size_desc));
   f.write(reinterpret_cast<const char*>(&ndesc), sizeof(ndesc));

   for(it = data.begin(); it !=data.end(); ++it)
   {
     (*it)->write_bin(f);
   }
   
   f.close();
   return true;
}
  
  static bool saveParams(const char *filename,
                         const std::vector<Vector3d> &t1s,
                         const std::vector<Vector3d> &t2s,
                         const std::vector<Vector3d> &normals);

  static bool getnlines(const char *filename, unsigned int &nlines, unsigned int &nwords);

};

#endif
