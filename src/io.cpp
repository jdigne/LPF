/**
 * This file is part of lpf
 *
 * @file io.cpp
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
#include "io.h"
#include "utils.h"
#include "utils_octree.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>

using namespace std;

IO::IO()
{

}


IO::~IO()
{

}


bool IO::readPoints(const char* filename, vector< Vector3d >& data)
{
  unsigned int nlines, nwords = 0;
  bool ok = getnlines(filename, nlines, nwords);
  
  if(! ok )
    return false;
  
  ifstream in;
  in.open(filename);

  double x, y, z;
  data.reserve(nlines);
  
  string line;
  getline(in,line);
  istringstream line_in(line);
  line_in >> x >> y >> z;
  
  for(size_t i = 0; i< nlines; ++i)
  { 
    data.push_back(Vector3d(x,y,z));
    
    getline(in,line);
    istringstream iss(line);
    iss >> x >> y >> z;
  }
  
  in.close();
  return true;
}

bool IO::readPoints(const char* filename, const double min_radius, vector< Vector3d >& data, Vector3d& octree_origin, double& size, unsigned int& depth)
{
  unsigned int nlines, nwords = 0;
  bool ok = getnlines(filename, nlines, nwords);
  
  if(! ok )
    return false;
  
  ifstream in;
  in.open(filename);

  double x, y, z;
  data.reserve(nlines);
  
  string line;
  getline(in,line);
  istringstream line_in(line);
  line_in >> x >> y >> z;
  
  octree_origin(0) = x;
  octree_origin(1) = y;
  octree_origin(2) = z;
  double lx = x;
  double ly = y;
  double lz = z;
  
  for(size_t i = 0; i< nlines; ++i)
  { 
    data.push_back(Vector3d(x,y,z));
    
    octree_origin(0) = octree_origin(0) < x ? octree_origin(0) : x;
    octree_origin(1) = octree_origin(1) < y ? octree_origin(1) : y;
    octree_origin(2) = octree_origin(2) < z ? octree_origin(2) : z;
    
    lx = lx > x ? lx : x;
    ly = ly > y ? ly : y;
    lz = lz > z ? lz : z;
    getline(in,line);
    istringstream iss(line);
    iss >> x >> y >> z;
  }
  
  lx -= octree_origin(0);
  ly -= octree_origin(1);
  lz -= octree_origin(2);
  
  size = lx > ly ? lx : ly;
  size = size > lz ? size : lz;
  
  double nsize = 1.2 * size;
  double margin=0;
  
  if(min_radius > 0)
  {
    depth = (unsigned int)ceil( log2( nsize / (min_radius) ));
    double adapted_size = pow2(depth) * min_radius;
    margin = 0.5 * (adapted_size - nsize);
    size = adapted_size;
  }
  else
  {
    margin = 0.1 * size;
  }
  
  octree_origin = octree_origin.array() - margin;
  in.close();
  return true;
}



bool IO::savePoints(const char* filename, const vector< Vector3d >& data)
{
  ofstream out;
  out.open(filename);
  if(!out)
    return false;
  out.precision( std::numeric_limits<double>::digits10 + 1);
  
  for(size_t i = 0; i < data.size(); ++i)
  {
    out << data[i](0) <<"\t"<<data[i](1)<<"\t"<<data[i](2)<< "\n";
  }
  out.close();
  return true;

}

bool IO::saveMatrix(const char* filename, const MatrixXd& M)
{
  
  ofstream out;
  out.open(filename);
  if(!out)
    return false;
  out.precision( std::numeric_limits<double>::digits10 + 1);
  out << M;
  return true;
}



bool IO::readParameterizedPoints(const char* filename,
                                 const double min_radius,
                                 std::vector< Vector3d >& data,
                                 std::vector< Vector3d >& t1s,
                                 std::vector< Vector3d >& t2s,
                                 std::vector< Vector3d >& normals,
                                 Vector3d &octree_origin,
                                 double &size,
                                 unsigned int &depth)
{
  unsigned int nlines, nwords = 0;
  bool ok = getnlines(filename, nlines, nwords);

  if(! ok )
    return false;
  
  if(nwords <12)
  {
    std::cout<<nwords<<" doubles per line, not enough to represent a point"
              <<"its tangent vectors and normal. Exiting"<<std::endl;
    return false;
  }

  ifstream in;
  in.open(filename);

  double x,y,z;
  double nx,ny,nz;
  double t1x,t1y,t1z;
  double t2x,t2y,t2z;
  data.reserve(nlines);
  t1s.reserve(nlines);
  t2s.reserve(nlines);
  normals.reserve(nlines);
  
  
  string line;
  getline(in,line);
  istringstream line_in(line);
  line_in >> x >> y >> z;
  line_in >> t1x >> t1y >> t1z;
  line_in >> t2x >> t2y >> t2z;
  line_in >> nx  >> ny  >> nz;
  
  octree_origin(0) = x;
  octree_origin(1) = y;
  octree_origin(2) = z;
  double lx = x;
  double ly = y;
  double lz = z;
  
  for(size_t i = 0; i< nlines; ++i)
  {
    data.push_back(Vector3d(x, y, z));
    t1s.push_back(Vector3d(t1x, t1y, t1z));
    t2s.push_back(Vector3d(t2x, t2y, t2z));
    normals.push_back(Vector3d(nx, ny, nz));
    
    octree_origin(0) = octree_origin(0) < x ? octree_origin(0) : x;
    octree_origin(1) = octree_origin(1) < y ? octree_origin(1) : y;
    octree_origin(2) = octree_origin(2) < z ? octree_origin(2) : z;
    
    lx = lx > x ? lx : x;
    ly = ly > y ? ly : y;
    lz = lz > z ? lz : z;
    
    getline(in,line);
    line_in.str(line);
    line_in >> x >> y >> z;
    line_in >> t1x >> t1y >> t1z;
    line_in >> t2x >> t2y >> t2z;
    line_in >> nx  >> ny  >> nz;   
  }
  in.close();
  
  size = lx - octree_origin.x();
  double sizet = ly - octree_origin.y();
  size = size > sizet ? size : sizet;
  
  sizet = lz - octree_origin.z();
  size = size > sizet ? size : sizet;
  
  octree_origin = octree_origin.array() - 0.05*size;
  size = 1.1*size;
  
  return true;
}


bool IO::getnlines(const char* filename, unsigned int& nlines, unsigned int& nwords)
{
  ifstream in;
  in.open(filename);
  
  if(!in)
  {
    cerr<<filename<<"could not be opened. Exiting..."<<endl;
    return false;
  }
  
  string line;
  getline(in,line);
  istringstream line_in(line);
  string word;
  nwords = 0;
  while (line_in>> word)
    nwords++;
  
  //count number of lines
  nlines = 1;
  while(getline(in,line))
    ++nlines;

  in.close();
  return true;
}

bool IO::saveParams(const char* filename,
                    const vector< Vector3d >& t1s,
                    const vector< Vector3d >& t2s,
                    const vector< Vector3d >& normals)
{
    std::ofstream f;
    f.open(filename, std::ios::out);
    if(!f)
    {
      std::cout<<"Could not open file "<<filename<<" for saving the parameterizations. Exiting."<<std::endl;
      return false;
    }
    f.precision( std::numeric_limits<double>::digits10 + 1);
    
    typename std::vector<Vector3d>::const_iterator t1;
    typename std::vector<Vector3d>::const_iterator t2;
    typename std::vector<Vector3d>::const_iterator n;
    for(t1 = t1s.begin(), t2 = t2s.begin(), n = normals.begin(); t1 != t1s.end(); ++t1, ++t2, ++n)
    {
      const Vector3d &tt1 = *t1; 
      const Vector3d &tt2 = *t2; 
      const Vector3d &nn = *n; 
      f << tt1(0) << "\t" << tt1(1) << "\t" << tt1(2) << "\t"
        << tt2(0) << "\t" << tt2(1) << "\t" << tt2(2) << "\t"
        <<  nn(0) << "\t" << nn(1)  << "\t" << nn(2)  << "\n";
    }
    f.close();
    return true;
}


bool IO::readSeedParams(const char* filename, vector< Vector3d >& t1s, vector< Vector3d >& t2s, vector< Vector3d >& normals)
{
  unsigned int nlines, nwords = 0;
  bool ok = getnlines(filename, nlines, nwords);

  if(! ok )
    return false;
  ifstream in;
  in.open(filename);
  
  double x,y,z;
  double nx,ny,nz;
  double t1x,t1y,t1z;
  double t2x,t2y,t2z;
  t1s.reserve(nlines);
  t2s.reserve(nlines);
  normals.reserve(nlines);
  
  string line;
  getline(in,line);
  istringstream line_in(line);
  line_in >> t1x >> t1y >> t1z;
  line_in >> t2x >> t2y >> t2z;
  line_in >> nx  >> ny  >> nz;
  
  for(size_t i = 0; i< nlines; ++i)
  {
    t1s.push_back(Vector3d(t1x, t1y, t1z));
    t2s.push_back(Vector3d(t2x, t2y, t2z));
    normals.push_back(Vector3d(nx, ny, nz));
    
    getline(in,line);
    line_in.str(line);
    line_in >> x >> y >> z;
    line_in >> t1x >> t1y >> t1z;
    line_in >> t2x >> t2y >> t2z;
    line_in >> nx  >> ny  >> nz;   
  }
  in.close();
  
  return true;
}

bool IO::readMatrix(const char* filename, MatrixXd &M)
{
  unsigned int nrows,ncols;
  bool ok = getnlines(filename, nrows, ncols);
  if(!ok)
    return false;
  
  M.resize(nrows, ncols);
  
  ifstream f;
  f.open(filename);
  for(int i = 0; i < nrows; ++i)
    for(int j = 0; j < ncols ; ++j)
    {
      f>>M(i,j);
    }
  f.close();
  
  return true;
}
