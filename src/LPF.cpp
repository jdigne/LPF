/**
 * This file is part of lpf
 *
 * @file LPF.cpp
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
#include "LPF.h"

#include <cmath>
#include <cassert>
#include <Eigen/Dense>
#include "utils.h"
#include "nanoflann.hpp"
#include "FLANNPointSet.h"

typedef nanoflann::KDTreeSingleIndexAdaptor<
 nanoflann::L2_Simple_Adaptor<double, FLANNPointSet > ,
 FLANNPointSet,
 3 /* dim */
 > KdTree;

using namespace std;
using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Matrix2d;
using Eigen::Vector2d;
using Eigen::JacobiSVD;

LPF::LPF()
{
  m_nbins = 0;
  m_irreg = false;
}

LPF::LPF(unsigned int nbins)
{
  m_nbins = nbins;
  m_v.resize(nbins,3);
}

LPF::LPF(std::vector<Vector3d>& qpoints, std::vector<Vector3d> &locpoints)
{
  m_nbins = qpoints.size();
  m_irreg = true;
  m_v.resize(m_nbins,3);
  m_v.setZero();
  
  //saving target zone
  m_target.insert(m_target.end(), locpoints.begin(), locpoints.end());
  
  FLANNPointSet set;
  set.setPoints(&locpoints);
  KdTree *kdtree = new KdTree(3, set, nanoflann::KDTreeSingleIndexAdaptorParams(10) );
  kdtree->buildIndex();

#pragma omp parallel for
  for(int i = 0; i < qpoints.size(); ++i)
  {
    nanoflann::KNNResultSet<double> resultSet(1);
    size_t ret_index;
    double out_dist_sqr;
    resultSet.init(&ret_index, &out_dist_sqr);

    double query_pt[3] = {qpoints[i].x(), qpoints[i].y(), qpoints[i].z()};
    kdtree->findNeighbors(resultSet, &query_pt[0],
                        nanoflann::SearchParams(10));

    Vector3d &p = locpoints[ret_index];
    m_v.row(i) = p - qpoints[i]; 
    
//    assert(!isnan( m_v(i,0) ));
//    assert(!isnan( m_v(i,1) ));
//    assert(!isnan( m_v(i,2) ));
  }
  delete kdtree;
}

LPF::~LPF()
{
  m_nbins = 0;
  m_irreg = false;
  m_v.setZero();
}


bool LPF::write(ofstream& f)
{
  f << m_v.col(0).transpose() << "\t" << m_v.col(1).transpose() << "\t" << m_v.col(2).transpose() << std::endl;
  
  return true;

}

bool LPF::writeBin(std::ofstream& f)
{
   f.write( (char *) m_v.data(), m_v.rows() * m_v.cols() * sizeof(m_v(0,0)) );
   //f.write(reinterpret_cast<const char*>(&(m_vx[0])), m_vx.size()*sizeof(m_vx[0]));
   //f.write(reinterpret_cast<const char*>(&(m_vy[0])), m_vy.size()*sizeof(m_vy[0]));
   //f.write(reinterpret_cast<const char*>(&(m_vz[0])), m_vz.size()*sizeof(m_vz[0]));
   return true;
}

bool LPF::read(unsigned int nbins, std::ifstream& f)
{
  m_nbins = nbins;
  //m_vx.assign(nbins, 0);
  //m_vy.assign(nbins, 0);
  //m_vz.assign(nbins, 0);
  m_v.resize(nbins,3);

  for(size_t i = 0; i < nbins; ++i)
  {
    f >> m_v(i,0);
  }
  for(size_t i = 0; i < nbins; ++i)
  {
    f >> m_v(i,1);
  }
  for(size_t i = 0; i < nbins; ++i)
  {
    f >> m_v(i,2);
  }
  return true;
}

bool LPF::readBin(unsigned int nbins, ifstream& f)
{
  m_v.resize(nbins,3);
  f.read(reinterpret_cast<char*>(m_v.data()), m_v.rows() * m_v.cols() * sizeof(m_v(0,0)) );
  return true;
}


bool LPF::reconstruct(const Vector3d& seed,
                      const Vector3d& t1,
                      const Vector3d& t2,
                      const Vector3d& n,
                      const double radius,
                      const vector< Vector3d >& qpoints,
                      vector< Vector3d >& recpoints)
{
  m_nbins = qpoints.size();
  recpoints.assign(qpoints.size(), Vector3d(0,0,0));

#pragma omp parallel for
  for(size_t i = 0; i < qpoints.size(); ++i)
  {
    double x,y,z;
    double dx = radius * (qpoints[i].x() + m_v(i,0));
    double dy = radius * (qpoints[i].y() + m_v(i,1));
    double dz = radius * (qpoints[i].z() + m_v(i,2));
    x = seed(0) + dx * t1(0) + dy * t2(0) + dz * n.x();
    y = seed(1) + dx * t1(1) + dy * t2(1) + dz * n.y();
    z = seed(2) + dx * t1(2) + dy * t2(2) + dz * n.z();
    recpoints[i] << x, y , z;
  }
  return true;
}


bool LPF::save_target_zone(const Vector3d& seed,
                      const Vector3d& t1,
                      const Vector3d& t2,
                      const Vector3d& n,
                      const double radius,
                      vector< Vector3d >& recpoints)
{
  recpoints.assign(m_target.size(), Vector3d(0,0,0));

#pragma omp parallel for
  for(size_t i = 0; i < m_target.size(); ++i)
  {
    Vector3d &p = m_target[i];
    double x = seed(0) + radius * (p(0) * t1(0) + p(1) * t2(0) + p(2) * n.x());
    double y = seed(1) + radius * (p(0) * t1(1) + p(1) * t2(1) + p(2) * n.y());
    double z = seed(2) + radius * (p(0) * t1(2) + p(1) * t2(2) + p(2) * n.z());
    recpoints[i] << x, y , z;
  }
  return true;
}


unsigned int LPF::getSizeDesc() const
{
  return m_v.rows()*3;
}

bool LPF::vectorize(size_t i, double &val) const
{
  if(i < m_nbins)
  {
    val = m_v(i,0);
    return true;
  }

  i = i - m_nbins;  
  if(i < m_nbins)
  {
    val = m_v(i,1);
    return true;
  }

  i = i - m_nbins;  
  if(i < m_nbins)
  {
    val = m_v(i,2);
    return true;
  }

  return false;
}

bool LPF::readMatrix(const VectorXd& X)
{
  if(X.size() % 3 != 0)
    return false;

  unsigned int size = X.size()/3;
  
  if(size != m_nbins)
    return false;
  
  for(size_t i = 0; i < size; ++i)
  {
      m_v(i,0) = X(i);
      m_v(i,1) = X(i+size);
      m_v(i,2) = X(i+2*size);
  }
  return true;
}


void LPF::recompute(const MatrixXd& qpoints)
{ 
  FLANNPointSet set;

  set.setPoints(&m_target);
  KdTree *kdtree = new KdTree(3, set, nanoflann::KDTreeSingleIndexAdaptorParams(10) );
  kdtree->buildIndex();

#pragma omp parallel for
  for(int i = 0; i < qpoints.rows(); ++i)
  {
    nanoflann::KNNResultSet<double> resultSet(1);
    size_t ret_index;
    double out_dist_sqr;
    resultSet.init(&ret_index, &out_dist_sqr);

    double query_pt[3] = {qpoints(i,0), qpoints(i,1), qpoints(i,2)};
    kdtree->findNeighbors(resultSet, &query_pt[0],
                        nanoflann::SearchParams(10));

    Vector3d &p = m_target[ret_index];
    m_v.row(i) = p.transpose() - qpoints.row(i); 
  }
  delete kdtree;
}


void LPF::recompute(const VectorXd& qpoints)
{
  FLANNPointSet set;

  
  set.setPoints(&m_target);
  KdTree *kdtree = new KdTree(3, set, nanoflann::KDTreeSingleIndexAdaptorParams(10) );
  kdtree->buildIndex();

#pragma omp parallel for
  for(int i = 0; i < m_v.rows(); ++i)
  {
    nanoflann::KNNResultSet<double> resultSet(1);
    size_t ret_index;
    double out_dist_sqr;
    resultSet.init(&ret_index, &out_dist_sqr);

    double query_pt[3] = {qpoints(i), qpoints(i + m_nbins), qpoints(i + 2*m_nbins)};
    kdtree->findNeighbors(resultSet, &query_pt[0],
                        nanoflann::SearchParams(10));

    Vector3d &p = m_target[ret_index];
    m_v.row(i) = p;
    m_v(i,0) -= query_pt[0];
    m_v(i,1) -= query_pt[1];
    m_v(i,2) -= query_pt[2];
  }
  delete kdtree;
}

void LPF::update(const VectorXd& val)
{
  int i,j,k;
  for(i = 0, j=m_nbins, k =2*m_nbins; i < m_nbins; ++i, ++j, ++k)
  {
    m_v(i,0) = val(i);
    m_v(i,1) = val(j);
    m_v(i,2) = val(k);
  }
}

void LPF::tovector(VectorXd& val) const
{
  int i,j,k;
  val.resize(3*m_nbins);
  for(i = 0, j=m_nbins, k =2*m_nbins; i < m_nbins; ++i, ++j, ++k)
  {
    val(i) = m_v(i,0);
    val(j) = m_v(i,1);
    val(k) = m_v(i,2);
  }
}



double LPF::vote(Eigen::Vector3d& point,
                  const Eigen::Vector3d seed,
                  const Eigen::Vector3d& t1,
                  const Eigen::Vector3d& t2,
                  const Eigen::Vector3d& normal,
                  const vector< Eigen::Vector3d >& qpoints,
                  const double radius)
{
  Vector3d loc = (point - seed)/radius;
  double xl = loc.dot(t1);
  double yl = loc.dot(t2);
  double zl = loc.dot(normal);
  double d3 = (point -seed).squaredNorm();
  
  assert(t1.dot(t2)<1e-8);
  assert(t1.dot(normal)<1e-8);
  assert(normal.dot(t2)<1e-8);
  
  loc(0) = xl;
  loc(1) = yl;
  loc(2) = zl;
  
  Vector3d temp(0, 0, 0);
  double sumw = 0.0;
  
  double sqmaxrad = 4*radius*radius /qpoints.size()*3.14;
  double invsigma2 = -0.5/sqmaxrad;
  
  for(int i = 0; i <qpoints.size(); ++i)
  {
    Vector3d t = m_v.row(i) + qpoints[i].transpose();
    double d2 = sqr(loc.x()- qpoints[i].x())+sqr(loc.y() - qpoints[i].y());
    if(d2 > sqmaxrad)
      continue;
    
    double w = exp(invsigma2*d2);
    sumw += w;
    temp += t *w;
  }
  
  if(sumw < 1e-16)
    return 0.0;
  
  temp = temp * radius/sumw;
  
  point = seed + (temp(0) * t1 + temp(1) * t2 + temp(2) * normal);
  return sumw;
}


void LPF::optimize_efficiency(vector< Vector3d >& locpoints,
                              const MatrixXd &genpat,
                              Vector3d &seed,
                              Vector3d& t1,
                              Vector3d& t2,
                              Vector3d& n,
                              double radius)
{
     double pareff0 = parallelism_measure();
     MatrixXd p  = m_v + genpat;
     
     //Minimize m_v norm
     Matrix3d R;
     Vector3d t;
     
     svd_transform_estimate_3(genpat, p, R, t);
     
     m_v = (p.rowwise() - t.transpose())*R - genpat;
     
     //update seed position
     seed = seed + radius*(t(0) * t1 + t(1) * t2 + t(2) * n);
     
     Vector3d t1n = R(0,0) * t1 + R(1,0) * t2 + R(2,0) * n;
     Vector3d t2n = R(0,1) * t1 + R(1,1) * t2 + R(2,1) * n;
     Vector3d nn  = R(0,2) * t1 + R(1,2) * t2 + R(2,2) * n;
     
     t1 = t1n;
     t2 = t2n;
     n = nn;
     
     update_target_zone(R, t);
     
     recompute(genpat);
     
     double pareff1 = parallelism_measure();
}

double LPF::parallelism_measure()
{
  return sqrt((m_v.col(0).cwiseProduct(m_v.col(0)) + m_v.col(1).cwiseProduct(m_v.col(1))).mean());
}

void LPF::update_target_zone(Matrix3d& R, Vector3d& t)
{
#pragma omp parallel for
  for(int i = 0; i < m_target.size(); ++i)
  {
    Vector3d v = m_target[i];
    m_target[i] = R.transpose()*(v - t);
  }
}


void LPF::update_target_zone_after_ICP(Matrix3d& R, Vector3d& t)
{
#pragma omp parallel for
  for(int i = 0; i < m_target.size(); ++i)
  {
    m_target[i] = R * m_target[i];
    m_target[i] += t;
  }
}

