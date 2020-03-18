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
/**
 * @file OctreeIterator.h
 * @brief defines methods for getting range neighborhoods in the octree
 * @author Julie Digne
 * @date 2012/10/12
 * @copyright This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OCTREE_ITERATOR_H
#define OCTREE_ITERATOR_H

#include<cstdlib>
#include <map>
#include <vector>
#include "Octree.h"
#include "OctreeNode.h"
#include <cmath>
#include<cassert>
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <list>

using Eigen::Vector3d;

/**
 * @class TOctreeIterator
 * @brief Defines methods to access range neighbors of points
 *
 * Templated class that contains all methods for neighborhoods
 * queries on an octree
 */
template<class T>
class TOctreeIterator
{
    public : //some typedefs

        typedef std::vector<T*> Neighbor_star_list;
        typedef std::set<T*> Exception_set;
        typedef std::map<double, T*> Neighbor_star_map;
        typedef std::vector<double> Distance_list;

    private ://class members
        /** @brief active depth*/
        unsigned int m_activeDepth;

        /** @brief active adius*/
        double m_radius;

        /** @brief square active radius*/
        double m_sqradius;

        /** @brief index of the set to retriev neighbors from
         *(must be 0,1 or 2)
         */
        unsigned int m_setIndex;

        /** @brief TOctree the iterator refers to*/
        TOctree<T> *m_octree;

    public ://constructor+destructor

        /** @brief constructor*/
        TOctreeIterator<T>();

        /** @brief constructor
         * @param octree octree containing the point set to process
         */
        TOctreeIterator<T>(TOctree<T> * octree);

        /** @brief destructor*/
        ~TOctreeIterator<T>();

    public : //accessors modifiers

        /** @brief set radius and active depth accordingly
         * @param radius radius to set
         * @return false if the radius was not coherent with the
         * structure (> bounding box size)
         */
        bool setR(double radius);

        /** @brief get active radius
         * @return active radius
         */
        double getR() const;

        /** @brief get square active radius
         * @return square radius
         */
        double getSquareR() const;

        /** @brief set depth and radius accordingly (smallest radius
         * for this depth)
         * @param depth depth to set
         * @return false if the depth was not coherent with the octree
         * depth
         */
        bool setDepth(unsigned int depth);

        /** @brief get active depth
         * @return activeDepth
         */
        unsigned int getDepth() const;

        /** @brief get index of the set
         * @return index of the set that is being searched for
         */
        unsigned int getSetIndex() const;

        /** @brief set index of the set
         * @param index of the set that is being searched for
         */
        void setSetIndex(unsigned int index);

    public : //getting neighbors

        /** @brief get star-neighbors of a given point
         *@param query query point
         *@param[out] neighbors list of neighbors to be filled
         * by the method
         *@return number of neighbors
         */
        unsigned int getNeighbors(const Vector3d &query,
                Neighbor_star_list &neighbors) const;

        /** @brief get star-neighbors of a given point
         *@param query query point
         *@param[out] neighbors list of neighbors to be filled
         * by the method
         *@param[out] distances list of neighbors distances in
         * the same order
         *@return number of neighbors
         */
        unsigned int getNeighbors(const Vector3d &query,
                Neighbor_star_list &neighbors,
                Distance_list &distances) const;

        /** @brief get star neighbors of a given point when
         * the node containing that point is known
         *@param query query point
         *@param query_node node containing the query point
         *@param[out] neighbors list of neighbors to be filled
         * by the method
         *@return number of neighbors
         */
        unsigned int getNeighbors(const Vector3d &query,
                TOctreeNode<T> *query_node,
                Neighbor_star_list &neighbors) const;

        /** @brief get star neighbors of a given point when
         * the node containing that point is known
         *@param query query point
         *@param query_node node containing the query point
         *@param[out] neighbors list of neighbors to be filled
         * by the method
         *@param[out] distances list of neighbors distances in
         * the same order
         *@return number of neighbors
         */
        unsigned int getNeighbors(const Vector3d &query,
                TOctreeNode<T> *query_node,
                Neighbor_star_list &neighbors,
                Distance_list &distances) const;

        /** @brief get neighbors of a given point sorted by their distances
         *@param query query point
         *@param[out] neighbors map of neighbors to be filled by the method
         *@return number of neighbors
         */
        unsigned int getSortedNeighbors(const Vector3d &query,
                Neighbor_star_map &neighbors) const;

        /** @brief get neighbors of a given point when the node containing
         * that point is known neighbors are sorted by their distances to
         * the query
         *@param query query point
         *@param node node containing the query point
         *@param[out] neighbors map of neighbors to be filled by the method
         *@return number of neighbors
         */
        unsigned int getSortedNeighbors(const Vector3d &query,
                TOctreeNode<T> *node,
                Neighbor_star_map &neighbors) const;

        /** @brief Look in a ball centered at query point if
         * there is any other
         * point than those given in the parameter set
         *@param query center point
         *@param exceptions set of elements that are allowed
         * in the neighborhood
         *@return false if the neighborhood contains other elements than
         *those in the exception set
         */
        bool containsOnly(const Vector3d &query,
                const Exception_set &exceptions) const;

        /** @brief Look in a ball centered at query point if there is any
         * other point than those given in the parameter set
         *@param query center point
         *@param query_node node containing the query point
         *@param exceptions set of elements that are allowed in the
         * neighborhood
         *@return false if the neighborhood contains other elements
         * than those in the exception set
         */
        bool containsOnly(const Vector3d &query, TOctreeNode<T> *query_node,
                const Exception_set &exceptions) const;

    private :

        /**
         * @brief explore a node to find neighbors of a point.
         * @param node (node to explore)
         * @param query_point (center of the neighborhood)
         * @param neighbors deque of points (at the end: neighbors of pt)
         */
        void explore(TOctreeNode<T> *node, const Vector3d &query_point,
                Neighbor_star_list &neighbors) const;

        /**
         * @brief explore a node to look at neighbors of a point. Stops if
         * one of those neighbors is not in the exception set
         * @param node (node to explore)
         * @param query_point (center of the neighborhood)
         * @param check false if the neighborhood contains other elements
         * @param exceptions set of elements that are allowed in the
         * neighborhood
         */
        void explore(TOctreeNode<T> *node, const Vector3d &query_point,
                const Exception_set &exceptions, bool &check) const;

        /**
         * @brief explore a node to find neighbors of a point.
         * @param node (node to explore)
         * @param query_point (center of the neighborhood)
         * @param neighbors deque of points (at the end: neighbors of pt)
         * @param distances distance to the neigbors
         */
        void explore(TOctreeNode<T> *node, const Vector3d &query_point,
                Neighbor_star_list &neighbors,
                Distance_list &distances) const;

        /** @brief explore a node to find neighbors of a point and
         * sort them according to their distance
         * @param node (node to explore)
         * @param query_point (center of the neighborhood)
         * @param neighbors map of points sorted by their distances
         * to the query point
         */
        void exploreSort(TOctreeNode<T> *node, const Vector3d &query_point,
                Neighbor_star_map &neighbors) const;

        /**
         * @brief follow the path given by loc-codes beginning at node
         * @param node pointer to the node where the path begins
         * @param xLocCode x locational code of the path
         * @param yLocCode y locational code of the path
         * @param zLocCode z locational code of the path
         * @param k max level to look for
         */
        void traverseToLevel(TOctreeNode<T> **node, unsigned int xLocCode,
                unsigned int yLocCode, unsigned int zLocCode,
                unsigned int k) const;

        /**
         * @brief get left neighbor code (along x-axis) of a cell
         * @param cell
         * @return unsigned int left neighbor xloc
         */
        unsigned int getXLeftCode(TOctreeNode<T> *cell) const ;

        /**
         * @brief get left neighbor code (along y-axis) of a cell
         * @param cell
         * @return unsigned int left neighbor yloc
         */
        unsigned int getYLeftCode(TOctreeNode<T> *cell) const;

        /**
         * @brief get left neighbor code (along z-axis) of a cell
         * @param cell
         * @return unsigned int left neighbor zloc
         */
        unsigned int getZLeftCode(TOctreeNode<T> *cell) const;

        /**
         * @brief get right neighbor code (along x-axis) of a cell
         * @param cell
         * @return unsigned int right neighbor xloc
         */
        unsigned int getXRightCode(TOctreeNode<T> *cell) const;

        /**
         * @brief get right neighbor code (along y-axis) of a cell
         * @param cell
         * @return unsigned int right neighbor yloc code
         */
        unsigned int getYRightCode(TOctreeNode<T> *cell) const;

        /**
         * @brief get right neighbor code (along z-axis) of a cell
         * @param cell cell
         * @return unsigned int right neighbor zloc
         */
        unsigned int getZRightCode(TOctreeNode<T> *cell) const;

        /** @brief compute the code of a given position (useful for
         * locating points and traversing the octree
         * @param point point to locate
         * @param[out] xcode x locational code
         * @param[out] ycode y locational code
         * @param[out] zcode z locational code
         */
        void computeCode(const Vector3d &point, unsigned int &codx,
                unsigned int &cody, unsigned int &codz) const;

        /** @brief return a cell containing the point at active depth
         * @param point to locate
         * @return the node containing the point at active depth
         */
        TOctreeNode<T>* locateVector3dNode(const Vector3d &point) const;
};

    template<class T>
TOctreeIterator<T>::TOctreeIterator()
{
    m_octree = NULL;
    m_setIndex = 0;
}


    template<class T>
TOctreeIterator<T>::TOctreeIterator(TOctree<T>* octree)
{
    m_octree = octree;
    m_activeDepth = 0;
    m_radius = m_octree->getSize() / ((double)pow2(octree->getDepth() + 1));//smallest cell half-size
    m_sqradius = m_radius * m_radius;
    m_setIndex = 0;
}


    template<class T>
TOctreeIterator<T>::~TOctreeIterator()
{
    m_octree = NULL;
    m_setIndex = 0;
    m_radius = 0;
    m_sqradius = 0;
    m_activeDepth = 0;
}

    template<class T>
bool TOctreeIterator<T>::setDepth(unsigned int depth)
{
   unsigned int octreedepth= m_octree->getDepth();
   
   if(depth <= octreedepth)
   {
        m_activeDepth = depth;
        //m_radius = m_octree->getSize() / ((double)pow2(depth));
        m_radius = m_octree->getSize() / ((double)pow2(octreedepth - depth + 1));
        m_sqradius = m_radius * m_radius;
        return true;
    }
    return false;
}

template<class T>
unsigned int TOctreeIterator<T>::getDepth() const
{
    return m_activeDepth;
}


template<class T>
double TOctreeIterator<T>::getR() const
{
    return m_radius;
}

    template<class T>
bool TOctreeIterator<T>::setR(double radius)
{
    if(radius < m_octree->getSize())
    {
        m_radius = radius;
        m_sqradius = m_radius * m_radius;
        m_activeDepth = (unsigned int)(m_octree->getDepth()
                - floor( log2( m_octree->getSize() / (2.0*m_radius) )));
        return true;
    }
    return false;
}

template<class T>
double TOctreeIterator<T>::getSquareR() const
{
    return m_sqradius;
}


template<class T>
unsigned int TOctreeIterator<T>::getSetIndex() const
{
    return m_setIndex;
}

    template<class T>
void TOctreeIterator<T>::setSetIndex(unsigned int index)
{
    m_setIndex = index;
}


template<class T>
unsigned int TOctreeIterator<T>::getNeighbors(const Vector3d& query,
        Neighbor_star_list& neighbors) const
{
    TOctreeNode<T> *node = locateVector3dNode(query);
    unsigned int n = getNeighbors(query, node, neighbors);
    return n;
}

template<class T>
unsigned int TOctreeIterator<T>::getNeighbors(const Vector3d& query,
        Neighbor_star_list& neighbors,
        Distance_list &distances) const
{
    TOctreeNode<T> *node = locateVector3dNode(query);
    assert(node->isInside(query.x(), query.y(), query.z()));
    unsigned int n = getNeighbors(query, node, neighbors, distances);
    return n;
}



template<class T>
unsigned int TOctreeIterator<T>::getNeighbors(const Vector3d& query,
        TOctreeNode<T>* query_node,
        Neighbor_star_list& neighbors) const
{
    Vector3d octree_origin = m_octree->getOrigin();
    Vector3d node_origin = query_node->getOrigin();
    double node_size = query_node->getSize();
    double octree_size = m_octree->getSize();

    if(query_node->getDepth() ==  m_activeDepth)
    {
        //find neighboring nodes
        std::list<unsigned int> xloc, yloc, zloc;
        xloc.push_back(query_node->getXLoc());
        yloc.push_back(query_node->getYLoc());
        zloc.push_back(query_node->getZLoc());

        if((query.x() - m_radius  < node_origin.x())
                &&(query.x() - m_radius > octree_origin.x()))
            xloc.push_back(getXLeftCode(query_node));
        if((query.x() + m_radius > node_origin.x() + node_size)
                && (query.x() + m_radius <octree_origin.x() + octree_size))
            xloc.push_back(getXRightCode(query_node));


        if((query.y() - m_radius < node_origin.y())
                &&(query.y() - m_radius >octree_origin.y()))
            yloc.push_back(getYLeftCode(query_node));
        if((query.y() + m_radius > node_origin.y() + node_size)
                && (query.y() + m_radius < octree_origin.y() + octree_size))
            yloc.push_back(getYRightCode(query_node));


        if((query.z() - m_radius  < node_origin.z())
                &&(query.z() - m_radius > octree_origin.z()))
            zloc.push_back(getZLeftCode(query_node));
        if((query.z() + m_radius > node_origin.z() +node_size)
                && (query.z() + m_radius <octree_origin.z() + octree_size))
            zloc.push_back(getZRightCode(query_node));

        //look inside neighboring node
        std::list<unsigned int>::iterator xi,yi,zi;
        xi = xloc.begin();
        while(xi != xloc.end())
        {
            yi = yloc.begin();
            while(yi != yloc.end())
            {
                zi = zloc.begin();
                while(zi != zloc.end())
                {
                    TOctreeNode<T> *node=m_octree->getRoot();
                    traverseToLevel(&node, *xi, *yi, *zi, m_activeDepth);
                    if((node!=NULL)&&(node->getDepth() == m_activeDepth))
                        explore(node, query, neighbors);
                    ++zi;
                }
                ++yi;
            }
            ++xi;
        }
    }
    else
    {
        unsigned int s = query_node->getDepth();
        std::list<unsigned int> xloc, yloc, zloc;
        xloc.push_back(query_node->getXLoc());
        yloc.push_back(query_node->getYLoc());
        zloc.push_back(query_node->getZLoc());
        if((query.x() - m_radius  < node_origin.x() )
                && (query.x() - m_radius > octree_origin.x() ))
            xloc.push_back(getXLeftCode(query_node));
        if((query.x() + m_radius > node_origin.x() + node_size)
                && (query.x() + m_radius < octree_origin.x() +octree_size))
            xloc.push_back(getXRightCode(query_node));

        if((query.y() - m_radius  < node_origin.y() )
                && (query.y() - m_radius >octree_origin.y() ))
            yloc.push_back(getYLeftCode(query_node));
        if((query.y() + m_radius > node_origin.y() + node_size)
                && (query.y() + m_radius < octree_origin.y()+octree_size))
            yloc.push_back(getYRightCode(query_node));

        if((query.z() - m_radius  < node_origin.z() )
                && (query.z() - m_radius > octree_origin.z()))
            zloc.push_back(getZLeftCode(query_node));
        if((query.z() + m_radius > node_origin.z() + node_size)
                && (query.z() + m_radius <octree_origin.z() + octree_size))
            zloc.push_back(getZRightCode(query_node));

        //look inside neighboring nodes
        std::list<unsigned int>::iterator xi,yi,zi;
        xi = xloc.begin();
        while(xi != xloc.end())
        {
            yi = yloc.begin();
            while(yi != yloc.end())
            {
                zi = zloc.begin();
                while(zi != zloc.end())
                {
                    TOctreeNode<T> *node=m_octree->getRoot();
                    traverseToLevel(&node,*xi,*yi,*zi,s);
                    if( (node != NULL )&&( node->getDepth() == s))
                        explore(node,query,neighbors);
                    ++zi;
                }
                ++yi;
            }
            ++xi;
        }
    }
    return (int)neighbors.size();
}


template<class T>
unsigned int TOctreeIterator<T>::getNeighbors(const Vector3d& query,
        TOctreeNode<T>* query_node,
        Neighbor_star_list& neighbors,
        Distance_list &distances) const
{
    Vector3d octree_origin = m_octree->getOrigin();
    Vector3d node_origin = query_node->getOrigin();
    double node_size = query_node->getSize();
    double octree_size = m_octree->getSize();

    if(query_node->getDepth() ==  m_activeDepth)
    {
        //find neighboring nodes
        std::list<unsigned int> xloc, yloc, zloc;
        xloc.push_back(query_node->getXLoc());
        yloc.push_back(query_node->getYLoc());
        zloc.push_back(query_node->getZLoc());

        if((query.x() - m_radius  < node_origin.x())
                &&(query.x() - m_radius > octree_origin.x()))
            xloc.push_back(getXLeftCode(query_node));
        if((query.x() + m_radius > node_origin.x() + node_size)
                && (query.x() + m_radius <octree_origin.x() + octree_size))
            xloc.push_back(getXRightCode(query_node));

        if((query.y() - m_radius < node_origin.y())
                &&(query.y() - m_radius >octree_origin.y()))
            yloc.push_back(getYLeftCode(query_node));
        if((query.y() + m_radius > node_origin.y() + node_size)
                && (query.y() + m_radius < octree_origin.y() + octree_size))
            yloc.push_back(getYRightCode(query_node));

        if((query.z() - m_radius  < node_origin.z())
                &&(query.z() - m_radius >octree_origin.z()))
            zloc.push_back(getZLeftCode(query_node));
        if((query.z() + m_radius > node_origin.z() +node_size)
                && (query.z() + m_radius <octree_origin.z() + octree_size))
            zloc.push_back(getZRightCode(query_node));

        //look inside neighboring node
        std::list<unsigned int>::iterator xi,yi,zi;
        xi = xloc.begin();
        while(xi != xloc.end())
        {
            yi = yloc.begin();
            while(yi != yloc.end())
            {
                zi = zloc.begin();
                while(zi != zloc.end())
                {
                    TOctreeNode<T> *node=m_octree->getRoot();
                    traverseToLevel(&node, *xi, *yi, *zi, m_activeDepth);
                    if((node!=NULL)&&(node->getDepth() == m_activeDepth))
                        explore(node, query, neighbors, distances);
                    ++zi;
                }
                ++yi;
            }
            ++xi;
        }
    }
    else
    {
        unsigned int s=query_node->getDepth();
        std::list<unsigned int> xloc, yloc, zloc;
        xloc.push_back(query_node->getXLoc());
        yloc.push_back(query_node->getYLoc());
        zloc.push_back(query_node->getZLoc());
        if((query.x() - m_radius  < node_origin.x())
                &&(query.x() - m_radius > octree_origin.x()))
            xloc.push_back(getXLeftCode(query_node));
        if((query.x() + m_radius > node_origin.x() + node_size)
                && (query.x() + m_radius < octree_origin.x() + octree_size))
            xloc.push_back(getXRightCode(query_node));

        if((query.y() - m_radius  < node_origin.y())
                && (query.y() - m_radius >octree_origin.y() ))
            yloc.push_back(getYLeftCode(query_node));
        if((query.y() + m_radius > node_origin.y() + node_size)
                && (query.y() + m_radius < octree_origin.y() + octree_size))
            yloc.push_back(getYRightCode(query_node));

        if((query.z() - m_radius  < node_origin.z() )
                && (query.z() - m_radius >octree_origin.z()))
            zloc.push_back(getZLeftCode(query_node));
        if((query.z() + m_radius > node_origin.z() +node_size)
                && (query.z() + m_radius <octree_origin.z() + octree_size))
            zloc.push_back(getZRightCode(query_node));

        //look inside neighboring nodes
        std::list<unsigned int>::iterator xi,yi,zi;
        xi = xloc.begin();
        while(xi != xloc.end())
        {
            yi = yloc.begin();
            while(yi != yloc.end())
            {
                zi = zloc.begin();
                while(zi != zloc.end())
                {
                    TOctreeNode<T> *node=m_octree->getRoot();
                    traverseToLevel(&node,*xi,*yi,*zi,s);
                    if( (node != NULL )&&( node->getDepth() == s))
                        explore(node,query,neighbors, distances);
                    ++zi;
                }
                ++yi;
            }
            ++xi;
        }
    }
    return (int)neighbors.size();
}


template<class T>
void TOctreeIterator<T>::traverseToLevel(TOctreeNode<T>** node,
        unsigned int xLocCode, unsigned int yLocCode,
        unsigned int zLocCode, unsigned int k) const
{
    int l=(*node)->getDepth()-1;

    while((*node)->getDepth()>k)
    {
        unsigned int childBranchBit=1<<l;
        unsigned int childIndex=((((xLocCode) & (childBranchBit)) >> l) << 2)
            + (((yLocCode & childBranchBit) >> l) << 1)
            + ( (zLocCode & childBranchBit) >> l);

        if((*node)->getChild(childIndex) != NULL)
        {
            *node=(*node)->getChild(childIndex);
            l--;
        }
        else
            break;
    }
}



template<class T>
void TOctreeIterator<T>::explore(TOctreeNode<T>* node,
        const Vector3d& query_point,
        Neighbor_star_list &neighbors) const
{
    if(node->getDepth() != 0)
    {
        for(unsigned int i=0;i<8;i++)
            if(node->getChild(i) != NULL)
                explore(node->getChild(i), query_point, neighbors);

    }
    else if(node->getNpts(m_setIndex) != 0)
    {
        typename std::deque<T>::iterator iter;
        for(iter = node->points_begin(m_setIndex);
                iter != node->points_end(m_setIndex); ++iter)
        {
           double dist = (query_point - *iter).squaredNorm();
            if(dist < m_sqradius)
                neighbors.push_back(&(*iter));
        }
    }
}


template<class T>
void TOctreeIterator<T>::explore(TOctreeNode<T>* node,
        const Vector3d& query_point,
        Neighbor_star_list &neighbors,
        Distance_list &distances) const
{
    if(node->getDepth() != 0)
    {
        for(unsigned int i = 0; i < 8; ++i)
            if(node->getChild(i) != NULL)
                explore(node->getChild(i), query_point, neighbors, distances);
    }
    else if(node->getNpts(m_setIndex) != 0)
    {
        typename std::deque<T>::iterator iter;
        for(iter = node->points_begin(m_setIndex);
                iter != node->points_end(m_setIndex);
                ++iter)
        {
            double dist = ( query_point - *iter).squaredNorm();
            if(dist < m_sqradius)
            {
                neighbors.push_back(&(*iter));
                distances.push_back(dist);
            }
        }
    }
}

template<class T>
unsigned int TOctreeIterator<T>::getSortedNeighbors(const Vector3d &query,
        Neighbor_star_map &neighbors) const
{
    TOctreeNode<T> *node = locateVector3dNode(query);
    unsigned int n = getSortedNeighbors(query, node, neighbors);
    return n;
}

template<class T>
unsigned int TOctreeIterator<T>::getSortedNeighbors(const Vector3d& query,
        TOctreeNode<T>* query_node,
        Neighbor_star_map &neighbors) const
{
    Vector3d octree_origin = m_octree->getOrigin();
    Vector3d node_origin = query_node->getOrigin();
    double node_size = query_node->getSize();
    double octree_size = m_octree->getSize();

    if(query_node->getDepth() ==  m_activeDepth)
    {
        //find neighboring nodes
        std::list<unsigned int> xloc, yloc, zloc;
        xloc.push_back(query_node->getXLoc());
        yloc.push_back(query_node->getYLoc());
        zloc.push_back(query_node->getZLoc());

        if((query.x() - m_radius  < node_origin.x())
                &&( query.x() - m_radius > octree_origin.x()))
            xloc.push_back(getXLeftCode(query_node));
        if((query.x() + m_radius > node_origin.x() + node_size)
                && (query.x() + m_radius <octree_origin.x() + octree_size))
            xloc.push_back(getXRightCode(query_node));

        if((query.y() - m_radius < node_origin.y())
                &&(query.y() - m_radius >octree_origin.y()))
            yloc.push_back(getYLeftCode(query_node));
        if((query.y() + m_radius > node_origin.y() + node_size)
                && (query.y() + m_radius < octree_origin.y() + octree_size))
            yloc.push_back(getYRightCode(query_node));

        if((query.z() - m_radius  < node_origin.z())
                &&(query.z() - m_radius >octree_origin.z()))
            zloc.push_back(getZLeftCode(query_node));
        if((query.z() + m_radius > node_origin.z() +node_size)
                && (query.z() + m_radius <octree_origin.z() + octree_size))
            zloc.push_back(getZRightCode(query_node));

        //look inside neighboring node
        std::list<unsigned int>::iterator xi,yi,zi;
        xi = xloc.begin();
        while(xi != xloc.end())
        {
            yi = yloc.begin();
            while(yi != yloc.end())
            {
                zi = zloc.begin();
                while(zi != zloc.end())
                {
                    TOctreeNode<T> *node=m_octree->getRoot();
                    traverseToLevel(&node, *xi, *yi, *zi, m_activeDepth);
                    if((node!=NULL)&&(node->getDepth() == m_activeDepth))
                        exploreSort(node, query, neighbors);
                    ++zi;
                }
                ++yi;
            }
            ++xi;
        }
    }
    else
    {
        unsigned int s=query_node->getDepth();
        std::list<unsigned int> xloc, yloc, zloc;
        xloc.push_back(query_node->getXLoc());
        yloc.push_back(query_node->getYLoc());
        zloc.push_back(query_node->getZLoc());
        if((query.x() - m_radius  < node_origin.x())
                &&(query.x() - m_radius > octree_origin.x()))
            xloc.push_back(getXLeftCode(query_node));
        if((query.x() + m_radius > node_origin.x() + node_size)
                && (query.x() + m_radius < octree_origin.x() + octree_size))
            xloc.push_back(getXRightCode(query_node));

        if((query.y() - m_radius  < node_origin.y())
                &&(query.y() - m_radius >octree_origin.y() ))
            yloc.push_back(getYLeftCode(query_node));
        if((query.y() + m_radius > node_origin.y() + node_size)
                && (query.y() + m_radius < octree_origin.y() + octree_size))
            yloc.push_back(getYRightCode(query_node));

        if((query.z() - m_radius  < node_origin.z() )
                &&(query.z() - m_radius >octree_origin.z()))
            zloc.push_back(getZLeftCode(query_node));
        if((query.z() + m_radius > node_origin.z() +node_size)
                && (query.z() + m_radius <octree_origin.z() + octree_size))
            zloc.push_back(getZRightCode(query_node));

        //look inside neighboring nodes
        std::list<unsigned int>::iterator xi,yi,zi;
        xi = xloc.begin();
        while(xi != xloc.end())
        {
            yi = yloc.begin();
            while(yi != yloc.end())
            {
                zi = zloc.begin();
                while(zi != zloc.end())
                {
                    TOctreeNode<T> *node=m_octree->getRoot();
                    traverseToLevel(&node,*xi,*yi,*zi,s);
                    if( (node != NULL )&&( node->getDepth() == s))
                        exploreSort(node, query, neighbors);
                    ++zi;
                }
                ++yi;
            }
            ++xi;
        }
    }
    return (int)neighbors.size();
}

template<class T>
void TOctreeIterator<T>::exploreSort(TOctreeNode<T>* node,
        const Vector3d& query_point,
        Neighbor_star_map &neighbors) const
{
    if(node->getDepth() != 0)
    {
        for(unsigned int i=0; i<8; ++i)
            if(node->getChild(i) != NULL)
                exploreSort(node->getChild(i), query_point, neighbors);
    }
    else if(node->getNpts(m_setIndex) != 0)
    {
        typename std::deque<T>::iterator iter;

        for(iter = node->points_begin(m_setIndex);
                iter != node->points_end(m_setIndex); ++iter)
        {
            double dist = (query_point - *iter).squaredNorm();
            if(dist < m_sqradius)
                neighbors.insert( std::pair<double, T*>(dist, &(*iter)) );
        }
    }
}



template<class T>
unsigned int TOctreeIterator<T>::getXLeftCode(TOctreeNode<T>* node) const
{
    return node->getXLoc()-0x00000001;
}

template<class T>
unsigned int TOctreeIterator<T>::getXRightCode(TOctreeNode<T>* node) const
{
    return node->getXLoc()+(unsigned int)pow2(node->getDepth());
}

template<class T>
unsigned int TOctreeIterator<T>::getYLeftCode(TOctreeNode<T>* node) const
{
    return node->getYLoc()-0x00000001;
}

template<class T>
unsigned int TOctreeIterator<T>::getYRightCode(TOctreeNode<T>* node) const
{
    return node->getYLoc()+(unsigned int)pow2(node->getDepth());
}

template<class T>
unsigned int TOctreeIterator<T>::getZLeftCode(TOctreeNode<T>* node) const
{
    return node->getZLoc()-0x00000001;
}

template<class T>
unsigned int TOctreeIterator<T>::getZRightCode(TOctreeNode<T>* node) const
{
    return node->getZLoc()+(unsigned int)pow2(node->getDepth());
}

template<class T>
void TOctreeIterator<T>::computeCode(const Vector3d& point,
        unsigned int& codx,
        unsigned int& cody,
        unsigned int& codz) const
{
    double multiplier = 1.0 / (m_octree->getSize()) * m_octree->getBinSize();
    const Vector3d &octree_origin = m_octree->getOrigin();
    codx=(unsigned int)((point.x() - octree_origin.x()) * multiplier);
    cody=(unsigned int)((point.y() - octree_origin.y()) * multiplier);
    codz=(unsigned int)((point.z() - octree_origin.z()) * multiplier);
}

template<class T>
TOctreeNode<T>* TOctreeIterator<T>::locateVector3dNode(const Vector3d& point) const
{
    unsigned int codx,cody,codz;
    computeCode(point, codx, cody, codz);

    TOctreeNode<T> *node = m_octree->getRoot();
    if(!node->isInside(point.x(), point.y(), point.z()))
      std::cout<<"not in root"<<std::endl;
    
    traverseToLevel(&node, codx, cody, codz, m_activeDepth);
    
    //if(!node->isInside(point))
    //{
    //  std::cout<<node->getXLoc()<<"\t"<<codx<<std::endl;
    //  std::cout<<node->getYLoc()<<"\t"<<cody<<std::endl;
    //  std::cout<<node->getZLoc()<<"\t"<<codz<<std::endl;
    //  std::cout<<"not in node "<<std::endl;
    //  std::cout<<point<<std::endl;
    //  std::cout<<std::endl;
    //  std::cout<<node->getOrigin()<<std::endl;
    //  std::cout<<std::endl;
    //  std::cout<<node->getOrigin().array() + node->getSize()<<std::endl;
    //  node = m_octree->getRoot();
    //traverseToLevel(&node, codx, cody, codz, m_activeDepth);
    //  node = node->getParent();
    //  std::cout<<"isinside "<<node->isInside(point)<<std::endl;
    //}
    //
    //assert(node->isInside(point));
    return node;
}

template<class T>
bool TOctreeIterator<T>::containsOnly(const Vector3d& query,
        const Exception_set& exceptions) const
{
    TOctreeNode<T> *node = locateVector3dNode(query);
    return containsOnly(query, node, exceptions);
}

template<class T>
bool TOctreeIterator<T>::containsOnly(const Vector3d& query,
        TOctreeNode< T >* query_node,
        const Exception_set& exceptions) const
{
    Vector3d octree_origin = m_octree->getOrigin();
    Vector3d node_origin = query_node->getOrigin();
    double node_size = query_node->getSize();
    double octree_size = m_octree->getSize();

    if(query_node->getDepth() ==  m_activeDepth)
    {
        //find neighboring nodes
        std::list<unsigned int> xloc, yloc, zloc;
        xloc.push_back(query_node->getXLoc());
        yloc.push_back(query_node->getYLoc());
        zloc.push_back(query_node->getZLoc());

        if((query.x() - m_radius  < node_origin.x())
                &&(query.x() - m_radius > octree_origin.x()))
            xloc.push_back(getXLeftCode(query_node));
        if((query.x() + m_radius > node_origin.x() + node_size)
                && (query.x() + m_radius <octree_origin.x() + octree_size))
            xloc.push_back(getXRightCode(query_node));

        if((query.y() - m_radius < node_origin.y())
                &&(query.y() - m_radius >octree_origin.y()))
            yloc.push_back(getYLeftCode(query_node));
        if((query.y() + m_radius > node_origin.y() + node_size)
                && (query.y() + m_radius < octree_origin.y() + octree_size))
            yloc.push_back(getYRightCode(query_node));

        if((query.z() - m_radius  < node_origin.z())
                &&(query.z() - m_radius >octree_origin.z()))
            zloc.push_back(getZLeftCode(query_node));
        if((query.z() + m_radius > node_origin.z() +node_size)
                && (query.z() + m_radius <octree_origin.z() + octree_size))
            zloc.push_back(getZRightCode(query_node));

        //look inside neighboring node
        std::list<unsigned int>::iterator xi,yi,zi;
        xi = xloc.begin();
        while (xi != xloc.end())
        {
            yi = yloc.begin();
            while (yi != yloc.end())
            {
                zi = zloc.begin();
                while (zi != zloc.end())
                {
                    TOctreeNode<T> *node=m_octree->getRoot();
                    traverseToLevel(&node, *xi, *yi, *zi, m_activeDepth);
                    bool ok = true;
                    if((node!=NULL)&&( node->getDepth() == m_activeDepth))
                    {
                        explore(node, query, exceptions, ok);
                        if(!ok)
                            return false;
                    }
                    ++zi;
                }
                ++yi;
            }
            ++xi;
        }
    }
    else
    {
        unsigned int s=query_node->getDepth();
        std::list<unsigned int> xloc, yloc, zloc;
        xloc.push_back(query_node->getXLoc());
        yloc.push_back(query_node->getYLoc());
        zloc.push_back(query_node->getZLoc());
        if((query.x() - m_radius  < node_origin.x())
                && (query.x() - m_radius > octree_origin.x()))
            xloc.push_back(getXLeftCode(query_node));
        if((query.x() + m_radius > node_origin.x() + node_size)
                && (query.x() + m_radius < octree_origin.x() + octree_size))
            xloc.push_back(getXRightCode(query_node));

        if((query.y() - m_radius  < node_origin.y())
                && (query.y() - m_radius >octree_origin.y() ))
            yloc.push_back(getYLeftCode(query_node));
        if((query.y() + m_radius > node_origin.y() + node_size)
                && (query.y() + m_radius < octree_origin.y() + octree_size))
            yloc.push_back(getYRightCode(query_node));

        if((query.z() - m_radius  < node_origin.z() )
                && (query.z() - m_radius >octree_origin.z()))
            zloc.push_back(getZLeftCode(query_node));
        if((query.z() + m_radius > node_origin.z() +node_size)
                && (query.z() + m_radius <octree_origin.z() + octree_size))
            zloc.push_back(getZRightCode(query_node));


        //look inside neighboring nodes
        std::list<unsigned int>::iterator xi,yi,zi;
        xi = xloc.begin();
        while(xi != xloc.end())
        {
            yi = yloc.begin();
            while(yi != yloc.end())
            {
                zi = zloc.begin();
                while(zi != zloc.end())
                {
                    TOctreeNode<T> *node=m_octree->getRoot();
                    traverseToLevel(&node,*xi,*yi,*zi,s);
                    bool ok = true;
                    if( (node != NULL )&&( node->getDepth() == s))
                        explore(node, query, exceptions, ok);
                    if(!ok)
                        return false;
                    ++zi;
                }
                ++yi;
            }
            ++xi;
        }
    }
    return true;
}



template<class T>
void TOctreeIterator<T>::explore(TOctreeNode<T> *node,
        const Vector3d &query_point,
        const Exception_set &exceptions,
        bool &check) const
{
    if(!check)
        return;

    if( node->getDepth() != 0)
    {
        unsigned int i = 0;
        while( (i<8) && (check))
        {
            if(node->getChild(i) != NULL)
                explore(node->getChild(i), query_point, exceptions, check);
            i++;
        }
    }
    else if(node->getNpts(m_setIndex) != 0)
    {
        typename std::deque<T>::iterator iter;
        for(iter = node->points_begin(m_setIndex);
                iter != node->points_end(m_setIndex); ++iter)
        {
            double sqdist = (query_point - *iter).squaredNorm();
            if((sqdist < m_sqradius)
                    && (exceptions.find(&(*iter)) == exceptions.end()))
            {
                check = false;
                return;
            }
        }
    }
}


#endif
