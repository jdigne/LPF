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
 * @file Octree.h
 * @brief declares an octree storing structure
 * @author Julie Digne
 * @date 2012/10/10
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

#ifndef OCTREE_H
#define OCTREE_H

#include "OctreeNode.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "utils_octree.h"

using Eigen::Vector3d;

/**
 * @class TOctree
 * @brief Data structure for storing and sorting the points
 *
 * Templated octree data structure permitting to sort the input points
 */
template<class T>
class TOctree
{
    public :
        typedef typename std::vector< std::vector<TOctreeNode<T>* > >
        OctreeNode_collection;

    public : //constructors/destructors
        /**
         * @brief default constructor
         */
        TOctree();

        /** @brief initialize an octree with given parameters
         * @param depth depth of the octree
         */
        TOctree(unsigned int depth);

        /** @brief initialize an octree with given parameter
         * @param origin origin of the octree
         * @param size  size of the loose bounding box
         * @param depth depth of the octree
         */
        TOctree(const Vector3d &origin, double size, unsigned int depth);

        /** @brief Destructor
        */
        ~TOctree();

    public : //accessors and modifiers

        /** @brief get octree depth
         * @return depth
         */
        unsigned int  getDepth() const;

        /** @brief set octree depth
         * @param depth
         */
        void setDepth(unsigned int depth);

        /** @brief get origign
         *  @return origin of the octree
         */
        const Vector3d& getOrigin() const;

        /** @brief get number of points
         * @return number of points
         */
        unsigned int getNpoints() const;

        /** @brief get side size of the octree
         * @return size of the octree
         */
        double getSize() const;

        /** @brief set size of the octree
         * @param size
         */
        void setSize(double size);

        /** @brief returns the size of the smallest cell
         * @return side length of the smallest cell
         */
        double getSmallestCellSize() const;

        /** @brief get the bin size
         * @return binsize
         */
        unsigned int getBinSize() const;

        /** @brief get root of the octree
         * @return root of the octree
         */
        TOctreeNode<T>* getRoot() const;

    public : //adding points

        /** @brief initialize the octree with the origin and size
         * @param origin origin of the octree
         * @param size side size
         **/
        void initialize(const Vector3d & origin, double size, unsigned int depth);

        /** @brief Adding an initial point to the octree
         * @param pt point to add
         */
        void addInitialPoint(T &pt);

        /** @brief Adding a point to the octree
         * @param pt point to add
         * @param index index of the list to add the point to (0,1,2)
         */
        void addPoint(T &pt, unsigned int index);

        /**
         * @brief Adding a batch of points to the octree
         * @param begin begin iterator of the batch
         * @param end end iterator of the batch
         * @return number of added points
         */
        template<class Iterator>
            unsigned int addInitialPoints(Iterator begin, Iterator end);

        /** @brief print the mean number of points per non empty cell
         * at each level
         */
        void printOctreeStat();

        /** @brief get all nodes at given depth
         * @param depth input depth
         * @param starting_node node to start from
         * @param[out] nodes vector of nodes
         **/
        void getNodes(unsigned int depth, TOctreeNode<T>* starting_node,
                std::vector< TOctreeNode<T>* > &nodes);

        /** @brief get nodes in separate buckets according to the parity
         * of the indices with respect to the three coordinates (8 buckets)
         * process each node separately
         * @param depth input depth
         * @param starting_node node to start from
         * @param[out] node_collection vector of vector of nodes (each 
         *   vector of nodes can be processed separately)
         */
        void getNodes(unsigned int depth, TOctreeNode<T>* starting_node, 
                std::vector< std::vector<TOctreeNode<T>* > > &node_collection); 

        /** @brief clear set i in all octree nodes
         * @param index index of the sets to clear
         */
        void clearSet(unsigned int index);

    protected :

        /** @brief Maximum depth of the octree*/
        unsigned int m_depth;

        /** @brief Number of points sorted in the octree*/
        unsigned int m_npoints;

        /** @brief Origin of the octree*/
        Vector3d m_origin;

        /** @brief size of the side of the octree*/
        double m_size;

        /** @brief BinSize
         * the bin size is an octree parameter that determines the locational 
         *code of each node
         */
        unsigned int m_binsize;

        /** @brief root of the octree*/
        TOctreeNode<T> *m_root;

        /** @brief number of non-empty cells per level*/
        std::vector<unsigned int> m_nb_non_empty_cells;
};

    template<class T>
TOctree<T>::TOctree()
{
    m_size = 0;
    m_depth = 0;
    m_binsize = 0;
    m_npoints = 0;
    m_origin = Vector3d();
    m_root = NULL;
}


    template<class T>
TOctree<T>::TOctree(unsigned int depth)
{
    m_size = 0;
    m_depth = depth;
    m_binsize = pow2(depth);
    m_npoints = 0;
    m_root = NULL;
    m_nb_non_empty_cells.assign(depth,0);
}


    template<class T>
TOctree<T>::TOctree(const Vector3d& origin, double size, unsigned int depth)
{
    m_size = size;
    m_depth = depth;
    m_binsize = pow2(depth);
    m_origin = origin;
    m_npoints = 0;
    m_nb_non_empty_cells.assign(depth,0);
    
    m_root = new TOctreeNode<T>(m_origin, m_size, m_depth);
    m_root->setXLoc(0);
    m_root->setYLoc(0);
    m_root->setZLoc(0);
    m_root->setParent(NULL);
}

    template<class T>
TOctree<T>::~TOctree()
{
    m_size = 0;
    m_depth = 0;
    m_binsize = 0;
    m_npoints = 0;
    m_origin = Vector3d();

    if(m_root != NULL)
    {
        delete m_root;
        m_root = NULL;
    }
    m_nb_non_empty_cells.clear();
}


    template<class T>
void TOctree<T>::initialize(const Vector3d& origin, double size, unsigned int depth)
{
    m_size = size;
    m_origin = origin;
    m_depth = depth;

    m_root = new TOctreeNode<T>(m_origin, m_size, m_depth);
    m_root->setXLoc(0);
    m_root->setYLoc(0);
    m_root->setZLoc(0);
    m_root->setParent(NULL);
}


template<class T>
unsigned int TOctree<T>::getDepth() const
{
    return m_depth;
}

    template<class T>
void TOctree<T>::setDepth(unsigned int depth)
{
    m_depth = depth;
    m_binsize = pow2(depth);
    m_nb_non_empty_cells.clear();
    m_nb_non_empty_cells.assign(depth,0);
}

template<class T>
unsigned int TOctree<T>::getNpoints() const
{
    return m_npoints;
}


template<class T>
double TOctree<T>::getSize() const
{
    return m_size;
}

template<class T>
void TOctree<T>::setSize(double size)
{
    m_size = size;
}

template<class T>
double TOctree<T>::getSmallestCellSize() const
{
    return m_size/((double)pow2(m_depth-1));
}

template<class T>
unsigned int TOctree<T>::getBinSize()
    const
{
    return m_binsize;
}

template<class T>
const Vector3d& TOctree<T>::getOrigin() const
{
    return m_origin;
}

template<class T>
TOctreeNode<T>* TOctree<T>::getRoot() const
{
    return m_root;
}

template<class T>
    template<class Iterator>
unsigned int TOctree<T>::addInitialPoints(Iterator begin, Iterator end)
{
    Iterator it = begin;
    while(it != end)
    {
        T &a = *it;
        assert(m_root->isInside(a.x(), a.y(), a.z()));
        addInitialPoint(a);
        ++it;
    }
    return m_npoints;
}

template<class T>
void TOctree<T>::addInitialPoint(T& pt)
{
    unsigned int codx=(unsigned int)((pt.x() - m_origin.x()) / m_size *
            m_binsize);
    unsigned int cody=(unsigned int)((pt.y() - m_origin.y()) / m_size *
            m_binsize);
    unsigned int codz=(unsigned int)((pt.z() - m_origin.z())/ m_size *
            m_binsize);
    TOctreeNode<T> *node=getRoot();
    unsigned int l=node->getDepth()-1;

    //traverse the octree until we reach a leaf
    while(node->getDepth() != 0)
    {
        unsigned int childBranchBit=1<<l;
        unsigned int x = ( ( codx & childBranchBit) >> l );
        unsigned int y = ( ( cody & childBranchBit) >> l );
        unsigned int z = ( ( codz & childBranchBit) >> l );
        unsigned int childIndex = (x<<2) + (y<<1) + z;

        if(node->getChild(childIndex) == NULL)
        {
            double childSize = node->getSize()/2.0;
            unsigned int childDepth = node->getDepth() - 1;
            Vector3d origin = node->getOrigin();
            Vector3d childOrigin = Vector3d( origin.x()  + x * childSize,
                    origin.y() + y * childSize,
                    origin.z() + z * childSize);

            TOctreeNode<T> *child = node->initializeChild(childIndex, 
                    childOrigin);

            child->setXLoc( node->getXLoc() + ( x<<(childDepth) ) );
            child->setYLoc( node->getYLoc() + ( y<<(childDepth) ) );
            child->setZLoc( node->getZLoc() + ( z<<(childDepth) ) );
            m_nb_non_empty_cells[childDepth] += 1;
        }
        node = node->getChild(childIndex);
        l--;
    }

    //add the point to the leaf.
    node->addInitialPoint(pt);
    m_npoints++;
}


template<class T>
void TOctree<T>::addPoint(T& pt, unsigned int index)
{
    unsigned int codx=(unsigned int)((pt.x() - m_origin.x()) / m_size *
            m_binsize);
    unsigned int cody=(unsigned int)((pt.y() - m_origin.y()) / m_size *
            m_binsize);
    unsigned int codz=(unsigned int)((pt.z() - m_origin.z())/ m_size *
            m_binsize);
    TOctreeNode<T> *node=getRoot();
    unsigned int l=node->getDepth()-1;

    //traverse the octree until we reach a leaf
    while(node->getDepth() != 0)
    {
        unsigned int childBranchBit=1<<l;
        unsigned int x = ( ( codx & childBranchBit) >> l );
        unsigned int y = ( ( cody & childBranchBit) >> l );
        unsigned int z = ( ( codz & childBranchBit) >> l );
        unsigned int childIndex = (x<<2) + (y<<1) + z;

#ifdef _OPENMP
        //critical section to avoid conflict when adding points to new cells
        //during parallel scale space iterations
#pragma omp critical
        {
#endif
            if(node->getChild(childIndex) == NULL)
            {
                double childSize = node->getSize()/2.0;
                unsigned int childDepth = node->getDepth() - 1;
                Vector3d origin = node->getOrigin();
                Vector3d childOrigin = Vector3d( origin.x()  + x * childSize,
                        origin.y() + y * childSize,
                        origin.z() + z * childSize);

                TOctreeNode<T> *child = node->initializeChild(childIndex, 
                        childOrigin);

                child->setXLoc( node->getXLoc() + ( x<<(childDepth) ) );
                child->setYLoc( node->getYLoc() + ( y<<(childDepth) ) );
                child->setZLoc( node->getZLoc() + ( z<<(childDepth) ) );
                m_nb_non_empty_cells[childDepth] += 1;
            }
#ifdef _OPENMP
        }
#endif
        node = node->getChild(childIndex);
        l--;
    }

    //add the point to the leaf.
    node->addPoint(pt,index);
}


template<class T>
void TOctree<T>::getNodes(unsigned int depth, TOctreeNode<T> *starting_node, 
        std::vector< TOctreeNode<T>* >& nodes)
{
    if(starting_node->getDepth() == depth)
        nodes.push_back(starting_node);
    else
    {
        for(int i = 0; i < 8; ++i )
            if(starting_node->getChild(i) != NULL)
                getNodes(depth, starting_node->getChild(i), nodes);
    }
}

    template<class T>
void TOctree<T>::getNodes(unsigned int depth,
        TOctreeNode<T> *starting_node,
        std::vector< std::vector<TOctreeNode<T>* > > &node_collection)
{
    std::vector< TOctreeNode<T>* > nodes;
    getNodes(depth, starting_node, nodes);

    for(int i = 0 ; i < 8; ++i)
    {
        std::vector<TOctreeNode<T>* > ov;
        node_collection.push_back(ov);
    }

    while(! nodes.empty())
    {
        TOctreeNode<T> *node = nodes.back();
        nodes.pop_back();
        unsigned int nchild = node->getNChild();
        node_collection[nchild].push_back(node); 
    }
}

    template<class T>
void TOctree<T>::printOctreeStat()
{
    double size = m_size;
    for(int i = m_depth-1; i >= 0; i--)
    {
        std::cout<<"level "<<i<<" : "<<size<<" ; mean number of points: "
            <<(double)m_npoints/((double)m_nb_non_empty_cells[i])
            <<std::endl;
        size = size / 2.0;
    }
}

    template<class T>
void TOctree<T>::clearSet(unsigned int index)
{
    if(index > 2)
        return;
    TOctreeNode<T> *node = getRoot();
    node->clearSet(index);
}


#endif
