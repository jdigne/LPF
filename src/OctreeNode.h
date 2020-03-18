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
 * @file OctreeNode.h
 * @brief defines an octree node
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

#ifndef OCTREENODE_H
#define OCTREENODE_H

#include <cstdlib>
#include <deque>
#include <iostream>
#include <fstream>
#include <cassert>
#include <Eigen/Dense>

using Eigen::Vector3d;

/**
 * @class TOctreeNode
 * @brief Implements a generic node for a generic octree
 *
 * Templated class implementing a node of the octree, leaf nodes
 * contain the input points.
 */
template<class T>
class TOctreeNode
{
    protected :

        /** @brief pointer to the parent of the node
        */
        TOctreeNode<T> *m_parent;

        /** @brief pointer to the eight children of the node
        */
        TOctreeNode<T> *m_child[8];

        /** @brief
         * child number of the node (depends on the relative location
         * of the node to the middle of its parent)
         \verbatim
           0-------4
          /|      /|
         2-------6 |
         | 1-----|-5
         |/      |/
         3-------7
axis:
x: along direction 0->4
y: along direction 0->2
z: along direction 0->1
         \endverbatim
         * REMARK: this convention is not important because
         * we are dealing with a cube. It will not affect neither
         * representation nor computations
         */
        unsigned int m_nchild;

        /** @brief number of points included in the node or in the
         * node's children
         */
        unsigned int m_npts;

        /** @brief origin of the node*/
        Vector3d m_origin;

        /** @brief level of the node*/
        unsigned int m_depth;

        /**
         * @brief x locational code
         */
        unsigned int m_xloc;

        /**
         * @brief y locational code
         */
        unsigned int m_yloc;

        /**
         * @brief z locational code
         */
        unsigned int m_zloc;

        /** @brief size of the node side*/
        double m_size;

        /** @brief list of points contained in the node
         * (empty if the node is not a leaf)
         * the first list are the input points contained in this node
         * the second and third lists are the filtered sets
         * (2 lists are necessary for buffering)
         */
        std::deque<T> m_points[3];

    public :
        /**
         *  @brief Default constructor initializes all variables
         */
        TOctreeNode();

        /**
         * @brief constructor initializes size, depth and origin
         * @param size;
         * @param origin
         * @param depth
         */
        TOctreeNode(Vector3d & origin, double size, unsigned int depth);

        /**
         * @brief  Destructor
         */
        ~TOctreeNode();

        /**
         * @brief set node size
         * @param size desired node size
         */
        void setSize (double size);

        /**
         * @brief get node size
         * @return size of the node
         */
        double getSize() const;

        /**@brief returns the number of points contained in the node
         * @param index index of the set to process
         * @return number of points contained in the node (if leaf),
         * 0 otherwise
         */
        unsigned int getNpts(unsigned int index) const;

        /**
         * @brief set child number (depends on the relative location of
         * the node to the middle of its parent)
         * @param a child number
         */
        void setNchild(unsigned int a);

        /** @brief get the child number
         * @return child number
         */
        const unsigned int& getNChild() const;

        /**
         * @brief set the origin of the node by a point structure
         * @param pt Origin point
         */
        void setOrigin (Vector3d &pt);

        /** @brief get the origin of a given node
         * @return origin
         */
        Vector3d getOrigin() const;

        /**
         * @brief set parent of a node
         * @param node pointer to the parent node
         */
        void setParent ( TOctreeNode<T> *node );

        /**
         * @brief get parent of a node
         * @return TOctreeNode* pointer to the parent node
         */
        TOctreeNode<T>* getParent() const;

        /** @brief get child of a node
         * @param index of the child
         * @return child node
         */
        TOctreeNode<T>* getChild(unsigned int index);

        /**
         * @brief set level of a node
         * @param l level
         */
        void setDepth ( unsigned int l );

        /** @brief get depth of the node
         * @return depth
         */
        unsigned int getDepth() const;

        /**
         * @brief check if a point given by its coordinates is inside a node
         * @param x coordinates of the point
         * @param y
         * @param z
         * @return true if the point is inside
         */
        bool isInside(double x,double y,double z) const;

        /**
         * @brief check if a point given by its coordinates is inside a node
         *@param p  point
         *@return true if the point is inside
         */
        bool isInside(const Vector3d &p) const;

        /**
         * @brief check if a point given by its coordinates is inside or
         * in a band around the node
         * @param p  point
         * @param d width of the band around the node
         * @return true if the point is inside
         */
        bool isInside(const Vector3d &p, double d) const;

        /**
         * @brief locationnal code method
         * @return unsigned int x locationnal code of the node
         */
        unsigned int getXLoc() const;

        /**
         * @brief locationnal code method
         * @return unsigned int y locationnal code of the node
         */
        unsigned int getYLoc() const;

        /**
         * @brief locationnal code method
         * @return unsigned int z locationnal code of the node
         */
        unsigned int getZLoc() const;

        /**
         * @brief locationnal code method
         * @param Xloc locationnal code of the node
         */
        void setXLoc(unsigned int Xloc);

        /**
         *  @brief locationnal code method
         * @param Yloc locationnal code of the node
         */
        void setYLoc(unsigned int Yloc);

        /**
         *  @brief locationnal code method
         * @param Zloc locationnal code of the node
         */
        void setZLoc(unsigned int Zloc);

        /** @brief get a pointer to the list of points
         * @param index index of the point list
         * @return pointer to the beginning of the list
         */
        typename std::deque<T>::iterator points_begin(unsigned int index);

        /** @brief get a pointer to the end of the list of points
         * @param index index of the point list
         * @return pointer to the end of 'points'
         */
        typename std::deque<T>::iterator points_end(unsigned int index);

        /** @brief get a const pointer to the list of points
         * @param index index of the point list
         * @return const pointer to the beginning of the list
         */
        typename std::deque<T>::const_iterator points_begin(
                unsigned int index) const;

        /** @brief get a const pointer to the end of the list of points
         * @param index index of the point list
         * @return const pointer to the end of 'points'
         */
        typename std::deque<T>::const_iterator points_end(
                unsigned int index) const;

        /** @brief add a point to the list of points included in the cell
         * PREREQUISITE: the node is a leaf in the octree
         * @param pt point to add
         */
        void addInitialPoint(T &pt);

        /** @brief add a point to the list of points included in the cell
         * PREREQUISITE: the node is a leaf in the octree
         * @param pt point to add
         * @param index index of the point list
         */
        void addPoint(T &pt,unsigned int index);

        /** @brief build the i^th child of the node
         * @param index child index
         * @param origin origin of the node
         * @return pointer to the created node
         */
        TOctreeNode<T>* initializeChild(unsigned int index, Vector3d origin);

        /** @brief clear point sets in all the children of the node
         * (and itself)
         * @param index index of the set to clear
         */
        void clearSet(unsigned int index);
};


    template<class T>
TOctreeNode<T>::TOctreeNode()
{
    for(int i = 0 ; i <8 ; i++)
        m_child[i] = NULL;
    m_parent = NULL;
    m_xloc = m_yloc = m_zloc =0;
    m_depth = 0;
    m_npts = 0;
    m_origin = Vector3d();
    m_size = 0.0;
}

    template<class T>
TOctreeNode<T>::TOctreeNode(Vector3d& origin, double size, unsigned int depth)
{
    for(int i = 0 ; i <8 ; i++)
        m_child[i] = NULL;
    m_parent = NULL;
    m_xloc = m_yloc = m_zloc =0;
    m_depth = depth;
    m_npts = 0;
    m_origin = origin;
    m_size = size;
}

    template<class T>
TOctreeNode<T>::~TOctreeNode()
{
    for(int i = 0; i < 3 ; ++i)
        m_points[i].clear();

    m_xloc = m_yloc = m_zloc =0;
    m_depth = 0;
    m_npts = 0;
    for(int i = 0; i<8 ; i++)
        delete m_child[i];
    m_parent = NULL;
    m_origin = Vector3d();
    m_size = 0.0;
}


template<class T>
unsigned int TOctreeNode<T>::getDepth() const
{
    return m_depth;
}

    template<class T>
void TOctreeNode<T>::setDepth(unsigned int l)
{
    m_depth = l;
}

template<class T>
double TOctreeNode<T>::getSize() const
{
    return m_size;
}

    template<class T>
void TOctreeNode<T>::setSize(double size)
{
    m_size = size;
}


template<class T>
unsigned int TOctreeNode<T>::getNpts(unsigned int index) const
{
    return m_points[index].size();
}


template<class T>
const unsigned int& TOctreeNode<T>::getNChild() const
{
    return m_nchild;
}

    template<class T>
void TOctreeNode<T>::setNchild(unsigned int a)
{
    m_nchild = a;
}


    template<class T>
void TOctreeNode<T>::setParent(TOctreeNode *parent)
{
    m_parent = parent;
}

template<class T>
TOctreeNode<T>* TOctreeNode<T>::getParent() const
{
    return m_parent;
}

    template<class T>
TOctreeNode<T>* TOctreeNode<T>::getChild(unsigned int index)
{
    unsigned int i = index % 8;
    return m_child[i];
}


template<class T>
unsigned int TOctreeNode<T>::getXLoc() const
{
    return m_xloc;
}

    template<class T>
void TOctreeNode<T>::setXLoc(unsigned int xloc)
{
    m_xloc = xloc;
}

template<class T>
unsigned int TOctreeNode<T>::getYLoc() const
{
    return m_yloc;
}

    template<class T>
void TOctreeNode<T>::setYLoc(unsigned int yloc)
{
    m_yloc = yloc;
}

template<class T>
unsigned int TOctreeNode<T>::getZLoc() const
{
    return m_zloc;
}

    template<class T>
void TOctreeNode<T>::setZLoc(unsigned int zloc)
{
    m_zloc = zloc;
}

template<class T>
bool TOctreeNode<T>::isInside(double x, double y, double z) const
{
    if ( ( x >= m_origin.x()-1e-16 )&&( x <= m_origin.x() + m_size+1e-16)
            && ( y >= m_origin.y()-1e-16 )&&( y <= m_origin.y() + m_size+1e-16)
            && ( z >= m_origin.z()-1e-16 )&&( z <= m_origin.z() + m_size+1e-16))
        return true;
    else
        return false;
}


template<class T>
bool TOctreeNode<T>::isInside(const Vector3d &p) const
{
    if ( ( p.x() >= m_origin.x()-1e-16 )&&( p.x() <= m_origin.x() + m_size + 1e-16)
            && ( p.y() >= m_origin.y() -1e-16)&&( p.y() <= m_origin.y() + m_size + 1e-16)
            && ( p.z() >= m_origin.z() -1e-16)&&( p.z() <= m_origin.z() + m_size + 1e-16))
        return true;
    else
        return false;
}

template<class T>
bool TOctreeNode<T>::isInside(const Vector3d &p, double d) const
{
    double offset = m_size + d;
    if ( ( p.x() >= m_origin.x() - d )&&( p.x() <= m_origin.x() + offset)
            && ( p.y() >= m_origin.y() - d )&&( p.y() <= m_origin.y() + offset)
            && ( p.z() >= m_origin.z() - d )&&( p.z() <= m_origin.z() + offset))
        return true;
    else
        return false;
}

    template<class T>
void TOctreeNode<T>::setOrigin(Vector3d& pt)
{
    m_origin = pt;
}


template<class T>
Vector3d TOctreeNode<T>::getOrigin() const
{
    return m_origin;
}

template<class T>
typename std::deque<T>::iterator TOctreeNode<T>::points_begin(unsigned int
        index)
{
    return m_points[index].begin();
}

    template<class T>
typename std::deque<T>::iterator TOctreeNode<T>::points_end(unsigned int index)
{
    return m_points[index].end();
}

template<class T>
typename std::deque<T>::const_iterator TOctreeNode<T>::points_begin(
        unsigned int index) const
{
    return m_points[index].begin();
}

template<class T>
typename std::deque<T>::const_iterator TOctreeNode<T>::points_end(
        unsigned int index) const
{
    return m_points[index].end();
}

    template<class T>
void TOctreeNode<T>::addInitialPoint(T &t)
{
    m_points[0].push_back(t);
    m_npts++;
}

    template<class T>
void TOctreeNode<T>::addPoint(T &t, unsigned int index)
{
    m_points[index].push_back(t);
}

    template<class T>
TOctreeNode< T >* TOctreeNode<T>::initializeChild(unsigned int index,
        Vector3d origin)
{
    double size = m_size/2.0;
    unsigned int depth= m_depth -1;
    m_child[index] = new TOctreeNode<T>(origin, size, depth);
    m_child[index]->setParent(this);
    m_child[index]->setNchild(index);

    return m_child[index];
}


    template<class T>
void TOctreeNode<T>::clearSet(unsigned int index)
{
    if(getDepth() == 0)
    {
        m_points[index].clear();
    }
    else
    {
        for(unsigned int i=0; i<8; ++i)
            if(getChild(i) != NULL)
                getChild(i)->clearSet(index);
    }
}


#endif
