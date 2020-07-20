/**
 * This file is part of lpf
 *
 * @file main.cpp
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
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <ctime>
#include "cmdLine.h"
#include "io.h"
#include "Parameterizer.h"
#include "DescriptorComputation.h"
#include "DescriptorAnalysis.h"
#include "DescriptorReconstruction.h"
#include "utils.h"
#include "Denoising.h"
#include <Eigen/Dense>

#include "LPF.h"
#include "Octree.h"
#include "OctreeIterator.h"
#include "OctreeNode.h"


using namespace std;
using Eigen::Vector3d;

typedef TOctree<Vector3d> Octree;
typedef TOctreeIterator<Vector3d> OctreeIterator;

/** @brief display the usage of the program
 */
void usage()
{
    std::cout<<"USAGE"<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<"descriptor INPUT OUTPUT -r <radius> -p -type"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"INPUT  (mandatory) ascii file containing the point set to analyze"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"SEEDS  (mandatory) ascii file containing the seeds where analysis will be performed"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"OUTPUT (mandatory) output file combined file"<<std::endl;
    std::cout<<"-r     <radius> (recommended) Radius that will be used to interpolate x,y positions. Otherwise $r=1$ will be used"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"-n     <nbins> number of bins for building the generic pattern. Alternatively use a file (-g)"<<std::endl;
    std::cout<<"-l     <lambda> penalty for the sparse optimization"<<std::endl;
    std::cout<<"-s     save immediately the descriptor (low memory consumption)"<<std::endl;
    std::cout<<"-A     <natoms> analyse the descriptors"<<std::endl;
    std::cout<<"-R     <natoms> analyse the descriptors and resample"<<std::endl;
    std::cout<<"-D     <natoms> analyse the descriptors and denoise"<<std::endl;
    std::cout<<"-v     verify the reconstruction from the descriptors and outputs the point cloud"<<std::endl;
    std::cout<<"       if and only if the descriptor can provide such a reconstruction (LPF) "<<std::endl;
    std::cout<<"       if -v and -s are both set, then -v will prevail"<<std::endl;
    std::cout<<"-N     <n> number of iterations for the dictionary computation. Default:10 "<<std::endl;
    std::cout<<"-g     <file> sets the generic pattern to a precise file"<<std::endl;
    std::cout<<"-d     <file> input dictionary file"<<std::endl;
    std::cout<<"-c     Compute parameterization instead of using a random initial seed orientation"<<std::endl;    
    std::cout<<"-p     <file> use a parameterization file for the seeds."<<std::endl;
}

/**
 * @brief main function for the LPF analysis
 * See usage() function for command line options
 */
int main(int argc, char **argv)
{
  CmdLine cmd;
  //handling command line options
  std::stringstream f;
  string infile, seedfile, outfile, genpatfile, seedparamfile, dictionaryfile;
  double radius = -1;
  int nbins = 0;
  int natoms = 0;
  int parallel_flag = -1;
  int lpf_flag = -1;
  int niter = 10;
  double lambda=0.2;

  cmd.add(make_option('r', radius));
  cmd.add(make_option('n', nbins));
  cmd.add(make_switch('s'));
  cmd.add(make_switch('v'));
  cmd.add(make_option('A', natoms));
  cmd.add(make_option('D', natoms));
  cmd.add(make_option('R', natoms));
  cmd.add(make_option('N', niter));
  cmd.add(make_option('g', genpatfile));
  cmd.add(make_option('p', seedparamfile));
  cmd.add(make_option('d', dictionaryfile));
  cmd.add(make_option('l',lambda));
  cmd.add(make_switch('c'));

  try {
      cmd.process(argc, argv);
  } catch(std::string str) {
      std::cerr << "Error: " << str << std::endl<<std::endl;
      usage();
      return EXIT_FAILURE;
  }

  if(argc < 3)
  {
      usage();
      return EXIT_FAILURE;
  }
  
  infile = argv[1];
  seedfile = argv[2];
  outfile = argv[3];

  std::cout<<"Descriptor program launched with parameters: "<<std::endl;
  std::cout << "Input file : " << infile << std::endl;
  std::cout << "Seed file : "  << seedfile << std::endl;
  std::cout << "Output file: " << outfile << std::endl;
  std::cout << "radius : " << radius << std::endl;
  std::cout << "number of bins : " << nbins<< std::endl;
  std::cout << "number of iterations : " << niter<< std::endl;
  std::cout << "number of dictionary atoms : " << natoms<< std::endl;
  
  if(cmd.used('d'))
    std::cout<<"Using input dictionary "<<dictionaryfile<<std::endl;
  
  if(cmd.used('p'))
    std::cout<<"Using input seed parameterization "<<seedparamfile<<std::endl;

  srand((unsigned int) time(0));


  if( (! cmd.used('r')) || (radius < 0) )
  {
    radius = 1.0;
    std::cout<<"No radius given. Using default radius =  "<<radius<<" ."<<std::endl;
  }
  
  if((! cmd.used('g')) && ( (! cmd.used('n')) || (nbins < 0)) )
  {
    nbins = 10;
    std::cout<<"No grid size given and no generic pattern from file. Using default grid size = 10 and a 2d grid (=> 100 points)."<<std::endl;
  }
  if(!cmd.used('l') || lambda < 0)
  {
      lambda = 0.2;
      std::cout << "lambda (sparsity favoring weight) : " << lambda<< std::endl;
  }

  time_t start, end, glob_start, glob_end;


  std::time(&glob_start);
  
  unsigned int npts, nwords;
  if(! IO::getnlines(infile.c_str(), npts, nwords))
    return EXIT_FAILURE;
  
  unsigned int nseeds, nseedwords;
  if(! IO::getnlines(seedfile.c_str(), nseeds, nseedwords))
    return EXIT_FAILURE;
  
  std::vector<Vector3d> points;
  std::vector<Vector3d> seeds;
  std::vector<Vector3d> normals;
  std::vector<Vector3d> t1s;
  std::vector<Vector3d> t2s;
  std::vector<Vector3d> queries;
  
  Vector3d origin;
  unsigned int depth;
  double size;
  
  //reading points and seeds
  std::time(&start);
  IO::readPoints(infile.c_str(), radius, points, origin, size, depth);
  IO::readPoints(seedfile.c_str(), seeds);
  std::time(&end);
  
  std::cout<<points.size()<<" points and "<<seeds.size()<<" seeds read in "<<difftime(end, start)<<" s."<<std::endl;
  
  std::time(&start);
  Octree *octree = new Octree(origin, size, depth);
  octree->addInitialPoints(points.begin(), points.end());
  std::time(&end);
  std::cout<<"Octree built in "<<difftime(end, start)<<" s."<<std::endl;
  
  
  if(! cmd.used('p'))
  {
    //compute parameterization
    std::time(&start);
    Parameterizer *param = new Parameterizer(radius, octree);
    
    if(cmd.used('c'))
      param->compute_tangent_params(seeds, t1s, t2s, normals);
    else
      param->compute_random_params(seeds, t1s, t2s, normals);
    
    delete param;
    std::time(&end);
    
    std::cout<<"Parameterization time "<<difftime(end, start)<<" s."<<std::endl;
    IO::saveParams("params.asc", t1s, t2s, normals);
  }
  else 
  {
    //reed parameterization
    IO::readSeedParams(seedparamfile.c_str(), t1s, t2s, normals);
  }
  
  
  //debug
  check_frames(t1s,t2s,normals);
  
  
  //compute generic pattern
  bool loadgenpat = false;
  if(cmd.used('g'))
  {
    std::cout<<"Reading a generic pattern from a file"<<std::endl;
    loadgenpat = read_qpoints(genpatfile.c_str(), queries);
  }
  if(! cmd.used('g') || ! loadgenpat)
  {
    std::cout<<"Generating a 2D grid"<<std::endl;
    generate_2d_grid(nbins, queries);
  }
  
  //compute descriptors
  DescriptorComputation desc(radius, queries, octree);
  
  
  MatrixXd D, alpha;
  if(cmd.used('d'))
  {
    //initialize the dictionary from the file dictionaryfile
    bool ok = IO::readMatrix(dictionaryfile.c_str(), D);
    
    if(!ok)
    {
      std::cout<<D.cols()<<"\t"<<D.rows()<<std::endl;
      std::cerr<<"File "<<dictionaryfile<<" could not be opened; continuing using random dictionary."<<std::endl;
      D.resize(0,0);
    }
    else if(D.cols() != natoms)
    {
      std::cerr<<"Inconsistent number of atoms or signal size between the dictionary file (" <<D.cols()<<" columns) and parameter (="<<natoms<<")."<<std::endl;
      std::cerr<<"Continuing with random dictionary."<<std::endl;
      D.resize(0,0);
    }
    else
    {
      std::cout<<"Dictionary ("<<D.cols()<<"x"<<D.rows()<<" matrix) successfully read from file."<<std::endl; 
    }
  }
  
  if(lpf_flag)
  {
    std::cout<<"Computing LPF descriptors"<<std::endl;

    if( (cmd.used('s')) && (! cmd.used('v')) )
    {
      std::time(&start);
      desc.computeAndSaveDescriptors<LPF>(seeds, points, t1s, t2s, normals, outfile.c_str());
      std::time(&end);
      std::cout<<"LPF Computation and Saving time "<<difftime(end, start)<<" s."<<std::endl;
    }
    else
    {
      if (cmd.used('A'))
      {
        std::cout<<"Analyzing"<<std::endl;
        std::vector<LPF*> descriptors;
        std::time(&start);
        desc.computeDescriptorsWithEfficiencyOptim(seeds, points, t1s, t2s, normals, descriptors);
        std::time(&end);
        std::cout<<"LPF Computation  time "<<difftime(end, start)<<" s."<<std::endl;
        
        std::cout<<"Perform descriptor analysis"<<std::endl;
        
        DescriptorAnalysis analysis(descriptors);
        analysis.setGenericPattern(queries);
        
        MatrixXd t(descriptors.size(), 3);
        t.setZero();
        unsigned int dic_update_niter = 10;
        
        
        std::cout<<"Learning a dictionary of size "<<natoms<<std::endl;
        std::time(&start);
        analysis.learnLPFdictionaryAndFramesICP(descriptors, D, alpha, t, t1s, t2s, normals, natoms, lambda, dic_update_niter, niter);
        std::time(&end);
        std::cout<<"Dictionary learning time "<<difftime(end, start)<<" s."<<std::endl;
        
        DescriptorReconstruction rec(radius);
        std::cout<<"Translating the seeds"<<std::endl;
        rec.adaptSeeds(seeds, t);
        
        while(! descriptors.empty())
        {
          if(descriptors.back() != NULL)
            delete descriptors.back();
          descriptors.pop_back();
        }
        
        IO::saveMatrix("dictionary.asc", D);
        IO::saveMatrix("coeffs.asc", alpha);
        IO::savePoints("genpat.asc", queries);
        IO::saveParams("params_with_learning.asc", t1s, t2s, normals);
        IO::savePoints("seeds_with_learning.asc", seeds);
      }
      else if(cmd.used('R') || cmd.used('D'))
      {
        std::vector<LPF*> descriptors;
        std::time(&start);
        desc.computeDescriptorsWithEfficiencyOptim(seeds, points, t1s, t2s, normals, descriptors);
        std::time(&end);
        std::cout<<"LPF Computation  time "<<difftime(end, start)<<" s."<<std::endl;
        
        std::cout<<"Perform descriptor analysis"<<std::endl;

        DescriptorAnalysis analysis(descriptors);
        analysis.setGenericPattern(queries);
 
        MatrixXd t(descriptors.size(), 3);
        t.setZero();
        unsigned int dic_update_niter = 10;
        
        std::cout<<"Learning a dictionary of size "<<natoms<<std::endl;
        std::time(&start);
        analysis.learnLPFdictionaryAndFramesICP(descriptors, D, alpha, t, t1s, t2s, normals, natoms, lambda, dic_update_niter, niter);
        std::time(&end);
        std::cout<<"Dictionary learning time "<<difftime(end, start)<<" s."<<std::endl;
        
        std::cout<<"Reconstructing the LPFs"<<std::endl;
        DescriptorReconstruction rec(radius);
        std::vector<LPF*> recdescs;
        rec.buildDescriptorsFromMatrix(D*alpha, recdescs);
        
        std::cout<<"Translating the seeds"<<std::endl;
        rec.adaptSeeds(seeds, t);

        
        if(cmd.used('D'))
        {
          std::cout<<"Denoising"<<std::endl;
          std::time(&start);
          Denoising denoiser(origin, size, depth);
          denoiser.denoise(points, radius, seeds, t1s, t2s, normals, queries, recdescs,0.5);
          std::time(&end);
          std::cout<<"Denoising time "<<difftime(end, start)<<" s."<<std::endl;
          IO::savePoints(outfile.c_str(), points);
          IO::savePoints("seeds.xyz", seeds);
        }
        else if (cmd.used('R'))
        {
          std::cout<<"Resampling"<<std::endl;
          std::vector< Vector3d > recpoints;
          double precision = radius/ std::sqrt(queries.size());
          std::time(&start);
          rec.buildOctree(origin, size, precision);
          rec.reconstructAndMerge(seeds, precision, t1s, t2s, normals, queries, recdescs, recpoints);
          std::time(&end);
          std::cout<<"Resampling time "<<difftime(end, start)<<" s."<<std::endl;
          IO::savePoints(outfile.c_str(), recpoints);
        }
        while(! recdescs.empty())
        {
          if(recdescs.back() != NULL)
            delete recdescs.back();
          recdescs.pop_back();
        }
        std::cout<<"done"<<std::endl;
        while(! descriptors.empty())
        {
          if(descriptors.back() != NULL)
            delete descriptors.back();
          descriptors.pop_back();
        }
      }
      else if(cmd.used('v'))
      {
        std::vector<LPF*> descriptors;
        std::time(&start);
        desc.computeDescriptorsWithEfficiencyOptim(seeds, points, t1s, t2s, normals, descriptors);
        IO::saveParams("params_with_optim.asc", t1s, t2s, normals);
        IO::savePoints("seeds_with_optim.asc", seeds);
        std::time(&end);
        std::cout<<"LPF Computation  time "<<difftime(end, start)<<" s."<<std::endl;
        std::cout<<"reconstructing the points from the LPFs"<<std::endl;
        DescriptorReconstruction rec(radius);
        rec.reconstructAndSave(seeds, t1s, t2s, normals, queries, descriptors, "uplusv.xyz");
        rec.reconstructAndSaveTargetZones(seeds, t1s, t2s, normals, descriptors, "target_zones.xyz");
        
        while(! descriptors.empty())
        {
          if(descriptors.back() != NULL)
            delete descriptors.back();
          descriptors.pop_back();
        }
      }
    }
  }

  std::time(&glob_end);
  

  std::cout<<"Total computation time "<<difftime(glob_end, glob_start)
             <<" s."<<std::endl;


  delete octree;
  return EXIT_SUCCESS;
}
