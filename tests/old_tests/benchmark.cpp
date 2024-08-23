#include <iostream>
#include "../fsgrid.hpp"
/*
  Copyright (C) 2016 Finnish Meteorological Institute
  Copyright (C) 2016 CSC -IT Center for Science 

  This file is part of fsgrid
 
  fsgrid is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  fsgrid is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY;
  without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with fsgrid.  If not, see <http://www.gnu.org/licenses/>.
*/

template<class T, int stencil> void timeit(std::array<int32_t, 3> globalSize, std::array<bool, 3> isPeriodic, int iterations){
   double t1,t2;   
   FsGridCouplingInformation gridCoupling;
   FsGrid<T ,stencil> testGrid(globalSize, MPI_COMM_WORLD, isPeriodic, gridCoupling);
   int rank,size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   testGrid.updateGhostCells();         

   MPI_Barrier(MPI_COMM_WORLD);   
   t1=MPI_Wtime();  
   for(int i = 0; i < iterations; i++) {
      testGrid.updateGhostCells();         
   }
   MPI_Barrier(MPI_COMM_WORLD);   
   t2=MPI_Wtime();
   if(rank==0)
      printf("%g s per update: nprocs %d, grid is %d x %d x %d, stencil %d, element size %ld \n", (t2 - t1)/iterations, size, globalSize[0], globalSize[1], globalSize[2], stencil, sizeof(*testGrid.get(0,0,0)));
}


int main(int argc, char** argv) {
   
   MPI_Init(&argc,&argv);

   int rank,size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   // Create a 8×8 Testgrid
   std::array<int32_t, 3> globalSize{2000,1000,1};
   std::array<bool, 3> isPeriodic{false,false,true};

   const int iterations = 200;

   timeit<std::array<double,1>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,2>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,4>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,8>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,16>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,32>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,64>, 2>(globalSize, isPeriodic, iterations);
   timeit<std::array<double,128>, 2>(globalSize, isPeriodic, iterations);
   
      
   MPI_Finalize();
   return 0;
}
