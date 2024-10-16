#pragma once

/*
  Copyright (C) 2016 Finnish Meteorological Institute
  Copyright (C) 2016-2024 CSC -IT Center for Science

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
#include "tools.hpp"

#include <array>
#include <cstdint>
#include <mpi.h>

namespace fsgrid {
struct CartesianGrid {
private:
   using FsSize_t = FsGridTools::FsSize_t;
   using FsIndex_t = FsGridTools::FsIndex_t;
   using LocalID = FsGridTools::LocalID;
   using GlobalID = FsGridTools::GlobalID;
   using Task_t = FsGridTools::Task_t;

public:
   //!< Global size of the simulation space, in cells
   std::array<FsSize_t, 3> globalSize = {};
   //!< Information about whether a given direction is periodic
   std::array<bool, 3> periodic = {};
   //!< Rank in the parent communicator
   int32_t parentRank = -1;
   //!< Size of the parent communicator
   int32_t parentSize = 0;
   //!< The number of fs ranks
   int32_t size = 0;
   //!< Colour that identifies fs ranks
   int32_t colourFs = MPI_UNDEFINED;
   //!< Colour that identifies non-fs ranks
   int32_t colourAux = MPI_UNDEFINED;
   //!< Number of tasks in each direction
   std::array<Task_t, 3> numTasksPerDim = {};
   //! MPI Cartesian communicator used in this grid
   MPI_Comm cartesian3d = MPI_COMM_NULL;
   //!< This task's rank in the communicator
   int32_t rank = -1;
   //!< This task's position in the 3d task grid
   std::array<Task_t, 3> taskPosition = {};
   //!< Local size of simulation space handled by this task (without ghost cells)
   std::array<FsIndex_t, 3> localSize = {};
   //!< Offset of the local coordinate system against the global one
   std::array<FsIndex_t, 3> localStart = {};
   //!< Local size of simulation space handled by this task (including ghost cells)
   std::array<FsIndex_t, 3> storageSize = {};

   CartesianGrid(std::array<FsSize_t, 3> globalSize, MPI_Comm parentComm, std::array<bool, 3> periodic,
                 const std::array<Task_t, 3>& decomposition, int32_t stencilSize)
       : globalSize(globalSize), periodic(periodic), parentRank(FsGridTools::getCommRank(parentComm)),
         parentSize(FsGridTools::getCommSize(parentComm)), size(FsGridTools::getNumFsGridProcs(parentSize)),
         colourFs(FsGridTools::computeColourFs(parentRank, size)),
         colourAux(FsGridTools::computeColourAux(parentRank, parentSize, size)),
         numTasksPerDim(FsGridTools::computeNumTasksPerDim(globalSize, decomposition, size, stencilSize)),
         cartesian3d(FsGridTools::createCartesianCommunicator(parentComm, colourFs, colourAux, parentRank,
                                                              numTasksPerDim, periodic)),
         rank(FsGridTools::getCartesianRank(colourFs, cartesian3d)),
         taskPosition(FsGridTools::getTaskPosition(cartesian3d)),
         localSize(FsGridTools::calculateLocalSize(globalSize, numTasksPerDim, taskPosition, rank, stencilSize)),
         localStart(FsGridTools::calculateLocalStart(globalSize, numTasksPerDim, taskPosition)),
         storageSize(FsGridTools::calculateStorageSize(globalSize, localSize, stencilSize)) {}

   CartesianGrid() {}

   void finalize() {
      if (cartesian3d != MPI_COMM_NULL) {
         FSGRID_MPI_CHECK(MPI_Comm_free(&cartesian3d), "Failed to free MPI comm");
      }
   }
};
}
