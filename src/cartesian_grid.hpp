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
#include <numeric>

namespace fsgrid {
struct CartesianGrid {
private:
   using FsSize_t = FsGridTools::FsSize_t;
   using FsIndex_t = FsGridTools::FsIndex_t;
   using LocalID = FsGridTools::LocalID;
   using GlobalID = FsGridTools::GlobalID;
   using Task_t = FsGridTools::Task_t;

   //!< Width of the ghost cells frame
   int32_t stencilSize = 0;
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

public:
   CartesianGrid(std::array<FsSize_t, 3> globalSize, MPI_Comm parentComm, std::array<bool, 3> periodic,
                 const std::array<Task_t, 3>& decomposition, int32_t stencilSize)
       : stencilSize(stencilSize), globalSize(globalSize), periodic(periodic),
         parentRank(FsGridTools::getCommRank(parentComm)), parentSize(FsGridTools::getCommSize(parentComm)),
         size(FsGridTools::getNumFsGridProcs(parentSize)), colourFs(FsGridTools::computeColourFs(parentRank, size)),
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

   /*! Returns the task responsible, and its localID for handling the cell with the given GlobalID
    * \param id GlobalID of the cell for which task is to be determined
    * \return a task for the grid's cartesian communicator
    */
   std::pair<int32_t, LocalID> getTaskForGlobalID(GlobalID id) {
      // Transform globalID to global cell coordinate
      std::array<FsIndex_t, 3> cell = FsGridTools::globalIDtoCellCoord(id, globalSize);

      // Find the index in the task grid this Cell belongs to
      std::array<int32_t, 3> taskIndex;
      for (int32_t i = 0; i < 3; i++) {
         int32_t n_per_task = globalSize[i] / numTasksPerDim[i];
         int32_t remainder = globalSize[i] % numTasksPerDim[i];

         if (cell[i] < remainder * (n_per_task + 1)) {
            taskIndex[i] = cell[i] / (n_per_task + 1);
         } else {
            taskIndex[i] = remainder + (cell[i] - remainder * (n_per_task + 1)) / n_per_task;
         }
      }

      // Get the task number from the communicator
      std::pair<int32_t, LocalID> retVal;
      FSGRID_MPI_CHECK(MPI_Cart_rank(cartesian3d, taskIndex.data(), &retVal.first),
                       "Unable to find FsGrid rank for global ID ", id, "(coordinates [", cell[0], ", ", cell[1], ", ",
                       cell[2], "])");

      // Determine localID of that cell within the target task
      std::array<FsIndex_t, 3> thatTasksStart;
      std::array<FsIndex_t, 3> thatTaskStorageSize;
      for (int32_t i = 0; i < 3; i++) {
         thatTasksStart[i] = FsGridTools::calcLocalStart(globalSize[i], numTasksPerDim[i], taskIndex[i]);
         thatTaskStorageSize[i] =
             FsGridTools::calcLocalSize(globalSize[i], numTasksPerDim[i], taskIndex[i]) + 2 * stencilSize;
      }

      retVal.second = 0;
      int32_t stride = 1;
      for (int32_t i = 0; i < 3; i++) {
         if (globalSize[i] <= 1) {
            // Collapsed dimension, doesn't contribute.
            retVal.second += 0;
         } else {
            retVal.second += stride * (cell[i] - thatTasksStart[i] + stencilSize);
            stride *= thatTaskStorageSize[i];
         }
      }

      return retVal;
   }

   /*! Transform global cell coordinates into the local domain.
    * If the coordinates are out of bounds, (-1,-1,-1) is returned.
    * \param x The cell's global x coordinate
    * \param y The cell's global y coordinate
    * \param z The cell's global z coordinate
    */
   std::array<FsIndex_t, 3> globalToLocal(FsSize_t x, FsSize_t y, FsSize_t z) {
      std::array<FsIndex_t, 3> retval{(FsIndex_t)x - localStart[0], (FsIndex_t)y - localStart[1],
                                      (FsIndex_t)z - localStart[2]};

      if (retval[0] >= localSize[0] || retval[1] >= localSize[1] || retval[2] >= localSize[2] || retval[0] < 0 ||
          retval[1] < 0 || retval[2] < 0) {
         return {-1, -1, -1};
      }

      return retval;
   }

   /*! Determine the cell's GlobalID from its local x,y,z coordinates
    * \param x The cell's task-local x coordinate
    * \param y The cell's task-local y coordinate
    * \param z The cell's task-local z coordinate
    */
   GlobalID getGlobalIDForCoords(int32_t x, int32_t y, int32_t z) {
      return x + localStart[0] + globalSize[0] * (y + localStart[1]) +
             globalSize[0] * globalSize[1] * (z + localStart[2]);
   }

   /*! Determine the cell's LocalID from its local x,y,z coordinates
    * \param x The cell's task-local x coordinate
    * \param y The cell's task-local y coordinate
    * \param z The cell's task-local z coordinate
    */
   LocalID getLocalIDForCoords(int32_t x, int32_t y, int32_t z) {
      return static_cast<int32_t>(globalSize[2] > 1) * storageSize[0] * storageSize[1] * (stencilSize + z) +
             static_cast<int32_t>(globalSize[1] > 1) * storageSize[0] * (stencilSize + y) +
             static_cast<int32_t>(globalSize[0] > 1) * (stencilSize + x);
   }

   std::array<FsIndex_t, 3>& getTaskPosition() { return taskPosition; }
   const std::array<FsIndex_t, 3>& getTaskPosition() const { return taskPosition; }

   /*! Get the size of the local domain handled by this grid.
    */
   std::array<FsIndex_t, 3>& getLocalSize() { return localSize; }
   const std::array<FsIndex_t, 3>& getLocalSize() const { return localSize; }

   /*! Get the start coordinates of the local domain handled by this grid.
    */
   std::array<FsIndex_t, 3>& getLocalStart() { return localStart; }
   const std::array<FsIndex_t, 3>& getLocalStart() const { return localStart; }

   /*! Get global size of the fsgrid domain
    */
   std::array<FsSize_t, 3>& getGlobalSize() { return globalSize; }
   const std::array<FsSize_t, 3>& getGlobalSize() const { return globalSize; }

   /*! Calculate global cell position (XYZ in global cell space) from local cell coordinates.
    *
    * \param x x-Coordinate, in cells
    * \param y y-Coordinate, in cells
    * \param z z-Coordinate, in cells
    *
    * \return Global cell coordinates
    */
   std::array<FsIndex_t, 3> getGlobalIndices(int64_t x, int64_t y, int64_t z) {
      return {static_cast<int32_t>(localStart[0] + x), static_cast<int32_t>(localStart[1] + y),
              static_cast<int32_t>(localStart[2] + z)};
   }

   /*! Get the rank of this CPU in the FsGrid communicator */
   int32_t getRank() const { return rank; }

   /*! Get the number of ranks in the FsGrid communicator */
   int32_t getSize() const { return numTasksPerDim[0] * numTasksPerDim[1] * numTasksPerDim[2]; }

   int32_t getStencilSize() const { return stencilSize; }

   /*! Get in which directions, if any, this grid is periodic */
   std::array<bool, 3>& getPeriodic() { return periodic; }
   const std::array<bool, 3>& getPeriodic() const { return periodic; }

   /*! Get the decomposition array*/
   std::array<Task_t, 3>& getDecomposition() { return numTasksPerDim; }
   const std::array<Task_t, 3>& getDecomposition() const { return numTasksPerDim; }

   std::array<Task_t, 3>& getStorageSize() { return storageSize; }
   const std::array<Task_t, 3>& getStorageSize() const { return storageSize; }

   size_t getNumTotalStored() const {
      return std::accumulate(storageSize.cbegin(), storageSize.cend(), 1, std::multiplies<>());
   }

   MPI_Comm getComm() const { return cartesian3d; }

   constexpr int32_t computeNeighbourIndex(int32_t x, int32_t y, int32_t z) const {
      return 13 - (x < 0) * 9 + (x >= localSize[0]) * 9 - (y < 0) * 3 + (y >= localSize[1]) * 3 - (z < 0) +
             (z >= localSize[2]);
   }
};
}
