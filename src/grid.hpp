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
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <ios>
#include <iostream>
#include <mpi.h>
#include <sstream>
#include <type_traits>
#include <vector>

/*! Simple cartesian, non-loadbalancing MPI Grid for use with the fieldsolver
 *
 * \param T datastructure containing the field in each cell which this grid manages
 * \param stencil ghost cell width of this grid
 */
template <typename T, int32_t stencil> class FsGrid {
   using FsSize_t = FsGridTools::FsSize_t;
   using FsIndex_t = FsGridTools::FsIndex_t;
   using LocalID = FsGridTools::LocalID;
   using GlobalID = FsGridTools::GlobalID;
   using Task_t = FsGridTools::Task_t;
   template <typename ArrayT> void swapArray(std::array<ArrayT, 3>& array) {
      ArrayT a = array[0];
      array[0] = array[2];
      array[2] = a;
   }

public:
   /*! Constructor for this grid.
    * \param globalSize Cell size of the global simulation domain.
    * \param MPI_Comm The MPI communicator this grid should use.
    * \param isPeriodic An array specifying, for each dimension, whether it is to be treated as periodic.
    */
   FsGrid(std::array<FsSize_t, 3> globalSize, MPI_Comm parentComm, std::array<bool, 3> isPeriodic,
          const std::array<Task_t, 3>& decomposition = {0, 0, 0})
       : globalSize(globalSize) {
      int32_t size;

      // Get parentComm info
      int32_t parentRank;
      FSGRID_MPI_CHECK(MPI_Comm_rank(parentComm, &parentRank), "Couldn't get rank from parent communicator");
      int32_t parentSize;
      FSGRID_MPI_CHECK(MPI_Comm_size(parentComm, &parentSize), "Couldn't get size from parent communicator");

      // If environment variable FSGRID_PROCS is set,
      // use that for determining the number of FS-processes
      size = parentSize;
      if (getenv("FSGRID_PROCS") != NULL) {
         const int32_t fsgridProcs = atoi(getenv("FSGRID_PROCS"));
         if (fsgridProcs > 0 && fsgridProcs < size)
            size = fsgridProcs;
      }

      std::array<Task_t, 3> emptyarr = {0, 0, 0};
      if (decomposition == emptyarr) {
         // If decomposition isn't pre-defined, heuristically choose a good domain decomposition for our field size
         ntasksPerDim = FsGridTools::computeDomainDecomposition(globalSize, size, stencil);
      } else {
         ntasksPerDim = decomposition;
         if (ntasksPerDim[0] * ntasksPerDim[1] * ntasksPerDim[2] != size) {
            std::cerr << "Given decomposition (" << ntasksPerDim[0] << " " << ntasksPerDim[1] << " " << ntasksPerDim[2]
                      << ") does not distribute to the number of tasks given" << std::endl;
            throw std::runtime_error("Given decomposition does not distribute to the number of tasks given");
         }
         ntasksPerDim[0] = decomposition[0];
         ntasksPerDim[1] = decomposition[1];
         ntasksPerDim[2] = decomposition[2];
      }

      // set private array
      periodic = isPeriodic;
      // set temporary int32_t arrays for MPI_Cart_create
      std::array<int32_t, 3> isPeriodicInt, ntasksInt;
      for (size_t i = 0; i < isPeriodic.size(); i++) {
         isPeriodicInt[i] = (int32_t)isPeriodic[i];
         ntasksInt[i] = (int32_t)ntasksPerDim[i];
      }

      // Create a temporary FS subcommunicator for the MPI_Cart_create
      int32_t colorFs = (parentRank < size) ? 1 : MPI_UNDEFINED;
      FSGRID_MPI_CHECK(MPI_Comm_split(parentComm, colorFs, parentRank, &comm1d),
                       "Couldn's split parent communicator to subcommunicators");

      if (colorFs != MPI_UNDEFINED) {
         // Create cartesian communicator. Note, that reorder is false so
         // ranks match the ones in parentComm
         FSGRID_MPI_CHECK(MPI_Cart_create(comm1d, 3, ntasksPerDim.data(), isPeriodicInt.data(), 0, &comm3d),
                          "Creating cartesian communicatior failed when attempting to create FsGrid!");
         FSGRID_MPI_CHECK(MPI_Comm_rank(comm3d, &rank), "Getting rank failed when attempting to create FsGrid!");
         FSGRID_MPI_CHECK(
             MPI_Cart_coords(comm3d, rank, 3, taskPosition.data()), "Rank ", rank,
             " unable to determine own position in cartesian communicator when attempting to create FsGrid!");
      }

      // Create a temporary aux subcommunicator for the (Aux) MPI_Cart_create
      int32_t colorAux =
          (parentRank > (parentSize - 1) % size) ? (parentRank - (parentSize % size)) / size : MPI_UNDEFINED;
      FSGRID_MPI_CHECK(MPI_Comm_split(parentComm, colorAux, parentRank, &comm1d_aux),
                       "Couldn's split parent communicator to auxiliary subcommunicators");

      int32_t rankAux;
      std::array<int32_t, 3> taskPositionAux;

      if (colorAux != MPI_UNDEFINED) {
         // Create an aux cartesian communicator corresponding to the comm3d (but shidted).
         FSGRID_MPI_CHECK(MPI_Cart_create(comm1d_aux, 3, ntasksPerDim.data(), isPeriodicInt.data(), 0, &comm3d_aux),
                          "Creating cartesian communicatior failed when attempting to create FsGrid!");
         FSGRID_MPI_CHECK(MPI_Comm_rank(comm3d_aux, &rankAux), "Getting rank failed when attempting to create FsGrid!");
         FSGRID_MPI_CHECK(
             MPI_Cart_coords(comm3d_aux, rankAux, 3, taskPositionAux.data()), "Rank ", rankAux,
             " unable to determine own position in cartesian communicator when attempting to create FsGrid!");
      }

#ifdef FSGRID_DEBUG
      {
         // All FS ranks send their true comm3d rank and taskPosition data to dest
         MPI_Request* request = new MPI_Request[(parentSize - 1) / size * 2 + 2];
         for (int32_t i = 0; i < (parentSize - 1) / size; i++) {
            int32_t dest =
                (colorFs != MPI_UNDEFINED) ? parentRank + i * size + (parentSize - 1) % size + 1 : MPI_PROC_NULL;
            if (dest >= parentSize)
               dest = MPI_PROC_NULL;
            FSGRID_MPI_CHECK(MPI_Isend(&rank, 1, MPI_INT, dest, 9274, parentComm, &request[2 * i]),
                             "Failed to send comm3d rank");
            FSGRID_MPI_CHECK(MPI_Isend(taskPosition.data(), 3, MPI_INT, dest, 9275, parentComm, &request[2 * i + 1]),
                             "Failed to send comm3d taskPosition");
         }

         // All Aux ranks receive the true comm3d rank and taskPosition data from
         // source and then compare that it matches their aux data
         std::array<int32_t, 3> taskPositionRecv;
         int32_t rankRecv;
         int32_t source = (colorAux != MPI_UNDEFINED)
                              ? parentRank - (parentRank - (parentSize % size)) / size * size - parentSize % size
                              : MPI_PROC_NULL;

         FSGRID_MPI_CHECK(
             MPI_Irecv(&rankRecv, 1, MPI_INT, source, 9274, parentComm, &request[(parentSize - 1) / size * 2]),
             "Failed to receive comm3d rank");
         FSGRID_MPI_CHECK(MPI_Irecv(taskPositionRecv.data(), 3, MPI_INT, source, 9275, parentComm,
                                    &request[(parentSize - 1) / size * 2 + 1]),
                          "Couldn't receive taskPosition");
         FSGRID_MPI_CHECK(MPI_Waitall((parentSize - 1) / size * 2 + 2, request, MPI_STATUS_IGNORE), "Waitall failed");

         if (colorAux != MPI_UNDEFINED) {
            if (rankRecv != rankAux || taskPositionRecv[0] != taskPositionAux[0] ||
                taskPositionRecv[1] != taskPositionAux[1] || taskPositionRecv[2] != taskPositionAux[2]) {
               std::cerr << "Rank: " << parentRank
                         << ". Aux cartesian communicator 'comm3d_aux' does not match with 'comm3d' !" << std::endl;
               throw std::runtime_error("FSGrid aux communicator setup failed.");
            }
         }
         delete[] request;
      }
#endif // FSGRID_DEBUG

      // Set correct task position for non-FS ranks
      if (colorFs == MPI_UNDEFINED) {
         for (int32_t i = 0; i < 3; i++) {
            taskPosition[i] = taskPositionAux[i];
         }
      }

      // Determine size of our local grid
      for (int32_t i = 0; i < 3; i++) {
         localSize[i] = FsGridTools::calcLocalSize(globalSize[i], ntasksPerDim[i], taskPosition[i]);
         localStart[i] = FsGridTools::calcLocalStart(globalSize[i], ntasksPerDim[i], taskPosition[i]);
      }

      if (localSize[0] == 0 || (globalSize[0] > stencil && localSize[0] < stencil) || localSize[1] == 0 ||
          (globalSize[1] > stencil && localSize[1] < stencil) || localSize[2] == 0 ||
          (globalSize[2] > stencil && localSize[2] < stencil)) {
         std::cerr << "FSGrid space partitioning leads to a space that is too small on Rank " << rank << "."
                   << std::endl;
         std::cerr << "Please run with a different number of Tasks, so that space is better divisible." << std::endl;
         throw std::runtime_error("FSGrid too small domains");
      }

      // If non-FS process, set rank to -1 and localSize to zero and return
      if (colorFs == MPI_UNDEFINED) {
         rank = -1;
         localSize[0] = 0;
         localSize[1] = 0;
         localSize[2] = 0;
         comm3d = comm3d_aux;
         comm3d_aux = MPI_COMM_NULL;
         return;
      }

      // Allocate the array of neighbours
      for (int32_t i = 0; i < size; i++) {
         neighbour_index.push_back(MPI_PROC_NULL);
      }
      for (int32_t i = 0; i < 27; i++) {
         neighbour[i] = MPI_PROC_NULL;
      }

      // Get the IDs of the 26 direct neighbours
      for (int32_t x = -1; x <= 1; x++) {
         for (int32_t y = -1; y <= 1; y++) {
            for (int32_t z = -1; z <= 1; z++) {
               std::array<Task_t, 3> neighPosition;

               /*
                * Figure out the coordinates of the neighbours in all three
                * directions
                */
               neighPosition[0] = taskPosition[0] + x;
               if (isPeriodic[0]) {
                  neighPosition[0] += ntasksPerDim[0];
                  neighPosition[0] %= ntasksPerDim[0];
               }

               neighPosition[1] = taskPosition[1] + y;
               if (isPeriodic[1]) {
                  neighPosition[1] += ntasksPerDim[1];
                  neighPosition[1] %= ntasksPerDim[1];
               }

               neighPosition[2] = taskPosition[2] + z;
               if (isPeriodic[2]) {
                  neighPosition[2] += ntasksPerDim[2];
                  neighPosition[2] %= ntasksPerDim[2];
               }

               /*
                * If those coordinates exist, figure out the responsible CPU
                * and store its rank
                */
               if (neighPosition[0] >= 0 && neighPosition[0] < ntasksPerDim[0] && neighPosition[1] >= 0 &&
                   neighPosition[1] < ntasksPerDim[1] && neighPosition[2] >= 0 && neighPosition[2] < ntasksPerDim[2]) {

                  // Calculate the rank
                  int32_t neighRank;
                  FSGRID_MPI_CHECK(MPI_Cart_rank(comm3d, neighPosition.data(), &neighRank), "Rank ", rank,
                                   " can't determine neighbour rank at position [", neighPosition[0], ", ",
                                   neighPosition[1], ", ", neighPosition[2], "]");

                  // Forward lookup table
                  neighbour[(x + 1) * 9 + (y + 1) * 3 + (z + 1)] = neighRank;

                  // Reverse lookup table
                  if (neighRank >= 0 && neighRank < size) {
                     neighbour_index[neighRank] = (char)((x + 1) * 9 + (y + 1) * 3 + (z + 1));
                  }
               } else {
                  neighbour[(x + 1) * 9 + (y + 1) * 3 + (z + 1)] = MPI_PROC_NULL;
               }
            }
         }
      }

      // Allocate local storage array
      size_t totalStorageSize = 1;
      for (int32_t i = 0; i < 3; i++) {
         if (globalSize[i] <= 1) {
            // Collapsed dimension => only one cell thick
            storageSize[i] = 1;
         } else {
            // Size of the local domain + 2* size for the ghost cell stencil
            storageSize[i] = (localSize[i] + stencil * 2);
         }
         totalStorageSize *= storageSize[i];
      }
      data.resize(totalStorageSize);

      MPI_Datatype mpiTypeT;
      FSGRID_MPI_CHECK(MPI_Type_contiguous(sizeof(T), MPI_BYTE, &mpiTypeT), "Failed to create a contiguous data type");
      for (int32_t x = -1; x <= 1; x++) {
         for (int32_t y = -1; y <= 1; y++) {
            for (int32_t z = -1; z <= 1; z++) {
               neighbourSendType[(x + 1) * 9 + (y + 1) * 3 + (z + 1)] = MPI_DATATYPE_NULL;
               neighbourReceiveType[(x + 1) * 9 + (y + 1) * 3 + (z + 1)] = MPI_DATATYPE_NULL;
            }
         }
      }

      // Compute send and receive datatypes
      // loop through the shifts in the different directions
      for (int32_t x = -1; x <= 1; x++) {
         for (int32_t y = -1; y <= 1; y++) {
            for (int32_t z = -1; z <= 1; z++) {
               std::array<int32_t, 3> subarraySize;
               std::array<int32_t, 3> subarrayStart;
               const int32_t shiftId = (x + 1) * 9 + (y + 1) * 3 + (z + 1);

               if ((storageSize[0] == 1 && x != 0) || (storageSize[1] == 1 && y != 0) ||
                   (storageSize[2] == 1 && z != 0) || (x == 0 && y == 0 && z == 0)) {
                  // skip flat dimension for 2 or 1D simulations, and self
                  neighbourSendType[shiftId] = MPI_DATATYPE_NULL;
                  neighbourReceiveType[shiftId] = MPI_DATATYPE_NULL;
                  continue;
               }

               subarraySize[0] = (x == 0) ? localSize[0] : stencil;
               subarraySize[1] = (y == 0) ? localSize[1] : stencil;
               subarraySize[2] = (z == 0) ? localSize[2] : stencil;

               if (x == 0 || x == -1)
                  subarrayStart[0] = stencil;
               else if (x == 1)
                  subarrayStart[0] = storageSize[0] - 2 * stencil;
               if (y == 0 || y == -1)
                  subarrayStart[1] = stencil;
               else if (y == 1)
                  subarrayStart[1] = storageSize[1] - 2 * stencil;
               if (z == 0 || z == -1)
                  subarrayStart[2] = stencil;
               else if (z == 1)
                  subarrayStart[2] = storageSize[2] - 2 * stencil;

               for (int32_t i = 0; i < 3; i++)
                  if (storageSize[i] == 1)
                     subarrayStart[i] = 0;

               std::array<int32_t, 3> swappedStorageSize = {(int32_t)storageSize[0], (int32_t)storageSize[1],
                                                            (int32_t)storageSize[2]};
               swapArray(swappedStorageSize);
               swapArray(subarraySize);
               swapArray(subarrayStart);
               FSGRID_MPI_CHECK(MPI_Type_create_subarray(3, swappedStorageSize.data(), subarraySize.data(),
                                                         subarrayStart.data(), MPI_ORDER_C, mpiTypeT,
                                                         &(neighbourSendType[shiftId])),
                                "Failed to create a subarray type");

               if (x == 1)
                  subarrayStart[0] = 0;
               else if (x == 0)
                  subarrayStart[0] = stencil;
               else if (x == -1)
                  subarrayStart[0] = storageSize[0] - stencil;
               if (y == 1)
                  subarrayStart[1] = 0;
               else if (y == 0)
                  subarrayStart[1] = stencil;
               else if (y == -1)
                  subarrayStart[1] = storageSize[1] - stencil;
               if (z == 1)
                  subarrayStart[2] = 0;
               else if (z == 0)
                  subarrayStart[2] = stencil;
               else if (z == -1)
                  subarrayStart[2] = storageSize[2] - stencil;
               for (int32_t i = 0; i < 3; i++)
                  if (storageSize[i] == 1)
                     subarrayStart[i] = 0;

               swapArray(subarrayStart);
               FSGRID_MPI_CHECK(MPI_Type_create_subarray(3, swappedStorageSize.data(), subarraySize.data(),
                                                         subarrayStart.data(), MPI_ORDER_C, mpiTypeT,
                                                         &(neighbourReceiveType[shiftId])),
                                "Failed to create a subarray type");
            }
         }
      }

      for (int32_t i = 0; i < 27; i++) {
         if (neighbourReceiveType[i] != MPI_DATATYPE_NULL)
            FSGRID_MPI_CHECK(MPI_Type_commit(&(neighbourReceiveType[i])), "Failed to commit MPI type");
         if (neighbourSendType[i] != MPI_DATATYPE_NULL)
            FSGRID_MPI_CHECK(MPI_Type_commit(&(neighbourSendType[i])), "Failed to commit MPI type");
      }
   }

   std::vector<T>& getData() { return this->data; }

   void copyData(FsGrid& other) {
      this->data = other.getData(); // Copy assignment
   }

   /*! Finalize instead of destructor, as the MPI calls fail after the main program called MPI_Finalize().
    *  Cleans up the cartesian communicator
    */
   void finalize() {
      // If not a non-FS process
      if (rank != -1) {
         for (int32_t i = 0; i < 27; i++) {
            if (neighbourReceiveType[i] != MPI_DATATYPE_NULL)
               FSGRID_MPI_CHECK(MPI_Type_free(&(neighbourReceiveType[i])), "Failed to free MPI type");
            if (neighbourSendType[i] != MPI_DATATYPE_NULL)
               FSGRID_MPI_CHECK(MPI_Type_free(&(neighbourSendType[i])), "Failed to free MPI type");
         }
      }

      if (comm3d != MPI_COMM_NULL)
         FSGRID_MPI_CHECK(MPI_Comm_free(&comm3d), "Failed to free MPI comm3d");
      if (comm3d_aux != MPI_COMM_NULL)
         FSGRID_MPI_CHECK(MPI_Comm_free(&comm3d_aux), "Failed to free MPI comm3d aux");
      if (comm1d != MPI_COMM_NULL)
         FSGRID_MPI_CHECK(MPI_Comm_free(&comm1d), "Failed to free MPI comm1d");
      if (comm1d_aux != MPI_COMM_NULL)
         FSGRID_MPI_CHECK(MPI_Comm_free(&comm1d_aux), "Failed to free MPI comm1d aux");
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
         int32_t n_per_task = globalSize[i] / ntasksPerDim[i];
         int32_t remainder = globalSize[i] % ntasksPerDim[i];

         if (cell[i] < remainder * (n_per_task + 1)) {
            taskIndex[i] = cell[i] / (n_per_task + 1);
         } else {
            taskIndex[i] = remainder + (cell[i] - remainder * (n_per_task + 1)) / n_per_task;
         }
      }

      // Get the task number from the communicator
      std::pair<int32_t, LocalID> retVal;
      FSGRID_MPI_CHECK(MPI_Cart_rank(comm3d, taskIndex.data(), &retVal.first),
                       "Unable to find FsGrid rank for global ID ", id, "(coordinates [", cell[0], ", ", cell[1], ", ",
                       cell[2], "])");

      // Determine localID of that cell within the target task
      std::array<FsIndex_t, 3> thatTasksStart;
      std::array<FsIndex_t, 3> thatTaskStorageSize;
      for (int32_t i = 0; i < 3; i++) {
         thatTasksStart[i] = FsGridTools::calcLocalStart(globalSize[i], ntasksPerDim[i], taskIndex[i]);
         thatTaskStorageSize[i] = FsGridTools::calcLocalSize(globalSize[i], ntasksPerDim[i], taskIndex[i]) + 2 * stencil;
      }

      retVal.second = 0;
      int32_t stride = 1;
      for (int32_t i = 0; i < 3; i++) {
         if (globalSize[i] <= 1) {
            // Collapsed dimension, doesn't contribute.
            retVal.second += 0;
         } else {
            retVal.second += stride * (cell[i] - thatTasksStart[i] + stencil);
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
      std::array<FsIndex_t, 3> retval;
      retval[0] = (FsIndex_t)x - localStart[0];
      retval[1] = (FsIndex_t)y - localStart[1];
      retval[2] = (FsIndex_t)z - localStart[2];

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
   GlobalID GlobalIDForCoords(int32_t x, int32_t y, int32_t z) {
      return x + localStart[0] + globalSize[0] * (y + localStart[1]) +
             globalSize[0] * globalSize[1] * (z + localStart[2]);
   }
   /*! Determine the cell's LocalID from its local x,y,z coordinates
    * \param x The cell's task-local x coordinate
    * \param y The cell's task-local y coordinate
    * \param z The cell's task-local z coordinate
    */
   LocalID LocalIDForCoords(int32_t x, int32_t y, int32_t z) {
      LocalID index = 0;
      if (globalSize[2] > 1) {
         index += storageSize[0] * storageSize[1] * (stencil + z);
      }
      if (globalSize[1] > 1) {
         index += storageSize[0] * (stencil + y);
      }
      if (globalSize[0] > 1) {
         index += stencil + x;
      }

      return index;
   }

   /*! Perform ghost cell communication.
    */
   void updateGhostCells() {

      if (rank == -1)
         return;

      // TODO, faster with simultaneous isends& ireceives?
      std::array<MPI_Request, 27> receiveRequests;
      std::array<MPI_Request, 27> sendRequests;

      for (int32_t i = 0; i < 27; i++) {
         receiveRequests[i] = MPI_REQUEST_NULL;
         sendRequests[i] = MPI_REQUEST_NULL;
      }

      for (int32_t x = -1; x <= 1; x++) {
         for (int32_t y = -1; y <= 1; y++) {
            for (int32_t z = -1; z <= 1; z++) {
               int32_t shiftId = (x + 1) * 9 + (y + 1) * 3 + (z + 1);
               int32_t receiveId = (1 - x) * 9 + (1 - y) * 3 + (1 - z);
               if (neighbour[receiveId] != MPI_PROC_NULL && neighbourSendType[shiftId] != MPI_DATATYPE_NULL) {
                  FSGRID_MPI_CHECK(MPI_Irecv(data.data(), 1, neighbourReceiveType[shiftId], neighbour[receiveId],
                                             shiftId, comm3d, &(receiveRequests[shiftId])),
                                   "Rank ", rank, " failed to receive data from neighbor ", receiveId, " with rank ",
                                   neighbour[receiveId]);
               }
            }
         }
      }

      for (int32_t x = -1; x <= 1; x++) {
         for (int32_t y = -1; y <= 1; y++) {
            for (int32_t z = -1; z <= 1; z++) {
               int32_t shiftId = (x + 1) * 9 + (y + 1) * 3 + (z + 1);
               int32_t sendId = shiftId;
               if (neighbour[sendId] != MPI_PROC_NULL && neighbourSendType[shiftId] != MPI_DATATYPE_NULL) {
                  FSGRID_MPI_CHECK(MPI_Isend(data.data(), 1, neighbourSendType[shiftId], neighbour[sendId], shiftId,
                                             comm3d, &(sendRequests[shiftId])),
                                   "Rank ", rank, " failed to send data to neighbor ", sendId, " with rank ",
                                   neighbour[sendId]);
               }
            }
         }
      }
      FSGRID_MPI_CHECK(MPI_Waitall(27, receiveRequests.data(), MPI_STATUSES_IGNORE),
                       "Synchronization at ghost cell update failed");
      FSGRID_MPI_CHECK(MPI_Waitall(27, sendRequests.data(), MPI_STATUSES_IGNORE),
                       "Synchronization at ghost cell update failed");
   }

   /*! Get the size of the local domain handled by this grid.
    */
   std::array<FsIndex_t, 3>& getLocalSize() { return localSize; }

   /*! Get the start coordinates of the local domain handled by this grid.
    */
   std::array<FsIndex_t, 3>& getLocalStart() { return localStart; }

   /*! Get global size of the fsgrid domain
    */
   std::array<FsSize_t, 3>& getGlobalSize() { return globalSize; }

   /*! Calculate global cell position (XYZ in global cell space) from local cell coordinates.
    *
    * \param x x-Coordinate, in cells
    * \param y y-Coordinate, in cells
    * \param z z-Coordinate, in cells
    *
    * \return Global cell coordinates
    */
   std::array<FsIndex_t, 3> getGlobalIndices(int64_t x, int64_t y, int64_t z) {
      std::array<FsIndex_t, 3> retval;
      retval[0] = localStart[0] + x;
      retval[1] = localStart[1] + y;
      retval[2] = localStart[2] + z;

      return retval;
   }

   /*! Get a reference to the field data in a cell
    * \param x x-Coordinate, in cells
    * \param y y-Coordinate, in cells
    * \param z z-Coordinate, in cells
    * \return A reference to cell data in the given cell
    */
   T* get(int32_t x, int32_t y, int32_t z) {

      // Keep track which neighbour this cell actually belongs to (13 = ourself)
      int32_t isInNeighbourDomain = 13;
      int32_t coord_shift[3] = {0, 0, 0};
      if (x < 0) {
         isInNeighbourDomain -= 9;
         coord_shift[0] = 1;
      }
      if (x >= localSize[0]) {
         isInNeighbourDomain += 9;
         coord_shift[0] = -1;
      }
      if (y < 0) {
         isInNeighbourDomain -= 3;
         coord_shift[1] = 1;
      }
      if (y >= localSize[1]) {
         isInNeighbourDomain += 3;
         coord_shift[1] = -1;
      }
      if (z < 0) {
         isInNeighbourDomain -= 1;
         coord_shift[2] = 1;
      }
      if (z >= localSize[2]) {
         isInNeighbourDomain += 1;
         coord_shift[2] = -1;
      }

      // Santiy-Check that the requested cell is actually inside our domain
      // TODO: ugh, this is ugly.
#ifdef FSGRID_DEBUG
      bool inside = true;
      if (localSize[0] <= 1 && !periodic[0]) {
         if (x != 0) {
            std::cerr << "x != 0 despite non-periodic x-axis with only one cell." << std::endl;
            inside = false;
         }
      } else {
         if (x < -stencil || x >= localSize[0] + stencil) {
            std::cerr << "x = " << x << " is outside of [ " << -stencil << ", " << localSize[0] + stencil << "[!"
                      << std::endl;
            inside = false;
         }
      }

      if (localSize[1] <= 1 && !periodic[1]) {
         if (y != 0) {
            std::cerr << "y != 0 despite non-periodic y-axis with only one cell." << std::endl;
            inside = false;
         }
      } else {
         if (y < -stencil || y >= localSize[1] + stencil) {
            std::cerr << "y = " << y << " is outside of [ " << -stencil << ", " << localSize[1] + stencil << "[!"
                      << std::endl;
            inside = false;
         }
      }

      if (localSize[2] <= 1 && !periodic[2]) {
         if (z != 0) {
            std::cerr << "z != 0 despite non-periodic z-axis with only one cell." << std::endl;
            inside = false;
         }
      } else {
         if (z < -stencil || z >= localSize[2] + stencil) {
            inside = false;
            std::cerr << "z = " << z << " is outside of [ " << -stencil << ", " << localSize[2] + stencil << "[!"
                      << std::endl;
         }
      }
      if (!inside) {
         std::cerr << "Out-of bounds access in FsGrid::get! Expect weirdness." << std::endl;
         return NULL;
      }
#endif // FSGRID_DEBUG

      if (isInNeighbourDomain != 13) {

         // Check if the corresponding neighbour exists
         if (neighbour[isInNeighbourDomain] == MPI_PROC_NULL) {
            // Neighbour doesn't exist, we must be an outer boundary cell
            // (or something is quite wrong)
            return NULL;

         } else if (neighbour[isInNeighbourDomain] == rank) {
            // For periodic boundaries, where the neighbour is actually ourself,
            // return our own actual cell instead of the ghost
            x += coord_shift[0] * localSize[0];
            y += coord_shift[1] * localSize[1];
            z += coord_shift[2] * localSize[2];
         }
         // Otherwise we return the ghost cell
      }
      LocalID index = LocalIDForCoords(x, y, z);

      return &data[index];
   }

   T* get(LocalID id) {
      if (id < 0 || (size_t)id > data.size()) {
         std::cerr << "Out-of-bounds access in FsGrid::get!" << std::endl
                   << "(LocalID = " << id << ", but storage space is " << data.size() << ". Expect weirdness."
                   << std::endl;
         return NULL;
      }
      return &data[id];
   }

   /*! Get the physical coordinates in the global simulation space for
    * the given cell.
    *
    * \param x local x-Coordinate, in cells
    * \param y local y-Coordinate, in cells
    * \param z local z-Coordinate, in cells
    */
   std::array<double, 3> getPhysicalCoords(int32_t x, int32_t y, int32_t z) {
      std::array<double, 3> coords;
      coords[0] = physicalGlobalStart[0] + (localStart[0] + x) * DX;
      coords[1] = physicalGlobalStart[1] + (localStart[1] + y) * DY;
      coords[2] = physicalGlobalStart[2] + (localStart[2] + z) * DZ;

      return coords;
   }

   /*! Debugging output helper function. Allows for nicely formatted printing
    * of grid contents. Since the grid data format is varying, the actual
    * printing should be done in a lambda passed to this function. Example usage
    * to plot |B|:
    *
    * perBGrid.debugOutput([](const std::array<Real, fsgrids::bfield::N_BFIELD>& a)->void{
    *     cerr << sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) << ", ";
    * });
    *
    * \param func Function pointer (or lambda) which is called with a cell reference,
    * in order. Use std::cerr in it to print desired value.
    */
   void debugOutput(void(func)(const T&)) {
      int32_t xmin = 0, xmax = 1;
      int32_t ymin = 0, ymax = 1;
      int32_t zmin = 0, zmax = 1;
      if (localSize[0] > 1) {
         xmin = -stencil;
         xmax = localSize[0] + stencil;
      }
      if (localSize[1] > 1) {
         ymin = -stencil;
         ymax = localSize[1] + stencil;
      }
      if (localSize[2] > 1) {
         zmin = -stencil;
         zmax = localSize[2] + stencil;
      }
      for (int32_t z = zmin; z < zmax; z++) {
         for (int32_t y = ymin; y < ymax; y++) {
            for (int32_t x = xmin; x < xmax; x++) {
               func(*get(x, y, z));
            }
            std::cerr << std::endl;
         }
         std::cerr << " - - - - - - - - " << std::endl;
      }
   }

   /*! Get the rank of this CPU in the FsGrid communicator */
   int32_t getRank() { return rank; }

   /*! Get the number of ranks in the FsGrid communicator */
   int32_t getSize() { return ntasksPerDim[0] * ntasksPerDim[1] * ntasksPerDim[2]; }

   /*! Get in which directions, if any, this grid is periodic */
   std::array<bool, 3>& getPeriodic() { return periodic; }

   /*! Perform an MPI_Allreduce with this grid's internal communicator
    * Function syntax is identical to MPI_Allreduce, except the final (communicator
    * argument will not be needed) */
   int32_t Allreduce(void* sendbuf, void* recvbuf, int32_t count, MPI_Datatype datatype, MPI_Op op) {

      // If a normal FS-rank
      if (rank != -1) {
         return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm3d);
      }
      // If a non-FS rank, no need to communicate
      else {
         int32_t datatypeSize;
         MPI_Type_size(datatype, &datatypeSize);
         for (int32_t i = 0; i < count * datatypeSize; ++i)
            ((char*)recvbuf)[i] = ((char*)sendbuf)[i];
         return MPI_ERR_RANK; // This is ok for a non-FS rank
      }
   }

   /*! Get the decomposition array*/
   std::array<Task_t, 3>& getDecomposition() { return ntasksPerDim; }

   std::string display() const {
      std::stringstream ss;
      std::ios defaultState(nullptr);
      defaultState.copyfmt(ss);

      auto pushContainerValues = [&ss, &defaultState](auto container, bool asByteStr = false, uint32_t nPerLine = 0,
                                                      uint32_t numTabs = 2) {
         nPerLine = nPerLine == 0 ? container.size() : nPerLine;
         const uint32_t numBytes = sizeof(decltype(container[0]));
         if (asByteStr) {
            ss << std::hex << std::setfill('0');
         }

         uint32_t i = 1;
         for (const auto& v : container) {
            if (asByteStr) {
               ss << "0x" << std::setw(2 * numBytes);
               if constexpr (std::is_integral_v<typename std::remove_reference<decltype(v)>::type>) {
                  ss << static_cast<int64_t>(
                      static_cast<std::make_unsigned_t<typename std::remove_reference<decltype(v)>::type>>(v));
               } else {
                  ss << v;
               }
            } else {
               ss << v;
            }

            if (i < container.size()) {
               ss << ", ";
            }
            if (i % nPerLine == 0 && i < container.size()) {
               ss << "\n";
               for (uint32_t j = 0; j < numTabs; j++) {
                  ss << "\t";
               }
            }

            i++;
         }

         ss.copyfmt(defaultState);
      };

      auto pushMPIComm = [&ss, &defaultState](auto prefix, auto comm, auto newliner) {
         ss << prefix;
         if (comm == MPI_COMM_NULL) {
            ss << "MPI_COMM_NULL";
         } else {
            int rank = 0;
            FSGRID_MPI_CHECK(MPI_Comm_rank(comm, &rank), "Failed to get rank from comm ", comm);
            int size = 0;
            FSGRID_MPI_CHECK(MPI_Comm_size(comm, &size), "Failed to get size from comm ", comm);
            ss << newliner;
            ss << "comm rank: ";
            if (rank != MPI_UNDEFINED) {
               ss << rank;
            } else {
               ss << "MPI_UNDEFINED";
            }
            ss << newliner;
            ss << "comm size: " << rank;

            MPI_Group group = MPI_GROUP_NULL;
            FSGRID_MPI_CHECK(MPI_Comm_group(comm, &group), "Failed to get group from comm ", comm);
            if (group != MPI_GROUP_NULL) {
               int rank = 0;
               FSGRID_MPI_CHECK(MPI_Group_rank(group, &rank), "Failed to get rank from group ", group);
               int size = 0;
               FSGRID_MPI_CHECK(MPI_Group_size(group, &size), "Failed to get size from group ", group);

               ss << newliner;
               ss << "group rank: ";
               if (rank != MPI_UNDEFINED) {
                  ss << rank;
               } else {
                  ss << "MPI_UNDEFINED";
               }
               ss << newliner;
               ss << "group size: " << size;
            }

            MPI_Group remotegroup = MPI_GROUP_NULL;
            FSGRID_MPI_CHECK(MPI_Comm_remote_group(comm, &remotegroup), "Failed to get remotegroup from comm ", comm);
            if (remotegroup != MPI_GROUP_NULL) {
               int rank = 0;
               FSGRID_MPI_CHECK(MPI_Group_rank(remotegroup, &rank), "Failed to get rank from remotegroup ",
                                remotegroup);
               int size = 0;
               FSGRID_MPI_CHECK(MPI_Group_size(remotegroup, &size), "Failed to get size from remotegroup ",
                                remotegroup);

               ss << newliner;
               ss << "remotegroup rank: ";
               if (rank != MPI_UNDEFINED) {
                  ss << rank;
               } else {
                  ss << "MPI_UNDEFINED";
               }
               ss << newliner;
               ss << "remotegroup size: " << size;
            }

            int remotesize = 0;
            FSGRID_MPI_CHECK(MPI_Comm_remote_size(comm, &remotesize), "Failed to get remotesize from comm ", comm);
            ss << newliner;
            ss << "remotesize: " << remotesize;

            int isInterComm = 0;
            FSGRID_MPI_CHECK(MPI_Comm_test_inter(comm, &isInterComm), "Failed to get intecomm flag from comm ", comm);
            ss << newliner;
            ss << "is intercomm: " << isInterComm;
         }
      };

      ss << "{";

      pushMPIComm("\n\tcomm1d: ", comm1d, "\n\t\t");
      pushMPIComm("\n\tcomm1d_aux: ", comm1d_aux, "\n\t\t");
      pushMPIComm("\n\tcomm3d: ", comm3d, "\n\t\t");
      pushMPIComm("\n\tcomm3d_aux: ", comm3d_aux, "\n\t\t");
      ss << "\n\trank: " << rank;
      ss << "\n\tneigbour: [\n\t\t";
      pushContainerValues(neighbour, true, 9);
      ss << "\n\t]";
      ss << "\n\tneigbour_index: [\n\t\t";
      pushContainerValues(neighbour_index, true, 9);
      ss << "\n\t]";
      ss << "\n\tntasksPerDim: [\n\t\t";
      pushContainerValues(ntasksPerDim);
      ss << "\n\t]";
      ss << "\n\ttaskPosition: [\n\t\t";
      pushContainerValues(taskPosition);
      ss << "\n\t]";
      ss << "\n\tperiodic: [\n\t\t";
      pushContainerValues(periodic);
      ss << "\n\t]";
      ss << "\n\tglobalSize: [\n\t\t";
      pushContainerValues(globalSize);
      ss << "\n\t]";
      ss << "\n\tlocalSize: [\n\t\t";
      pushContainerValues(localSize);
      ss << "\n\t]";
      ss << "\n\tlocalStart: [\n\t\t";
      pushContainerValues(localStart);
      ss << "\n\t]";
      ss << "\n\tneigbourSendType: [";
      for (const auto& v : getMPITypes(true)) {
         ss << "\n\t\t" << v.display("\n\t\t");
      }
      ss << "\n\t]";
      ss << "\n\tneighbourReceiveType: [";
      for (const auto& v : getMPITypes(false)) {
         ss << "\n\t\t" << v.display("\n\t\t");
      }
      ss << "\n\t]";
      ss << "\n\tinfo on data:";
      ss << "\n\t\tcapacity: " << data.capacity();
      ss << "\n\t\tsize: " << data.size();
      ss << "\n\t\tdata.front: [\n\t\t\t";
      pushContainerValues(data.front());
      ss << "\n\t\t]";
      ss << "\n\t\tdata.back: [\n\t\t\t";
      pushContainerValues(data.back());
      ss << "\n\t\t]";

      ss << "\n}";

      return ss.str();
   }

   struct MPITypeMetaData {
      int combiner = -1;
      std::vector<int> integers;
      std::vector<MPI_Aint> addresses;
      std::vector<MPI_Datatype> dataTypes;

      std::string display(std::string newliner) const {
         std::stringstream ss;
         std::ios defaultState(nullptr);
         defaultState.copyfmt(ss);

         auto pushContainerValues = [&ss, &defaultState, &newliner](auto container, bool asByteStr = false,
                                                                    uint32_t nPerLine = 0, uint32_t numTabs = 2) {
            nPerLine = nPerLine == 0 ? container.size() : nPerLine;
            const uint32_t numBytes = sizeof(decltype(container[0]));
            if (asByteStr) {
               ss << std::hex << std::setfill('0') << std::uppercase;
            }

            uint32_t i = 1;
            for (const auto& v : container) {
               if (asByteStr) {
                  ss << "0x" << std::setw(2 * numBytes);
                  if constexpr (std::is_integral_v<typename std::remove_reference<decltype(v)>::type>) {
                     ss << static_cast<int64_t>(
                         static_cast<std::make_unsigned_t<typename std::remove_reference<decltype(v)>::type>>(v));
                  } else {
                     ss << v;
                  }
               } else {
                  ss << v;
               }

               if (i < container.size()) {
                  ss << ", ";
               }
               if (i % nPerLine == 0 && i < container.size()) {
                  ss << newliner;
                  for (uint32_t j = 0; j < numTabs; j++) {
                     ss << "\t";
                  }
               }

               i++;
            }

            ss.copyfmt(defaultState);
         };

         ss << "{";
         ss << newliner << "\tcombiner :" << combiner;
         ss << newliner << "\tintegers: [" << newliner << "\t\t";
         pushContainerValues(integers, false, 9);
         ss << newliner << "\t]";
         ss << newliner << "\taddresses: [" << newliner << "\t\t";
         pushContainerValues(addresses, true, 9);
         ss << newliner << "\t]";
         ss << newliner << "\tdata types: [" << newliner << "\t\t";
         pushContainerValues(dataTypes, true, 9);
         ss << newliner << "\t]";
         ss << newliner << "}";

         return ss.str();
      }
   };

   std::array<MPITypeMetaData, 27> getMPITypes(bool send) const {
      const auto typeVec = send ? neighbourSendType : neighbourReceiveType;
      std::array<MPITypeMetaData, 27> metadatas;
      for (size_t i = 0; i < typeVec.size(); i++) {
         const auto mpiType = typeVec[i];

         if (mpiType == MPI_DATATYPE_NULL) {
            continue;
         }

         int numIntegers = 0;
         int numAddresses = 0;
         int numDataTypes = 0;
         FSGRID_MPI_CHECK(
             MPI_Type_get_envelope(mpiType, &numIntegers, &numAddresses, &numDataTypes, &metadatas[i].combiner),
             "Failed to get envelope for type ", mpiType);

         metadatas[i].integers.resize(numIntegers);
         metadatas[i].addresses.resize(numAddresses);
         metadatas[i].dataTypes.resize(numDataTypes);
         FSGRID_MPI_CHECK(MPI_Type_get_contents(mpiType, numIntegers, numAddresses, numDataTypes,
                                                metadatas[i].integers.data(), metadatas[i].addresses.data(),
                                                metadatas[i].dataTypes.data()),
                          "Failed to get type contents for type ", mpiType);
      }

      return metadatas;
   }

   /*! Physical grid spacing and physical coordinate space start.
    */
   double DX = 0.0;
   double DY = 0.0;
   double DZ = 0.0;
   std::array<double, 3> physicalGlobalStart = {};

private:
   //! MPI Cartesian communicator used in this grid
   MPI_Comm comm1d = MPI_COMM_NULL;
   MPI_Comm comm1d_aux = MPI_COMM_NULL;
   MPI_Comm comm3d = MPI_COMM_NULL;
   MPI_Comm comm3d_aux = MPI_COMM_NULL;
   //!< This task's rank in the communicator
   int32_t rank;
   //!< Tasks of the 26 neighbours (plus ourselves)
   std::array<int32_t, 27> neighbour = {};
   //!< Lookup table from rank to index in the neighbour array
   std::vector<char> neighbour_index = {};

   // We have, fundamentally, two different coordinate systems we're dealing with:
   // 1) Task grid in the MPI_Cartcomm

   //!< Number of tasks in each direction
   std::array<Task_t, 3> ntasksPerDim = {};
   //!< This task's position in the 3d task grid
   std::array<Task_t, 3> taskPosition = {};

   // 2) Cell numbers in global and local view

   //!< Information about whether a given direction is periodic
   std::array<bool, 3> periodic = {};
   //!< Global size of the simulation space, in cells
   std::array<FsSize_t, 3> globalSize = {};
   //!< Local size of simulation space handled by this task (without ghost cells)
   std::array<FsIndex_t, 3> localSize = {};
   //!< Local size of simulation space handled by this task (including ghost cells)
   std::array<FsIndex_t, 3> storageSize = {};
   //!< Offset of the local coordinate system against the global one
   std::array<FsIndex_t, 3> localStart = {};
   //!< Datatype for sending data
   std::array<MPI_Datatype, 27> neighbourSendType = {};
   //!< Datatype for receiving data
   std::array<MPI_Datatype, 27> neighbourReceiveType = {};

   //! Actual storage of field data
   std::vector<T> data = {};
};
