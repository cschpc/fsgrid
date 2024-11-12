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
#include "coordinates.hpp"
#include "tools.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <sstream>
#include <type_traits>
#include <vector>

namespace fsgrid_detail {
using FsSize_t = FsGridTools::FsSize_t;
using FsIndex_t = FsGridTools::FsIndex_t;
using LocalID = FsGridTools::LocalID;
using GlobalID = FsGridTools::GlobalID;
using Task_t = FsGridTools::Task_t;

// Assumes x, y and z to belong to set [-1, 0, 1]
// returns a value in (inclusive) range [0, 26]
constexpr static uint32_t xyzToLinear(int32_t x, int32_t y, int32_t z) {
   return static_cast<uint32_t>((x + 1) * 9 + (y + 1) * 3 + (z + 1));
}

// These assume i to be in (inclusive) range [0, 26]
// returns a value from the set [-1, 0, 1]
constexpr static int32_t linearToX(uint32_t i) { return static_cast<int32_t>(i) / 9 - 1; }
constexpr static int32_t linearToY(uint32_t i) { return (static_cast<int32_t>(i) % 9) / 3 - 1; }
constexpr static int32_t linearToZ(uint32_t i) { return static_cast<int32_t>(i) % 3 - 1; }

static std::array<int32_t, 27> mapNeigbourIndexToRank(const std::array<Task_t, 3>& taskPosition,
                                                      const std::array<Task_t, 3>& numTasksPerDim,
                                                      const std::array<bool, 3>& periodic, MPI_Comm comm,
                                                      int32_t rank) {
   auto calculateNeighbourRank = [&](uint32_t neighbourIndex) {
      auto calculateNeighbourPosition = [&](uint32_t neighbourIndex, uint32_t i) {
         const auto pos3D =
             i == 0 ? linearToX(neighbourIndex) : (i == 1 ? linearToY(neighbourIndex) : linearToZ(neighbourIndex));
         const auto nonPeriodicPos = taskPosition[i] + pos3D;
         return periodic[i] ? (nonPeriodicPos + numTasksPerDim[i]) % numTasksPerDim[i] : nonPeriodicPos;
      };

      const std::array<Task_t, 3> neighbourPosition = {
          calculateNeighbourPosition(neighbourIndex, 0),
          calculateNeighbourPosition(neighbourIndex, 1),
          calculateNeighbourPosition(neighbourIndex, 2),
      };

      const bool taskPositionWithinLimits = numTasksPerDim[0] > neighbourPosition[0] && neighbourPosition[0] >= 0 &&
                                            numTasksPerDim[1] > neighbourPosition[1] && neighbourPosition[1] >= 0 &&
                                            numTasksPerDim[2] > neighbourPosition[2] && neighbourPosition[2] >= 0;

      if (taskPositionWithinLimits) {
         int32_t neighbourRank;
         FSGRID_MPI_CHECK(MPI_Cart_rank(comm, neighbourPosition.data(), &neighbourRank), "Rank ", rank,
                          " can't determine neighbour rank at position [", neighbourPosition[0], ", ",
                          neighbourPosition[1], ", ", neighbourPosition[2], "]");
         return neighbourRank;
      } else {
         return MPI_PROC_NULL;
      }
   };

   std::array<int32_t, 27> ranks;
   if (rank == -1) {
      ranks.fill(MPI_PROC_NULL);
   } else {
      std::generate(ranks.begin(), ranks.end(),
                    [&calculateNeighbourRank, n = 0u]() mutable { return calculateNeighbourRank(n++); });
   }
   return ranks;
}

static std::vector<char> mapNeighbourRankToIndex(const std::array<int32_t, 27>& indexToRankMap, int32_t numRanks) {
   std::vector<char> indices(static_cast<size_t>(numRanks), MPI_PROC_NULL);
   std::for_each(indexToRankMap.cbegin(), indexToRankMap.cend(), [&indices, &numRanks, n = 0u](int32_t rank) mutable {
      if (rank >= 0 && numRanks > rank) {
         indices[static_cast<FsSize_t>(rank)] = static_cast<int8_t>(n);
      }
      n++;
   });
   return indices;
}

static int32_t getFSCommSize(int32_t parentCommSize) {
   const auto envVar = getenv("FSGRID_PROCS");
   const auto fsgridProcs = envVar != NULL ? atoi(envVar) : 0;
   return parentCommSize > fsgridProcs && fsgridProcs > 0 ? fsgridProcs : parentCommSize;
}

static int32_t getCommRank(MPI_Comm parentComm) {
   int32_t parentRank;
   FSGRID_MPI_CHECK(MPI_Comm_rank(parentComm, &parentRank), "Couldn't get rank from parent communicator");
   return parentRank;
}

static int32_t getCommSize(MPI_Comm parentComm) {
   int32_t parentCommSize;
   FSGRID_MPI_CHECK(MPI_Comm_size(parentComm, &parentCommSize), "Couldn't get size from parent communicator");
   return parentCommSize;
}

constexpr static int32_t computeColourFs(int32_t parentRank, int32_t numRanks) {
   return (parentRank < numRanks) ? 1 : MPI_UNDEFINED;
}

constexpr static int32_t computeColourAux(int32_t parentRank, int32_t parentCommSize, int32_t numRanks) {
   return (parentRank > (parentCommSize - 1) % numRanks) ? (parentRank - (parentCommSize % numRanks)) / numRanks
                                                         : MPI_UNDEFINED;
}

static MPI_Comm createCartesianCommunicator(MPI_Comm parentComm, const std::array<Task_t, 3>& numTasksPerDim,
                                            const std::array<bool, 3>& isPeriodic) {
   const auto parentRank = getCommRank(parentComm);
   const auto parentSize = getCommSize(parentComm);
   const auto numFsRanks = getFSCommSize(parentSize);
   const auto colourFs = computeColourFs(parentRank, numFsRanks);
   const auto colour = colourFs != MPI_UNDEFINED ? colourFs : computeColourAux(parentRank, parentSize, numFsRanks);

   MPI_Comm comm;
   FSGRID_MPI_CHECK(MPI_Comm_split(parentComm, colour, parentRank, &comm),
                    "Couldn's split parent communicator to subcommunicators");

   const std::array<int32_t, 3> pi = {
       isPeriodic[0],
       isPeriodic[1],
       isPeriodic[2],
   };
   MPI_Comm comm3d;
   FSGRID_MPI_CHECK(MPI_Cart_create(comm, 3, numTasksPerDim.data(), pi.data(), 0, &comm3d),
                    "Creating cartesian communicatior failed when attempting to create FsGrid!");
   FSGRID_MPI_CHECK(MPI_Comm_free(&comm), "Failed to free MPI comm");

   return comm3d;
}

static int32_t getCartesianRank(MPI_Comm parentComm, MPI_Comm cartesianComm) {
   const auto parentRank = getCommRank(parentComm);
   const auto parentSize = getCommSize(parentComm);
   const auto numFsRanks = getFSCommSize(parentSize);
   const auto colourFs = computeColourFs(parentRank, numFsRanks);
   return colourFs != MPI_UNDEFINED ? getCommRank(cartesianComm) : -1;
}

static std::array<int32_t, 3> getTaskPosition(MPI_Comm comm) {
   std::array<int32_t, 3> taskPos;
   const int rank = getCommRank(comm);
   FSGRID_MPI_CHECK(MPI_Cart_coords(comm, rank, taskPos.size(), taskPos.data()), "Rank ", rank,
                    " unable to determine own position in cartesian communicator when attempting to create FsGrid!");
   return taskPos;
}

template <typename T>
static std::array<MPI_Datatype, 27> generateMPITypes(const std::array<FsIndex_t, 3>& storageSize,
                                                     const std::array<FsIndex_t, 3>& localSize, int32_t stencilSize,
                                                     bool generateForSend) {
   MPI_Datatype baseType;
   FSGRID_MPI_CHECK(MPI_Type_contiguous(sizeof(T), MPI_BYTE, &baseType), "Failed to create a contiguous data type");
   const std::array<int32_t, 3> reverseStorageSize = {
       storageSize[2],
       storageSize[1],
       storageSize[0],
   };
   std::array<MPI_Datatype, 27> types;
   types.fill(MPI_DATATYPE_NULL);

   for (uint32_t i = 0; i < 27; i++) {
      const auto x = linearToX(i);
      const auto y = linearToY(i);
      const auto z = linearToZ(i);

      const bool self = x == 0 && y == 0 && z == 0;
      const bool flatX = storageSize[0] == 1 && x != 0;
      const bool flatY = storageSize[1] == 1 && y != 0;
      const bool flatZ = storageSize[2] == 1 && z != 0;
      const bool skip = flatX || flatY || flatZ || self;

      if (skip) {
         continue;
      }

      const std::array<int32_t, 3> reverseSubarraySize = {
          (z == 0) ? localSize[2] : stencilSize,
          (y == 0) ? localSize[1] : stencilSize,
          (x == 0) ? localSize[0] : stencilSize,
      };

      const std::array<int32_t, 3> reverseSubarrayStart = [&]() {
         if (generateForSend) {
            return std::array<int32_t, 3>{
                storageSize[2] == 1 ? 0 : (z == 1 ? storageSize[2] - 2 * stencilSize : stencilSize),
                storageSize[1] == 1 ? 0 : (y == 1 ? storageSize[1] - 2 * stencilSize : stencilSize),
                storageSize[0] == 1 ? 0 : (x == 1 ? storageSize[0] - 2 * stencilSize : stencilSize),
            };
         } else {
            return std::array<int32_t, 3>{
                storageSize[2] == 1 ? 0 : (z == -1 ? storageSize[2] - stencilSize : (z == 0 ? stencilSize : 0)),
                storageSize[1] == 1 ? 0 : (y == -1 ? storageSize[1] - stencilSize : (y == 0 ? stencilSize : 0)),
                storageSize[0] == 1 ? 0 : (x == -1 ? storageSize[0] - stencilSize : (x == 0 ? stencilSize : 0)),
            };
         }
      }();

      FSGRID_MPI_CHECK(MPI_Type_create_subarray(3, reverseStorageSize.data(), reverseSubarraySize.data(),
                                                reverseSubarrayStart.data(), MPI_ORDER_C, baseType, &(types[i])),
                       "Failed to create a subarray type");
      FSGRID_MPI_CHECK(MPI_Type_commit(&(types[i])), "Failed to commit MPI type");
   }

   FSGRID_MPI_CHECK(MPI_Type_free(&baseType), "Couldn't free the basetype used to create the sendTypes");

   return types;
}
} // namespace fsgrid_detail

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

public:
   /*! Constructor for this grid.
    * \param globalSize Cell size of the global simulation domain.
    * \param MPI_Comm The MPI communicator this grid should use.
    * \param periodic An array specifying, for each dimension, whether it is to be treated as periodic.
    */
   FsGrid(const std::array<FsSize_t, 3>& globalSize, MPI_Comm parentComm, const std::array<bool, 3>& periodic,
          const std::array<double, 3>& physicalGridSpacing, const std::array<double, 3>& physicalGlobalStart,
          const std::array<Task_t, 3>& decomposition = {0, 0, 0})
       : comm3d(fsgrid_detail::createCartesianCommunicator(
             parentComm,
             fsgrid_detail::computeNumTasksPerDim(globalSize, decomposition,
                                                  fsgrid_detail::getFSCommSize(fsgrid_detail::getCommSize(parentComm)),
                                                  stencil),
             periodic)),
         rank(fsgrid_detail::getCartesianRank(parentComm, comm3d)),
         coordinates(physicalGridSpacing, physicalGlobalStart, globalSize, periodic, decomposition,
                     fsgrid_detail::getTaskPosition(comm3d),
                     fsgrid_detail::getFSCommSize(fsgrid_detail::getCommSize(parentComm)), stencil),
         neighbourIndexToRank(fsgrid_detail::mapNeigbourIndexToRank(
             fsgrid_detail::getTaskPosition(comm3d), coordinates.numTasksPerDim, periodic, comm3d, rank)),
         neighbourRankToIndex(fsgrid_detail::mapNeighbourRankToIndex(
             neighbourIndexToRank, fsgrid_detail::getFSCommSize(fsgrid_detail::getCommSize(parentComm)))),
         neighbourSendType(
             fsgrid_detail::generateMPITypes<T>(coordinates.storageSize, coordinates.localSize, stencil, true)),
         neighbourReceiveType(
             fsgrid_detail::generateMPITypes<T>(coordinates.storageSize, coordinates.localSize, stencil, false)),
         data(rank == -1
                  ? 0ul
                  : static_cast<size_t>(std::accumulate(coordinates.storageSize.cbegin(),
                                                        coordinates.storageSize.cend(), 1, std::multiplies<>()))) {}

   ~FsGrid() { finalize(); }

   void finalize() noexcept {
      // If not a non-FS process
      if (rank != -1) {
         for (size_t i = 0; i < 27; i++) {
            if (neighbourReceiveType[i] != MPI_DATATYPE_NULL)
               FSGRID_MPI_CHECK(MPI_Type_free(&(neighbourReceiveType[i])), "Failed to free MPI type");
            if (neighbourSendType[i] != MPI_DATATYPE_NULL)
               FSGRID_MPI_CHECK(MPI_Type_free(&(neighbourSendType[i])), "Failed to free MPI type");
         }
      }

      if (comm3d != MPI_COMM_NULL)
         FSGRID_MPI_CHECK(MPI_Comm_free(&comm3d), "Failed to free MPI comm3d");
   }

   // ============================
   // Data access functions
   // ============================
   std::vector<T>& getData() { return data; }
   const std::vector<T>& getData() const { return data; }

   bool localIdInBounds(LocalID id) const { return 0 <= id && (size_t)id < data.size(); }

   T* get(LocalID id) {
      FSGRID_DEBUG_ASSERT(localIdInBounds(id), "Out-of bounds access in FsGrid::get!", "(LocalID = ", id,
                          ", but storage space is ", data.size(), ". Expect weirdness.");
      return localIdInBounds(id) ? &data[static_cast<size_t>(id)] : nullptr;
   }

   /*! Get a reference to the field data in a cell
    * \param x x-Coordinate, in cells
    * \param y y-Coordinate, in cells
    * \param z z-Coordinate, in cells
    * \return A reference to cell data in the given cell
    */
   T* get(FsIndex_t x, FsIndex_t y, FsIndex_t z) { return get(localIDFromCellCoordinates(x, y, z)); }

   // ============================
   // Coordinate change functions
   // - Redirected to Coordinates' implementation
   // ============================
   template <typename... Args> auto globalIDFromLocalCoordinates(Args... args) const {
      return coordinates.globalIDFromLocalCoordinates(args...);
   }
   template <typename... Args> auto localIDFromLocalCoordinates(Args... args) const {
      return coordinates.localIDFromLocalCoordinates(args...);
   }
   template <typename... Args> auto globalToLocal(Args... args) const { return coordinates.globalToLocal(args...); }
   template <typename... Args> auto localToGlobal(Args... args) const { return coordinates.localToGlobal(args...); }
   template <typename... Args> auto getPhysicalCoords(Args... args) const {
      return coordinates.getPhysicalCoords(args...);
   }
   template <typename... Args> auto physicalToGlobal(Args... args) const {
      return coordinates.physicalToGlobal(args...);
   }
   template <typename... Args> auto physicalToFractionalGlobal(Args... args) const {
      return coordinates.physicalToFractionalGlobal(args...);
   }

   /*! Compute the local id from cell coordinates (these include ghost cells)
    * \param x x-Coordinate, in cells
    * \param y y-Coordinate, in cells
    * \param z z-Coordinate, in cells
    * \return local id of the cell
    */
   LocalID localIDFromCellCoordinates(FsIndex_t x, FsIndex_t y, FsIndex_t z) const {
      FSGRID_DEBUG_ASSERT(coordinates.cellIndicesAreWithinBounds(x, y, z), "Out-of bounds access in FsGrid::get!");
      const auto neighbourIndex = coordinates.neighbourIndexFromCellCoordinates(x, y, z);
      const auto neighbourRank = neighbourIndexToRank[neighbourIndex];
      const auto isSelf = neighbourRank == rank;

      FSGRID_DEBUG_ASSERT(isSelf || neighbourRank != MPI_PROC_NULL,
                          "Trying to access data from a non-existing neighbour");

      const auto neighbourIsSelf = neighbourIndex != 13 && isSelf;
      const auto id = neighbourIsSelf ? coordinates.localIDFromLocalCoordinates(coordinates.shiftCellIndices(x, y, z))
                                      : coordinates.localIDFromLocalCoordinates(x, y, z);

      return coordinates.cellIndicesAreWithinBounds(x, y, z) && (isSelf || neighbourRank != MPI_PROC_NULL)
                 ? id
                 : std::numeric_limits<LocalID>::min();
   }

   LocalID localIDFromCellCoordinates(const std::array<FsIndex_t, 3>& indices) const {
      return localIDFromCellCoordinates(indices[0], indices[1], indices[2]);
   }

   // ============================
   // Getters
   // ============================

   /*! Get the size of the local domain handled by this grid.
    */
   const auto& getLocalSize() const {
      constexpr static std::array zero{0, 0, 0};
      return rank == -1 ? zero : coordinates.localSize;
   }

   /*! Get the start coordinates of the local domain handled by this grid.
    */
   const auto& getLocalStart() const { return coordinates.localStart; }

   /*! Get global size of the fsgrid domain
    */
   const auto& getGlobalSize() const { return coordinates.globalSize; }

   /*! Get the rank of this CPU in the FsGrid communicator */
   Task_t getRank() const { return rank; }

   /*! Get the number of ranks in the FsGrid communicator */
   Task_t getNumFsRanks() const {
      return coordinates.numTasksPerDim[0] * coordinates.numTasksPerDim[1] * coordinates.numTasksPerDim[2];
   }

   /*! Get in which directions, if any, this grid is periodic */
   const auto& getPeriodic() const { return coordinates.periodic; }

   /*! Get the decomposition array*/
   const auto& getDecomposition() const { return coordinates.numTasksPerDim; }

   /*! Get the physical grid spacing array*/
   const auto& getGridSpacing() const { return coordinates.physicalGridSpacing; }

   // ============================
   // MPI functions
   // ============================

   /*! Perform ghost cell communication.
    */
   template <typename D> void updateGhostCells(D& data) {
      if (rank == -1) {
         return;
      }

      std::array<MPI_Request, 27> receiveRequests;
      std::array<MPI_Request, 27> sendRequests;
      receiveRequests.fill(MPI_REQUEST_NULL);
      sendRequests.fill(MPI_REQUEST_NULL);

      for (uint32_t shiftId = 0; shiftId < 27; shiftId++) {
         const auto receiveId = 26 - shiftId;
         const auto receiveFrom = neighbourIndexToRank[receiveId];
         const auto sendType = neighbourSendType[shiftId];
         const auto receiveType = neighbourReceiveType[shiftId];
         // Is this a bug? Should the check be on receiveType, not sendType? It has been like this since 2016
         if (receiveFrom != MPI_PROC_NULL && sendType != MPI_DATATYPE_NULL) {
            FSGRID_MPI_CHECK(
                MPI_Irecv(data.data(), 1, receiveType, receiveFrom, shiftId, comm3d, &(receiveRequests[shiftId])),
                "Rank ", rank, " failed to receive data from neighbor ", receiveId, " with rank ", receiveFrom);
         }
      }

      for (uint32_t shiftId = 0; shiftId < 27; shiftId++) {
         const auto sendTo = neighbourIndexToRank[shiftId];
         const auto sendType = neighbourSendType[shiftId];
         if (sendTo != MPI_PROC_NULL && sendType != MPI_DATATYPE_NULL) {
            FSGRID_MPI_CHECK(MPI_Isend(data.data(), 1, sendType, sendTo, shiftId, comm3d, &(sendRequests[shiftId])),
                             "Rank ", rank, " failed to send data to neighbor ", shiftId, " with rank ", sendTo);
         }
      }

      FSGRID_MPI_CHECK(MPI_Waitall(27, receiveRequests.data(), MPI_STATUSES_IGNORE),
                       "Synchronization at ghost cell update failed");
      FSGRID_MPI_CHECK(MPI_Waitall(27, sendRequests.data(), MPI_STATUSES_IGNORE),
                       "Synchronization at ghost cell update failed");
   }

   void updateGhostCells() { updateGhostCells(data); }

   /*! Perform an MPI_Allreduce with this grid's internal communicator
    * Function syntax is identical to MPI_Allreduce, except the final (communicator
    * argument will not be needed) */
   int32_t Allreduce(void* sendbuf, void* recvbuf, int32_t count, MPI_Datatype datatype, MPI_Op op) const {
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

   /*! Returns the task responsible for handling the cell with the given GlobalID
    * \param id GlobalID of the cell for which task is to be determined
    * \return a task for the grid's cartesian communicator
    */
   Task_t getTaskForGlobalID(GlobalID id) const {
      const auto taskPos = coordinates.globalIdToTaskPos(id);
      Task_t taskID = -1;
      FSGRID_MPI_CHECK(MPI_Cart_rank(comm3d, taskPos.data(), &taskID), "Unable to find FsGrid rank for global ID ", id);

      return taskID;
   }
private:
   //! MPI Cartesian communicator used in this grid
   MPI_Comm comm3d = MPI_COMM_NULL;
   //!< This task's rank in the communicator
   const int32_t rank = 0;

   //!< A container for the coordinates of the fsgrid
   const Coordinates coordinates = {};

   //!< Lookup table from index to rank in the neighbour array (includes self)
   const std::array<int32_t, 27> neighbourIndexToRank = {
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
   };
   //!< Lookup table from rank to index in the neighbour array
   const std::vector<char> neighbourRankToIndex = {};
   //!< Datatype for sending data
   std::array<MPI_Datatype, 27> neighbourSendType = {};
   //!< Datatype for receiving data
   std::array<MPI_Datatype, 27> neighbourReceiveType = {};

   //! Actual storage of field data
   std::vector<T> data = {};
};
