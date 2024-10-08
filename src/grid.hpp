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
#include "fsstencil.hpp"
#include "tools.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <vector>

namespace fsgrid_detail {
using namespace fsgrid;

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
         mpiCheck(MPI_Cart_rank(comm, neighbourPosition.data(), &neighbourRank), "Rank ", rank,
                  " can't determine neighbour rank at position [", neighbourPosition[0], ", ", neighbourPosition[1],
                  ", ", neighbourPosition[2], "]");
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

static int32_t getCommRank(MPI_Comm parentComm) {
   int32_t parentRank = -1;
   mpiCheck(MPI_Comm_rank(parentComm, &parentRank), "Couldn't get rank from parent communicator");
   return parentRank;
}

static MPI_Comm createCartesianCommunicator(MPI_Comm parentComm, const std::array<Task_t, 3>& numTasksPerDim,
                                            const std::array<bool, 3>& isPeriodic, int32_t numProcs) {
   const auto parentRank = getCommRank(parentComm);
   const auto colour = (parentRank < numProcs) ? 1 : MPI_UNDEFINED;

   MPI_Comm comm = MPI_COMM_NULL;
   mpiCheck(MPI_Comm_split(parentComm, colour, parentRank, &comm),
            "Couldn's split parent communicator to subcommunicators");

   const std::array<int32_t, 3> pi = {
       isPeriodic[0],
       isPeriodic[1],
       isPeriodic[2],
   };

   MPI_Comm comm3d = MPI_COMM_NULL;
   if (comm != MPI_COMM_NULL) {
      mpiCheck(MPI_Cart_create(comm, 3, numTasksPerDim.data(), pi.data(), 0, &comm3d),
               "Creating cartesian communicatior failed when attempting to create FsGrid!");

      mpiCheck(MPI_Comm_free(&comm), "Failed to free MPI comm");
   }

   return comm3d;
}

static int32_t getCartesianRank(MPI_Comm cartesianComm) {
   return cartesianComm != MPI_COMM_NULL ? getCommRank(cartesianComm) : -1;
}

static std::array<int32_t, 3> getTaskPosition(MPI_Comm comm) {
   std::array<int32_t, 3> taskPos{-1, -1, -1};
   if (comm != MPI_COMM_NULL) {
      const int rank = getCommRank(comm);
      mpiCheck(MPI_Cart_coords(comm, rank, taskPos.size(), taskPos.data()), "Rank ", rank,
               " unable to determine own position in cartesian communicator when attempting to create FsGrid!");
   }
   return taskPos;
}

template <typename T>
static std::array<MPI_Datatype, 27> generateMPITypes(const std::array<FsIndex_t, 3>& storageSize,
                                                     const std::array<FsIndex_t, 3>& localSize, int32_t stencilSize,
                                                     bool generateForSend) {
   MPI_Datatype baseType;
   mpiCheck(MPI_Type_contiguous(sizeof(T), MPI_BYTE, &baseType), "Failed to create a contiguous data type");
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

      mpiCheck(MPI_Type_create_subarray(3, reverseStorageSize.data(), reverseSubarraySize.data(),
                                        reverseSubarrayStart.data(), MPI_ORDER_C, baseType, &(types[i])),
               "Failed to create a subarray type");
      mpiCheck(MPI_Type_commit(&(types[i])), "Failed to commit MPI type");
   }

   mpiCheck(MPI_Type_free(&baseType), "Couldn't free the basetype used to create the sendTypes");

   return types;
}

static std::vector<int32_t> taskPosToTask(MPI_Comm parentComm, MPI_Comm cartesianComm,
                                          const std::array<Task_t, 3>& numTasksPerDim) {
   std::vector<int32_t> tasks(static_cast<size_t>(numTasksPerDim[0] * numTasksPerDim[1] * numTasksPerDim[2]));
   if (cartesianComm != MPI_COMM_NULL) {
      size_t i = 0;
      for (auto x = 0; x < numTasksPerDim[0]; x++) {
         for (auto y = 0; y < numTasksPerDim[1]; y++) {
            for (auto z = 0; z < numTasksPerDim[2]; z++) {
               const std::array coords = {x, y, z};
               mpiCheck(MPI_Cart_rank(cartesianComm, coords.data(), &tasks[i++]),
                        "Unable to get rank from cartesian communicator");
            }
         }
      }
   }

   mpiCheck(MPI_Bcast(static_cast<void*>(tasks.data()), static_cast<int32_t>(tasks.size()), MPI_INT, 0, parentComm),
            "Unable to broadcast task pos array");

   return tasks;
}

static BitMask32 makeNeigbourBitMask(int32_t rank, const std::array<int32_t, 27>& neighbourIndexToRank) {
   auto getNeighbourBit = [&rank, &neighbourIndexToRank](uint32_t neighbourIndex) {
      const auto neighbourRank = neighbourIndexToRank[neighbourIndex];
      const auto neighbourIsSelf = neighbourRank == rank;
      return static_cast<uint32_t>(neighbourIsSelf) << neighbourIndex;
   };

   // Self is index 13, leave that bit to zero
   uint32_t bits = 0;
   for (auto i = 0u; i < 13u; i++) {
      bits |= getNeighbourBit(i);
   }

   for (auto i = 14u; i < 27u; i++) {
      bits |= getNeighbourBit(i);
   }

   return BitMask32(bits);
}

static BitMask32 makeNeigbourIsNullBitMask(const std::array<int32_t, 27>& neighbourIndexToRank) {
   auto getNeighbourBit = [&neighbourIndexToRank](uint32_t neighbourIndex) {
      const auto neighbourRank = neighbourIndexToRank[neighbourIndex];
      return static_cast<uint32_t>(neighbourRank == MPI_PROC_NULL) << neighbourIndex;
   };

   // Self is index 13, leave that bit to zero
   uint32_t bits = 0;
   for (auto i = 0u; i < 13u; i++) {
      bits |= getNeighbourBit(i);
   }

   for (auto i = 14u; i < 27u; i++) {
      bits |= getNeighbourBit(i);
   }

   return BitMask32(bits);
}

static std::array<int32_t, 3> computeStencilMultipliers(const Coordinates& coordinates) {
   return {
       static_cast<int32_t>(coordinates.globalSize[0] > 1),
       static_cast<int32_t>(coordinates.globalSize[1] > 1) * coordinates.storageSize[0],
       static_cast<int32_t>(coordinates.globalSize[2] > 1) * coordinates.storageSize[0] * coordinates.storageSize[1],
   };
}

static int32_t computeStencilOffset(const Coordinates& coordinates) {
   const auto muls = computeStencilMultipliers(coordinates);
   return coordinates.numGhostCells * (muls[0] + muls[1] + muls[2]);
}
} // namespace fsgrid_detail

/*! Simple cartesian, non-loadbalancing MPI Grid for use with the fieldsolver
 *
 * \param T datastructure containing the field in each cell which this grid manages
 * \param stencil ghost cell width of this grid
 */
namespace fsgrid {
using namespace fsgrid_detail;

template <typename T, int32_t stencil> class FsGrid {
public:
   /*! Constructor for this grid.
    * \param globalSize Cell size of the global simulation domain.
    * \param MPI_Comm The MPI communicator this grid should use.
    * \param periodic An array specifying, for each dimension, whether it is to be treated as periodic.
    */
   FsGrid(const std::array<FsSize_t, 3>& globalSize, MPI_Comm parentComm, int32_t numProcs,
          const std::array<bool, 3>& periodic, const std::array<double, 3>& physicalGridSpacing,
          const std::array<double, 3>& physicalGlobalStart, const std::array<Task_t, 3>& decomposition = {0, 0, 0})
       : numProcs(numProcs),
         comm3d(createCartesianCommunicator(
             parentComm, computeNumTasksPerDim(globalSize, decomposition, numProcs, stencil), periodic, numProcs)),
         rank(getCartesianRank(comm3d)), coordinates(physicalGridSpacing, physicalGlobalStart, globalSize, periodic,
                                                     decomposition, getTaskPosition(comm3d), numProcs, stencil),
         tasks(taskPosToTask(parentComm, comm3d, coordinates.numTasksPerDim)),
         neighbourIndexToRank(
             mapNeigbourIndexToRank(getTaskPosition(comm3d), coordinates.numTasksPerDim, periodic, comm3d, rank)),
         neighbourRankToIndex(mapNeighbourRankToIndex(neighbourIndexToRank, numProcs)),
         stencilConstants(coordinates.localSize, computeStencilMultipliers(coordinates),
                          computeStencilOffset(coordinates), stencil, makeNeigbourBitMask(rank, neighbourIndexToRank),
                          makeNeigbourIsNullBitMask(neighbourIndexToRank)),
         neighbourSendType(generateMPITypes<T>(coordinates.storageSize, coordinates.localSize, stencil, true)),
         neighbourReceiveType(generateMPITypes<T>(coordinates.storageSize, coordinates.localSize, stencil, false)),
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
               mpiCheck(MPI_Type_free(&(neighbourReceiveType[i])), "Failed to free MPI type");
            if (neighbourSendType[i] != MPI_DATATYPE_NULL)
               mpiCheck(MPI_Type_free(&(neighbourSendType[i])), "Failed to free MPI type");
         }
      }

      if (comm3d != MPI_COMM_NULL)
         mpiCheck(MPI_Comm_free(&comm3d), "Failed to free MPI comm3d");
   }

   // ============================
   // Data access functions
   // ============================
   std::vector<T>& getData() { return data; }
   const std::vector<T>& getData() const { return data; }

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

   /*! Returns the task responsible for handling the cell with the given GlobalID
    * \param id GlobalID of the cell for which task is to be determined
    * \return a task for the grid's cartesian communicator
    */
   Task_t getTaskForGlobalID(GlobalID id) const {
      const auto taskPos = coordinates.globalIdToTaskPos(id);
      const int32_t i = taskPos[0] * (coordinates.numTasksPerDim[1] * coordinates.numTasksPerDim[2]) +
                        taskPos[1] * coordinates.numTasksPerDim[2] + taskPos[2];
      return tasks[static_cast<size_t>(i)];
   }

   FsStencil makeStencil(int32_t x, int32_t y, int32_t z) const { return FsStencil(x, y, z, stencilConstants); }

   // ============================
   // Getters
   // ============================
   const auto& getLocalSize() const { return coordinates.localSize; }
   const auto& getLocalStart() const { return coordinates.localStart; }
   const auto& getGlobalSize() const { return coordinates.globalSize; }
   Task_t getRank() const { return rank; }
   Task_t getNumFsRanks() const { return numProcs; }
   const auto& getPeriodic() const { return coordinates.periodic; }
   const auto& getDecomposition() const { return coordinates.numTasksPerDim; }
   const auto& getGridSpacing() const { return coordinates.physicalGridSpacing; }

   // ============================
   // MPI functions
   // ============================

   // TODO: make a version that takes in a structure containing:
   // - span of data
   // - MPI send type
   // - MPI receive type
   /*! Perform ghost cell communication.
    */
   template <typename D> void updateGhostCells(D& data) {
      if (comm3d == MPI_COMM_NULL) {
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
            mpiCheck(MPI_Irecv(data.data(), 1, receiveType, receiveFrom, shiftId, comm3d, &(receiveRequests[shiftId])),
                     "Rank ", rank, " failed to receive data from neighbor ", receiveId, " with rank ", receiveFrom);
         }
      }

      for (uint32_t shiftId = 0; shiftId < 27; shiftId++) {
         const auto sendTo = neighbourIndexToRank[shiftId];
         const auto sendType = neighbourSendType[shiftId];
         if (sendTo != MPI_PROC_NULL && sendType != MPI_DATATYPE_NULL) {
            mpiCheck(MPI_Isend(data.data(), 1, sendType, sendTo, shiftId, comm3d, &(sendRequests[shiftId])), "Rank ",
                     rank, " failed to send data to neighbor ", shiftId, " with rank ", sendTo);
         }
      }

      mpiCheck(MPI_Waitall(27, receiveRequests.data(), MPI_STATUSES_IGNORE),
               "Synchronization at ghost cell update failed");
      mpiCheck(MPI_Waitall(27, sendRequests.data(), MPI_STATUSES_IGNORE),
               "Synchronization at ghost cell update failed");
   }

   void updateGhostCells() { updateGhostCells(data); }

   /*! Perform an MPI_Allreduce with this grid's internal communicator
    * Function syntax is identical to MPI_Allreduce, except the final (communicator
    * argument will not be needed) */
   int32_t Allreduce(void* sendbuf, void* recvbuf, int32_t count, MPI_Datatype datatype, MPI_Op op) const {
      // If a normal FS-rank
      if (comm3d != MPI_COMM_NULL) {
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

   template <typename Lambda>
   void parallel_for(Lambda loop_body) {
      // Using raw pointer for gridDims;
      // Workaround intel compiler bug in collapsed openmp loops
      // see https://github.com/fmihpc/vlasiator/commit/604c81142729c5025a0073cd5dc64a24882f1675
      const FsIndex_t* gridDims = &coordinates.localSize[0];

      #pragma omp parallel
      {
         #pragma omp for collapse(2)
         for (FsIndex_t k=0; k<gridDims[2]; k++) {
            for (FsIndex_t j=0; j<gridDims[1]; j++) {
               for (FsIndex_t i=0; i<gridDims[0]; i++) {
                  const auto s = makeStencil(i, j, k);
                  auto tech = getData()[s.center()];
                  auto sysBoundaryFlag = tech.sysBoundaryFlag;
                  auto sysBoundaryLayer = tech.sysBoundaryLayer;
                  loop_body(s, sysBoundaryFlag, sysBoundaryLayer);
               }
            }
         }
      }
   }
private:
   //! How many fieldsolver processes there are
   const int32_t numProcs = 0;
   //! MPI Cartesian communicator used in this grid
   MPI_Comm comm3d = MPI_COMM_NULL;
   //!< This task's rank in the communicator
   const int32_t rank = 0;

   //!< A container for the coordinates of the fsgrid
   const Coordinates coordinates = {};

   //!< Task position to task mapping
   const std::vector<int32_t> tasks = {};

   //!< Lookup table from index to rank in the neighbour array (includes self)
   const std::array<int32_t, 27> neighbourIndexToRank = {
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
   };
   //!< Lookup table from rank to index in the neighbour array
   const std::vector<char> neighbourRankToIndex = {};

   //!< Type containing data computed from FsGrid values that are constant for all stencils
   const StencilConstants stencilConstants = {};

   //!< Datatype for sending data
   std::array<MPI_Datatype, 27> neighbourSendType = {};
   //!< Datatype for receiving data
   std::array<MPI_Datatype, 27> neighbourReceiveType = {};

   //! Actual storage of field data
   std::vector<T> data = {};
};
} // namespace fsgrid
