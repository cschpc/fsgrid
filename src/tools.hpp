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
#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <stdexcept>
#include <stdint.h>
#include <vector>

#define FSGRID_MPI_CHECK(status, ...) FsGridTools::writeToCerrAndThrowIfFailed(status != MPI_SUCCESS, __VA_ARGS__)

namespace FsGridTools {
// Size type for global array indices
typedef uint32_t FsSize_t;
// Size type for global/local array indices, incl. possible negative values
typedef int32_t FsIndex_t;
typedef int64_t LocalID;
typedef int64_t GlobalID;
typedef int32_t Task_t;

//! Helper function: calculate size of the local coordinate space for the given dimension
// \param numCells Number of cells in the global Simulation, in this dimension
// \param numTasks Total number of tasks in this dimension
// \param taskIndex This task's position in this dimension
// \return Number of cells for this task's local domain (actual cells, not counting ghost cells)
static FsIndex_t calcLocalSize(FsSize_t numCells, Task_t numTasks, Task_t taskIndex) {
   const FsIndex_t nPerTask = numCells / numTasks;
   const FsIndex_t remainder = numCells % numTasks;
   return nPerTask + (taskIndex < remainder);
}

//! Helper function: calculate position of the local coordinate space for the given dimension
// \param numCells number of cells
// \param numTasks Total number of tasks in this dimension
// \param taskIndex This task's position in this dimension
// \return Cell number at which this task's domains cells start (actual cells, not counting ghost cells)
static FsIndex_t calcLocalStart(FsSize_t numCells, Task_t numTasks, Task_t taskIndex) {
   const FsIndex_t remainder = numCells % numTasks;
   return taskIndex * calcLocalSize(numCells, numTasks, taskIndex) + (taskIndex >= remainder) * remainder;
}

//! Helper function: given a global cellID, calculate the global cell coordinate from it.
// This is then used do determine the task responsible for this cell, and the
// local cell index in it.
static std::array<FsIndex_t, 3> globalIDtoCellCoord(GlobalID id, const std::array<FsSize_t, 3>& globalSize) {
   // We're returning FsIndex_t, which is int32_t. globalSizes are FsSize_t, which are uint32_t.
   // Need to check that the globalSizes aren't larger than maximum of int32_t
   const bool globalSizeOverflows = globalSize[0] > std::numeric_limits<FsIndex_t>::max() ||
                                    globalSize[1] > std::numeric_limits<FsIndex_t>::max() ||
                                    globalSize[2] > std::numeric_limits<FsIndex_t>::max();
   const bool idNegative = id < 0;

   // To avoid overflow, do this instead of checking for id < product of globalSize
   // The number of divisions stays the same anyway
   const GlobalID idPerGs0 = id / globalSize[0];
   const GlobalID idPerGs0PerGs1 = idPerGs0 / globalSize[1];
   const bool idTooLarge = idPerGs0PerGs1 >= globalSize[2];

   const bool badInput = idTooLarge || idNegative || globalSizeOverflows;

   if (badInput) {
      // For bad input, return bad output
      return {
          std::numeric_limits<FsIndex_t>::min(),
          std::numeric_limits<FsIndex_t>::min(),
          std::numeric_limits<FsIndex_t>::min(),
      };
   } else {
      return {
          (FsIndex_t)(id % globalSize[0]),
          (FsIndex_t)(idPerGs0 % globalSize[1]),
          (FsIndex_t)(idPerGs0PerGs1 % globalSize[2]),
      };
   }
}

//! Helper function to optimize decomposition of this grid over the given number of tasks
static std::array<Task_t, 3> computeDomainDecomposition(const std::array<FsSize_t, 3>& globalSize, Task_t nProcs,
                                                        int32_t stencilSize = 1) {
   const std::array minDomainSize = {
       globalSize[0] == 1 ? 1 : stencilSize,
       globalSize[1] == 1 ? 1 : stencilSize,
       globalSize[2] == 1 ? 1 : stencilSize,
   };
   const std::array maxDomainSize = {
       std::min(nProcs, static_cast<Task_t>(globalSize[0] / minDomainSize[0])),
       std::min(nProcs, static_cast<Task_t>(globalSize[1] / minDomainSize[1])),
       std::min(nProcs, static_cast<Task_t>(globalSize[2] / minDomainSize[2])),
   };

   int64_t minimumCost = std::numeric_limits<int64_t>::max();
   std::array dd = {1, 1, 1};
   for (Task_t i = 1; i <= maxDomainSize[0]; i++) {
      for (Task_t j = 1; j <= maxDomainSize[1]; j++) {
         const Task_t k = nProcs / (i * j);
         if (k == 0) {
            break;
         }

         // No need to optimize an incompatible DD, also checks for missing remainders
         if (i * j * k != nProcs || k > maxDomainSize[2]) {
            continue;
         }

         const std::array processBox = {
             calcLocalSize(globalSize[0], i, 0),
             calcLocalSize(globalSize[1], j, 0),
             calcLocalSize(globalSize[2], k, 0),
         };

         // clang-format off
         const int64_t baseCost = (i > 1 ? processBox[1] * processBox[2] : 0)
                                + (j > 1 ? processBox[0] * processBox[2] : 0)
                                + (k > 1 ? processBox[0] * processBox[1] : 0);
         const int64_t neighborMultiplier = (i != 1 && j != 1 && k != 1) * 13
                                          + (i == 1 && j != 1 && k != 1) * 4
                                          + (i != 1 && j == 1 && k != 1) * 4
                                          + (i != 1 && j != 1 && k == 1) * 4;
         // clang-format on
         const int64_t cost = baseCost * neighborMultiplier;
         if (cost < minimumCost) {
            minimumCost = cost;
            dd = {i, j, k};
         }
      }
   }

   if (minimumCost == std::numeric_limits<int64_t>::max() || (Task_t)(dd[0] * dd[1] * dd[2]) != nProcs) {
      throw std::runtime_error("FSGrid computeDomainDecomposition failed");
   }

   return dd;
}

template <typename... Args> void cerrArgs(Args...);
template <> inline void cerrArgs() {}
template <typename T> void cerrArgs(T t) { std::cerr << t << "\n"; }
template <typename Head, typename... Tail> void cerrArgs(Head head, Tail... tail) {
   cerrArgs(head);
   cerrArgs(tail...);
}

// Recursively write all arguments to cerr if status is not success, then throw
template <typename... Args> void writeToCerrAndThrowIfFailed(bool failed, Args... args) {
   if (failed) {
      cerrArgs(args...);
      throw std::runtime_error("Unrecoverable error encountered in FsGrid, consult cerr for more information");
   }
}

static int32_t getCommRank(MPI_Comm comm) {
   int32_t rank;
   FSGRID_MPI_CHECK(MPI_Comm_rank(comm, &rank), "Couldn't get rank from communicator");
   return rank;
}

static int32_t getCommSize(MPI_Comm comm) {
   int32_t size;
   FSGRID_MPI_CHECK(MPI_Comm_size(comm, &size), "Couldn't get size from communicator");
   return size;
}

static int32_t getNumFsGridProcs(int32_t parentCommSize) {
   const auto envVar = getenv("FSGRID_PROCS");
   const auto fsgridProcs = envVar != NULL ? atoi(envVar) : 0;
   return parentCommSize > fsgridProcs && fsgridProcs > 0 ? fsgridProcs : parentCommSize;
}

constexpr static int32_t computeColourFs(int32_t parentRank, int32_t numRanks) {
   return (parentRank < numRanks) ? 1 : MPI_UNDEFINED;
}

constexpr static int32_t computeColourAux(int32_t parentRank, int32_t parentCommSize, int32_t numRanks) {
   return (parentRank > (parentCommSize - 1) % numRanks) ? (parentRank - (parentCommSize % numRanks)) / numRanks
                                                         : MPI_UNDEFINED;
}

static std::array<Task_t, 3> computeNumTasksPerDim(std::array<FsSize_t, 3> globalSize,
                                                   const std::array<Task_t, 3>& decomposition, int32_t numRanks,
                                                   int32_t stencilSize) {
   const bool allZero = decomposition[0] == 0 && decomposition[1] == 0 && decomposition[2] == 0;
   if (allZero) {
      return computeDomainDecomposition(globalSize, numRanks, stencilSize);
   }

   const bool incorrectDistribution = decomposition[0] * decomposition[1] * decomposition[2] != numRanks;
   if (incorrectDistribution) {
      std::cerr << "Given decomposition (" << decomposition[0] << " " << decomposition[1] << " " << decomposition[2]
                << ") does not distribute to the number of tasks given" << std::endl;
      throw std::runtime_error("Given decomposition does not distribute to the number of tasks given");
   }

   return decomposition;
}

static MPI_Comm createCartesianCommunicator(MPI_Comm parentComm, int32_t colourFs, int32_t colourAux,
                                            int32_t parentRank, const std::array<Task_t, 3>& numTasksPerDim,
                                            const std::array<bool, 3>& isPeriodic) {
   MPI_Comm comm = MPI_COMM_NULL;
   if (colourFs != MPI_UNDEFINED) {
      FSGRID_MPI_CHECK(MPI_Comm_split(parentComm, colourFs, parentRank, &comm),
                       "Couldn's split parent communicator to subcommunicators");
   } else {
      FSGRID_MPI_CHECK(MPI_Comm_split(parentComm, colourAux, parentRank, &comm),
                       "Couldn's split parent communicator to subcommunicators");
   }
   const std::array<int32_t, 3> pi = {
       isPeriodic[0],
       isPeriodic[1],
       isPeriodic[2],
   };
   MPI_Comm comm3d = MPI_COMM_NULL;
   FSGRID_MPI_CHECK(MPI_Cart_create(comm, 3, numTasksPerDim.data(), pi.data(), 0, &comm3d),
                    "Creating cartesian communicatior failed when attempting to create FsGrid!");
   FSGRID_MPI_CHECK(MPI_Comm_free(&comm), "Failed to free MPI comm");

   return comm3d;
}

static int32_t getCartesianRank(int32_t colourFs, MPI_Comm comm) {
   return colourFs != MPI_UNDEFINED ? getCommRank(comm) : -1;
}

static std::array<Task_t, 3> getTaskPosition(MPI_Comm comm) {
   std::array<Task_t, 3> taskPos;
   const int32_t rank = getCommRank(comm);
   FSGRID_MPI_CHECK(MPI_Cart_coords(comm, rank, 3, taskPos.data()), "Rank ", rank,
                    " unable to determine own position in cartesian communicator when attempting to create FsGrid!");
   return taskPos;
}

constexpr static bool localSizeTooSmall(std::array<FsSize_t, 3> globalSize, std::array<FsIndex_t, 3> localSize,
                                        int32_t stencilSize) {
   const bool anyLocalIsZero = localSize[0] == 0 || localSize[1] == 0 || localSize[2] == 0;
   const bool stencilSizeBoundedByGlobalAndLocalSizes = (globalSize[0] > stencilSize && stencilSize > localSize[0]) ||
                                                        (globalSize[1] > stencilSize && stencilSize > localSize[1]) ||
                                                        (globalSize[2] > stencilSize && stencilSize > localSize[2]);

   return anyLocalIsZero || stencilSizeBoundedByGlobalAndLocalSizes;
}

static std::array<FsIndex_t, 3> calculateLocalSize(const std::array<FsSize_t, 3>& globalSize,
                                                   const std::array<Task_t, 3>& numTasksPerDim,
                                                   const std::array<Task_t, 3>& taskPosition, int32_t rank,
                                                   int32_t stencilSize) {
   std::array localSize = {
       calcLocalSize(globalSize[0], numTasksPerDim[0], taskPosition[0]),
       calcLocalSize(globalSize[1], numTasksPerDim[1], taskPosition[1]),
       calcLocalSize(globalSize[2], numTasksPerDim[2], taskPosition[2]),
   };

   if (localSizeTooSmall(globalSize, localSize, stencilSize)) {
      std::cerr << "FSGrid space partitioning leads to a space that is too small on Rank " << rank << "." << std::endl;
      std::cerr << "Please run with a different number of Tasks, so that space is better divisible." << std::endl;
      throw std::runtime_error("FSGrid too small domains");
   }

   return rank == -1 ? std::array{0, 0, 0} : localSize;
}

static std::array<FsIndex_t, 3> calculateLocalStart(const std::array<FsSize_t, 3>& globalSize,
                                                    const std::array<Task_t, 3>& numTasksPerDim,
                                                    const std::array<Task_t, 3>& taskPosition) {
   return {
       calcLocalStart(globalSize[0], numTasksPerDim[0], taskPosition[0]),
       calcLocalStart(globalSize[1], numTasksPerDim[1], taskPosition[1]),
       calcLocalStart(globalSize[2], numTasksPerDim[2], taskPosition[2]),
   };
}

static std::array<FsIndex_t, 3> calculateStorageSize(const std::array<FsSize_t, 3>& globalSize,
                                                     const std::array<Task_t, 3>& localSize, int32_t stencilSize) {
   return {
       globalSize[0] <= 1 ? 1 : localSize[0] + stencilSize * 2,
       globalSize[1] <= 1 ? 1 : localSize[1] + stencilSize * 2,
       globalSize[2] <= 1 ? 1 : localSize[2] + stencilSize * 2,
   };
}

// Assumes x, y and z to belong to set [-1, 0, 1]
// returns a value in (inclusive) range [0, 26]
constexpr static int32_t xyzToLinear(int32_t x, int32_t y, int32_t z) { return (x + 1) * 9 + (y + 1) * 3 + (z + 1); }

// These assume i to be in (inclusive) range [0, 26]
// returns a value from the set [-1, 0, 1]
constexpr static int32_t linearToX(int32_t i) { return i / 9 - 1; }
constexpr static int32_t linearToY(int32_t i) { return (i % 9) / 3 - 1; }
constexpr static int32_t linearToZ(int32_t i) { return i % 3 - 1; }

static std::array<int32_t, 27> mapNeigbourIndexToRank(const std::array<Task_t, 3>& taskPosition,
                                                      const std::array<Task_t, 3>& numTasksPerDim,
                                                      const std::array<bool, 3>& periodic, MPI_Comm comm,
                                                      int32_t rank) {
   auto calculateNeighbourRank = [&](int32_t neighbourIndex) {
      auto calculateNeighbourPosition = [&](int32_t neighbourIndex, uint32_t i) {
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
                    [&calculateNeighbourRank, n = 0]() mutable { return calculateNeighbourRank(n++); });
   }
   return ranks;
}

static std::vector<char> mapNeighbourRankToIndex(const std::array<int32_t, 27>& indexToRankMap, FsSize_t numRanks) {
   std::vector<char> indices(numRanks, MPI_PROC_NULL);
   std::for_each(indexToRankMap.cbegin(), indexToRankMap.cend(), [&indices, &numRanks, n = 0](auto rank) mutable {
      if (numRanks > rank && rank >= 0) {
         indices[rank] = n;
      }
      n++;
   });
   return indices;
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

   for (int32_t i = 0; i < 27; i++) {
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

} // namespace FsGridTools
