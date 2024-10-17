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

namespace fsgrid_detail {
using FsSize_t = FsGridTools::FsSize_t;
using FsIndex_t = FsGridTools::FsIndex_t;
using LocalID = FsGridTools::LocalID;
using GlobalID = FsGridTools::GlobalID;
using Task_t = FsGridTools::Task_t;

static std::array<Task_t, 3> computeNumTasksPerDim(std::array<FsSize_t, 3> globalSize,
                                                   const std::array<Task_t, 3>& decomposition, int32_t numRanks,
                                                   int32_t numGhostCells) {
   const bool allZero = decomposition[0] == 0 && decomposition[1] == 0 && decomposition[2] == 0;
   if (allZero) {
      return FsGridTools::computeDomainDecomposition(globalSize, numRanks, numGhostCells);
   }

   const bool incorrectDistribution = decomposition[0] * decomposition[1] * decomposition[2] != numRanks;
   if (incorrectDistribution) {
      std::cerr << "Given decomposition (" << decomposition[0] << " " << decomposition[1] << " " << decomposition[2]
                << ") does not distribute to the number of tasks given" << std::endl;
      throw std::runtime_error("Given decomposition does not distribute to the number of tasks given");
   }

   return decomposition;
}

constexpr static bool localSizeTooSmall(std::array<FsSize_t, 3> globalSize, std::array<FsIndex_t, 3> localSize,
                                        int32_t numGhostCells) {
   const bool anyLocalIsZero = localSize[0] == 0 || localSize[1] == 0 || localSize[2] == 0;
   const bool bounded = (globalSize[0] > static_cast<uint32_t>(numGhostCells) && numGhostCells > localSize[0]) ||
                        (globalSize[1] > static_cast<uint32_t>(numGhostCells) && numGhostCells > localSize[1]) ||
                        (globalSize[2] > static_cast<uint32_t>(numGhostCells) && numGhostCells > localSize[2]);

   return anyLocalIsZero || bounded;
}

static std::array<FsIndex_t, 3> calculateLocalSize(const std::array<FsSize_t, 3>& globalSize,
                                                   const std::array<Task_t, 3>& numTasksPerDim,
                                                   const std::array<Task_t, 3>& taskPosition, int rank,
                                                   int32_t numGhostCells) {
   std::array localSize = {
       FsGridTools::calcLocalSize(globalSize[0], numTasksPerDim[0], taskPosition[0]),
       FsGridTools::calcLocalSize(globalSize[1], numTasksPerDim[1], taskPosition[1]),
       FsGridTools::calcLocalSize(globalSize[2], numTasksPerDim[2], taskPosition[2]),
   };

   if (localSizeTooSmall(globalSize, localSize, numGhostCells)) {
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
       FsGridTools::calcLocalStart(globalSize[0], numTasksPerDim[0], taskPosition[0]),
       FsGridTools::calcLocalStart(globalSize[1], numTasksPerDim[1], taskPosition[1]),
       FsGridTools::calcLocalStart(globalSize[2], numTasksPerDim[2], taskPosition[2]),
   };
}

static std::array<FsIndex_t, 3> calculateStorageSize(const std::array<FsSize_t, 3>& globalSize,
                                                     const std::array<Task_t, 3>& localSize, int32_t numGhostCells) {
   return {
       globalSize[0] <= 1 ? 1 : localSize[0] + numGhostCells * 2,
       globalSize[1] <= 1 ? 1 : localSize[1] + numGhostCells * 2,
       globalSize[2] <= 1 ? 1 : localSize[2] + numGhostCells * 2,
   };
}
} // namespace fsgrid_detail

struct Coordinates {
private:
   using FsSize_t = FsGridTools::FsSize_t;
   using FsIndex_t = FsGridTools::FsIndex_t;
   using LocalID = FsGridTools::LocalID;
   using GlobalID = FsGridTools::GlobalID;
   using Task_t = FsGridTools::Task_t;

public:
   Coordinates() {}
   Coordinates(const std::array<double, 3>& physicalGridSpacing, const std::array<double, 3>& physicalGlobalStart,
               std::array<FsSize_t, 3> globalSize, std::array<bool, 3> periodic,
               const std::array<Task_t, 3>& decomposition, const std::array<Task_t, 3>& taskPosition, int32_t numRanks,
               Task_t rank, int32_t numGhostCells)
       : physicalGridSpacing(physicalGridSpacing), physicalGlobalStart(physicalGlobalStart), globalSize(globalSize),
         periodic(periodic),
         numTasksPerDim(fsgrid_detail::computeNumTasksPerDim(globalSize, decomposition, numRanks, numGhostCells)),
         localSize(fsgrid_detail::calculateLocalSize(globalSize, numTasksPerDim, taskPosition, rank, numGhostCells)),
         localStart(fsgrid_detail::calculateLocalStart(globalSize, numTasksPerDim, taskPosition)),
         storageSize(fsgrid_detail::calculateStorageSize(globalSize, localSize, numGhostCells)) {}

   const std::array<double, 3> physicalGridSpacing = {};
   const std::array<double, 3> physicalGlobalStart = {};
   const std::array<FsSize_t, 3> globalSize = {};
   const std::array<bool, 3> periodic = {};
   const std::array<Task_t, 3> numTasksPerDim = {};
   const std::array<FsIndex_t, 3> localSize = {};
   const std::array<FsIndex_t, 3> localStart = {};
   const std::array<FsIndex_t, 3> storageSize = {};
};
