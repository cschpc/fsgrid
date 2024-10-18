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

#include <cmath>

namespace fsgrid_detail {
using FsSize_t = FsGridTools::FsSize_t;
using FsIndex_t = FsGridTools::FsIndex_t;
using LocalID = FsGridTools::LocalID;
using GlobalID = FsGridTools::GlobalID;
using Task_t = FsGridTools::Task_t;

constexpr static std::array<Task_t, 3> computeNumTasksPerDim(std::array<FsSize_t, 3> globalSize,
                                                             const std::array<Task_t, 3>& decomposition,
                                                             int32_t numRanks, int32_t numGhostCells) {
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

constexpr static std::array<FsIndex_t, 3> calculateLocalSize(const std::array<FsSize_t, 3>& globalSize,
                                                             const std::array<Task_t, 3>& numTasksPerDim,
                                                             const std::array<Task_t, 3>& taskPosition,
                                                             int32_t numGhostCells) {
   std::array localSize = {
       FsGridTools::calcLocalSize(globalSize[0], numTasksPerDim[0], taskPosition[0]),
       FsGridTools::calcLocalSize(globalSize[1], numTasksPerDim[1], taskPosition[1]),
       FsGridTools::calcLocalSize(globalSize[2], numTasksPerDim[2], taskPosition[2]),
   };

   if (localSizeTooSmall(globalSize, localSize, numGhostCells)) {
      std::cerr << "FSGrid space partitioning leads to a space that is too small.\n";
      std::cerr << "Please run with a different number of Tasks, so that space is better divisible." << std::endl;
      throw std::runtime_error("FSGrid too small domains");
   }

   return localSize;
}

constexpr static std::array<FsIndex_t, 3> calculateLocalStart(const std::array<FsSize_t, 3>& globalSize,
                                                              const std::array<Task_t, 3>& numTasksPerDim,
                                                              const std::array<Task_t, 3>& taskPosition) {
   return {
       FsGridTools::calcLocalStart(globalSize[0], numTasksPerDim[0], taskPosition[0]),
       FsGridTools::calcLocalStart(globalSize[1], numTasksPerDim[1], taskPosition[1]),
       FsGridTools::calcLocalStart(globalSize[2], numTasksPerDim[2], taskPosition[2]),
   };
}

constexpr static std::array<FsIndex_t, 3> calculateStorageSize(const std::array<FsSize_t, 3>& globalSize,
                                                               const std::array<Task_t, 3>& localSize,
                                                               int32_t numGhostCells) {
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
   constexpr Coordinates() {}
   constexpr Coordinates(const std::array<double, 3>& physicalGridSpacing,
                         const std::array<double, 3>& physicalGlobalStart, const std::array<FsSize_t, 3>& globalSize,
                         const std::array<bool, 3>& periodic, const std::array<Task_t, 3>& decomposition,
                         const std::array<Task_t, 3>& taskPosition, int32_t numRanks, int32_t numGhostCells)
       : numGhostCells(numGhostCells), physicalGridSpacing(physicalGridSpacing),
         physicalGlobalStart(physicalGlobalStart), globalSize(globalSize), periodic(periodic),
         numTasksPerDim(fsgrid_detail::computeNumTasksPerDim(globalSize, decomposition, numRanks, numGhostCells)),
         localSize(fsgrid_detail::calculateLocalSize(globalSize, numTasksPerDim, taskPosition, numGhostCells)),
         localStart(fsgrid_detail::calculateLocalStart(globalSize, numTasksPerDim, taskPosition)),
         storageSize(fsgrid_detail::calculateStorageSize(globalSize, localSize, numGhostCells)) {}

   /*! Determine the cell's GlobalID from its local x,y,z coordinates
    * \param x The cell's task-local x coordinate
    * \param y The cell's task-local y coordinate
    * \param z The cell's task-local z coordinate
    */
   constexpr GlobalID globalIDFromLocalCoordinates(FsIndex_t x, FsIndex_t y, FsIndex_t z) const {
      // Perform casts to avoid overflow
      const std::array<FsSize_t, 3> global = localToGlobal(x, y, z);
      const auto xcontrib = global[0];
      const auto ycontrib = static_cast<GlobalID>(globalSize[0]) * static_cast<GlobalID>(global[1]);
      const auto zcontrib = static_cast<GlobalID>(globalSize[0]) * static_cast<GlobalID>(globalSize[1]) *
                            static_cast<GlobalID>(global[2]);
      return xcontrib + ycontrib + zcontrib;
   }

   /*! Determine the cell's LocalID from its local x,y,z coordinates
    * \param x The cell's task-local x coordinate
    * \param y The cell's task-local y coordinate
    * \param z The cell's task-local z coordinate
    */
   constexpr LocalID localIDFromLocalCoordinates(FsIndex_t x, FsIndex_t y, FsIndex_t z) const {
      // Perform casts to avoid overflow
      const auto xcontrib = static_cast<LocalID>(globalSize[0] > 1) * static_cast<LocalID>(numGhostCells + x);
      const auto ycontrib = static_cast<LocalID>(globalSize[1] > 1) * static_cast<LocalID>(storageSize[0]) *
                            static_cast<LocalID>(numGhostCells + y);
      const auto zcontrib = static_cast<LocalID>(globalSize[2] > 1) * static_cast<LocalID>(storageSize[0]) *
                            static_cast<LocalID>(storageSize[1]) * static_cast<LocalID>(numGhostCells + z);

      return xcontrib + ycontrib + zcontrib;
   }

   constexpr LocalID localIDFromLocalCoordinates(const std::array<FsIndex_t, 3>& indices) const {
      return localIDFromLocalCoordinates(indices[0], indices[1], indices[2]);
   }

   /*! Transform global cell coordinates into the local domain.
    * If the coordinates are out of bounds, (-1,-1,-1) is returned.
    * \param x The cell's global x coordinate
    * \param y The cell's global y coordinate
    * \param z The cell's global z coordinate
    */
   constexpr std::array<FsIndex_t, 3> globalToLocal(FsSize_t x, FsSize_t y, FsSize_t z) const {
      // Perform this check before doing the subtraction to avoid cases of underflow and overflow
      // Particularly for the first three checks:
      // - casting the localStart to unsigned and then doing the subtraction might cause underflow
      // - casting the global coordinate to signed might overflow, due to global being too large to fit to the signed
      // type
      bool outOfBounds = x < static_cast<FsSize_t>(localStart[0]);
      outOfBounds |= y < static_cast<FsSize_t>(localStart[1]);
      outOfBounds |= z < static_cast<FsSize_t>(localStart[2]);
      outOfBounds |= x >= static_cast<FsSize_t>(localSize[0]) + static_cast<FsSize_t>(localStart[0]);
      outOfBounds |= y >= static_cast<FsSize_t>(localSize[1]) + static_cast<FsSize_t>(localStart[1]);
      outOfBounds |= z >= static_cast<FsSize_t>(localSize[2]) + static_cast<FsSize_t>(localStart[2]);

      if (outOfBounds) {
         return {-1, -1, -1};
      } else {
         // This neither over nor underflows as per the checks above
         return {
             static_cast<FsIndex_t>(x - static_cast<FsSize_t>(localStart[0])),
             static_cast<FsIndex_t>(y - static_cast<FsSize_t>(localStart[1])),
             static_cast<FsIndex_t>(z - static_cast<FsSize_t>(localStart[2])),
         };
      }
   }

   /*! Calculate global cell position (XYZ in global cell space) from local cell coordinates.
    *
    * \param x x-Coordinate, in cells
    * \param y y-Coordinate, in cells
    * \param z z-Coordinate, in cells
    *
    * \return Global cell coordinates
    */
   constexpr std::array<FsSize_t, 3> localToGlobal(FsIndex_t x, FsIndex_t y, FsIndex_t z) const {
      // Cast both before adding to avoid overflow
      return {
          static_cast<FsSize_t>(localStart[0]) + static_cast<FsSize_t>(x),
          static_cast<FsSize_t>(localStart[1]) + static_cast<FsSize_t>(y),
          static_cast<FsSize_t>(localStart[2]) + static_cast<FsSize_t>(z),
      };
   }

   /*! Get the physical coordinates in the global simulation space for
    * the given cell.
    *
    * \param x local x-Coordinate, in cells
    * \param y local y-Coordinate, in cells
    * \param z local z-Coordinate, in cells
    */
   constexpr std::array<double, 3> getPhysicalCoords(FsIndex_t x, FsIndex_t y, FsIndex_t z) const {
      return {
          physicalGlobalStart[0] + (localStart[0] + x) * physicalGridSpacing[0],
          physicalGlobalStart[1] + (localStart[1] + y) * physicalGridSpacing[1],
          physicalGlobalStart[2] + (localStart[2] + z) * physicalGridSpacing[2],
      };
   }

   /*! Get the global cell coordinates for the given physical coordinates.
    *
    * \param x physical x-Coordinate
    * \param y physical y-Coordinate
    * \param z physical z-Coordinate
    */
   constexpr std::array<FsSize_t, 3> physicalToGlobal(double x, double y, double z) const {
      return {
          static_cast<FsSize_t>((x - physicalGlobalStart[0]) / physicalGridSpacing[0]),
          static_cast<FsSize_t>((y - physicalGlobalStart[1]) / physicalGridSpacing[1]),
          static_cast<FsSize_t>((z - physicalGlobalStart[2]) / physicalGridSpacing[2]),
      };
   }

   /*! Get the (fractional) global cell coordinates for the given physical coordinates.
    *
    * \param x physical x-Coordinate
    * \param y physical y-Coordinate
    * \param z physical z-Coordinate
    */
   constexpr std::array<double, 3> physicalToFractionalGlobal(double x, double y, double z) const {
      const auto global = physicalToGlobal(x, y, z);
      return {
          (x - physicalGlobalStart[0]) / physicalGridSpacing[0] - global[0],
          (y - physicalGlobalStart[1]) / physicalGridSpacing[1] - global[1],
          (z - physicalGlobalStart[2]) / physicalGridSpacing[2] - global[2],
      };
   }

   constexpr std::array<FsIndex_t, 3> globalIdToTaskPos(GlobalID id) const {
      const std::array<FsIndex_t, 3> cell = FsGridTools::globalIDtoCellCoord(id, globalSize);

      auto computeIndex = [&](uint32_t i) {
         const FsIndex_t nPerTask = static_cast<FsIndex_t>(globalSize[i] / static_cast<FsSize_t>(numTasksPerDim[i]));
         const FsIndex_t nPerTaskPlus1 = nPerTask + 1;
         const FsIndex_t remainder = static_cast<FsIndex_t>(globalSize[i] % static_cast<FsSize_t>(numTasksPerDim[i]));

         return cell[i] < remainder * nPerTaskPlus1 ? cell[i] / nPerTaskPlus1
                                                    : remainder + (cell[i] - remainder * nPerTaskPlus1) / nPerTask;
      };

      return {
          computeIndex(0),
          computeIndex(1),
          computeIndex(2),
      };
   }

   /*! Compute the neighbourIndex (in range [0, 26]) from the (local) cell indices. Note that the inputs are not local
    * coordinates, but rather coordinates that include the ghost cells as well.
    *
    * \param x cell x-Coordinate
    * \param y cell y-Coordinate
    * \param z cell z-Coordinate
    */
   constexpr FsSize_t neighbourIndexFromCellCoordinates(FsIndex_t x, FsIndex_t y, FsIndex_t z) const {
      return 13 - static_cast<FsSize_t>(x < 0) * 9 + static_cast<FsSize_t>(x >= localSize[0]) * 9 -
             static_cast<FsSize_t>(y < 0) * 3 + static_cast<FsSize_t>(y >= localSize[1]) * 3 -
             static_cast<FsSize_t>(z < 0) + static_cast<FsSize_t>(z >= localSize[2]);
   }

   /*! Compute the (shifted) local indices from the (local) cell indices. Note that the inputs are not local
    * coordinates, but rather coordinates that include the ghost cells as well.
    *
    * \param x cell x-Coordinate
    * \param y cell y-Coordinate
    * \param z cell z-Coordinate
    */
   constexpr std::array<FsIndex_t, 3> shiftCellIndices(FsIndex_t x, FsIndex_t y, FsIndex_t z) const {
      return {
          x + static_cast<FsIndex_t>(x < 0) * localSize[0] - static_cast<FsIndex_t>(x >= localSize[0]) * localSize[0],
          y + static_cast<FsIndex_t>(y < 0) * localSize[1] - static_cast<FsIndex_t>(y >= localSize[1]) * localSize[1],
          z + static_cast<FsIndex_t>(z < 0) * localSize[2] - static_cast<FsIndex_t>(z >= localSize[2]) * localSize[2],
      };
   }

   /*! Check that the (local) cell indices are within our domain. Note that the inputs are not local
    * coordinates, but rather coordinates that include the ghost cells as well.
    *
    * \param x cell x-Coordinate
    * \param y cell y-Coordinate
    * \param z cell z-Coordinate
    */
   constexpr bool cellIndicesAreWithinBounds(FsIndex_t x, FsIndex_t y, FsIndex_t z) const {
      const std::array minimums{
          localSize[0] > 1 || periodic[0] ? -numGhostCells : 0,
          localSize[1] > 1 || periodic[1] ? -numGhostCells : 0,
          localSize[2] > 1 || periodic[2] ? -numGhostCells : 0,
      };
      const std::array maximums{
          localSize[0] > 1 || periodic[0] ? localSize[0] + numGhostCells : 1,
          localSize[1] > 1 || periodic[1] ? localSize[1] + numGhostCells : 1,
          localSize[2] > 1 || periodic[2] ? localSize[2] + numGhostCells : 1,
      };

      return minimums[0] <= x && x < maximums[0] && minimums[1] <= y && y < maximums[1] && minimums[2] <= z &&
             z < maximums[2];
   }
   constexpr auto cellIndicesAreWithinBounds(const std::array<FsIndex_t, 3>& indices) const {
      return cellIndicesAreWithinBounds(indices[0], indices[1], indices[2]);
   }

   // =======================
   // Variables
   // =======================
   const int32_t numGhostCells = 0;
   const std::array<double, 3> physicalGridSpacing = {};
   const std::array<double, 3> physicalGlobalStart = {};
   const std::array<FsSize_t, 3> globalSize = {};
   const std::array<bool, 3> periodic = {};
   const std::array<Task_t, 3> numTasksPerDim = {};
   const std::array<FsIndex_t, 3> localSize = {};
   const std::array<FsIndex_t, 3> localStart = {};
   const std::array<FsIndex_t, 3> storageSize = {};
};
