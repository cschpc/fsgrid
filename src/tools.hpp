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
#include <cassert>
#include <iostream>
#include <limits>
#include <stdint.h>
#include <vector>

namespace FsGridTools {
   typedef uint32_t FsSize_t; // Size type for global array indices
   typedef int32_t FsIndex_t; // Size type for global/local array indices, incl. possible negative values

   typedef int64_t LocalID;
   typedef int64_t GlobalID;
   typedef int Task_t;

   //! Helper function: calculate size of the local coordinate space for the given dimension
   // \param numCells Number of cells in the global Simulation, in this dimension
   // \param nTasks Total number of tasks in this dimension
   // \param taskIndex This task's position in this dimension
   // \return Number of cells for this task's local domain (actual cells, not counting ghost cells)
   static FsIndex_t calcLocalSize(FsSize_t numCells, Task_t nTasks, Task_t taskIndex) {
      const FsIndex_t nPerTask = numCells / nTasks;
      const FsIndex_t remainder = numCells % nTasks;
      return taskIndex < remainder ? nPerTask + 1 : nPerTask;
   }

   //! Helper function: calculate position of the local coordinate space for the given dimension
   // \param numCells number of cells
   // \param nTasks Total number of tasks in this dimension
   // \param taskIndex This task's position in this dimension
   // \return Cell number at which this task's domains cells start (actual cells, not counting ghost cells)
   static FsIndex_t calcLocalStart(FsSize_t numCells, Task_t nTasks, Task_t taskIndex) {
      const FsIndex_t n_per_task = numCells / nTasks;
      const FsIndex_t remainder = numCells % nTasks;
      return taskIndex * calcLocalSize(numCells, nTasks, taskIndex) + (taskIndex >= remainder) * remainder;
   }

   //! Helper function: given a global cellID, calculate the global cell coordinate from it.
   // This is then used do determine the task responsible for this cell, and the
   // local cell index in it.
   static std::array<FsIndex_t, 3> globalIDtoCellCoord(GlobalID id, std::array<FsSize_t, 3>& globalSize) {

      // Transform globalID to global cell coordinate
      std::array<FsIndex_t, 3> cell;

      assert(id >= 0);
      assert(id < globalSize[0] * globalSize[1] * globalSize[2]);

      int stride = 1;
      for (int i = 0; i < 3; i++) {
         cell[i] = (id / stride) % globalSize[i];
         stride *= globalSize[i];
      }

      return cell;
   }

   //! Helper function to optimize decomposition of this grid over the given number of tasks
   static void computeDomainDecomposition(const std::array<FsSize_t, 3>& GlobalSize, Task_t nProcs,
                                          std::array<Task_t, 3>& processDomainDecomposition, int rank,
                                          int stencilSize = 1, int verbose = 0) {
      std::array<FsSize_t, 3> systemDim;
      std::array<FsSize_t, 3> processBox;
      std::array<FsSize_t, 3> minDomainSize;
      int64_t optimValue = std::numeric_limits<int64_t>::max();
      std::vector<std::pair<int64_t, std::array<Task_t, 3>>> scored_decompositions;
      scored_decompositions.push_back(std::pair<int64_t, std::array<Task_t, 3>>(optimValue, {0, 0, 0}));
      for (int i = 0; i < 3; i++) {
         systemDim[i] = GlobalSize[i];
         if (GlobalSize[i] == 1) {
            // In 2D simulation domains, the "thin" dimension can be a single cell thick.
            minDomainSize[i] = 1;
         } else {
            // Otherwise, it needs to be at least as large as our ghost
            // stencil, so that ghost communication remains consistent.
            minDomainSize[i] = stencilSize;
         }
      }
      processDomainDecomposition = {1, 1, 1};
      for (Task_t i = 1; i <= std::min(nProcs, (Task_t)(GlobalSize[0] / minDomainSize[0])); i++) {
         for (Task_t j = 1; j <= std::min(nProcs, (Task_t)(GlobalSize[1] / minDomainSize[1])); j++) {
            if (i * j > nProcs) {
               break;
            }
            Task_t k = nProcs / (i * j);
            // No need to optimize an incompatible DD, also checks for missing remainders
            if (i * j * k != nProcs) {
               continue;
            }
            if (k > (Task_t)(GlobalSize[2] / minDomainSize[2])) {
               continue;
            }

            processBox[0] = calcLocalSize(systemDim[0], i, 0);
            processBox[1] = calcLocalSize(systemDim[1], j, 0);
            processBox[2] = calcLocalSize(systemDim[2], k, 0);

            int64_t value = (i > 1 ? processBox[1] * processBox[2] : 0) + (j > 1 ? processBox[0] * processBox[2] : 0) +
                            (k > 1 ? processBox[0] * processBox[1] : 0);

            // account for singular domains
            if (i != 1 && j != 1 && k != 1) {
               value *= 13; // 26 neighbours to communicate to
            }
            if (i == 1 && j != 1 && k != 1) {
               value *= 4; // 8 neighbours to communicate to
            }
            if (i != 1 && j == 1 && k != 1) {
               value *= 4; // 8 neighbours to communicate to
            }
            if (i != 1 && j != 1 && k == 1) {
               value *= 4; // 8 neighbours to communicate to
            }
            // else: 2 neighbours to communicate to, no need to adjust

            if (value <= optimValue) {
               optimValue = value;
               if (value < scored_decompositions.back().first) {
                  scored_decompositions.clear();
               }
               scored_decompositions.push_back(std::pair<int64_t, std::array<Task_t, 3>>(value, {i, j, k}));
            }
         }
      }

      if (rank == 0 && verbose) {
         std::cout << "(FSGRID) Number of equal minimal-surface decompositions found: " << scored_decompositions.size()
                   << "\n";
         for (auto kv : scored_decompositions) {
            std::cout << "(FSGRID) Decomposition " << kv.second[0] << "," << kv.second[1] << "," << kv.second[2] << " "
                      << " for processBox size " << systemDim[0] / kv.second[0] << " " << systemDim[1] / kv.second[1]
                      << " " << systemDim[2] / kv.second[2] << "\n";
         }
      }

      // Taking the first scored_decomposition (smallest X decomposition)
      processDomainDecomposition[0] = scored_decompositions[0].second[0];
      processDomainDecomposition[1] = scored_decompositions[0].second[1];
      processDomainDecomposition[2] = scored_decompositions[0].second[2];

      if (optimValue == std::numeric_limits<int64_t>::max() ||
          (Task_t)(processDomainDecomposition[0] * processDomainDecomposition[1] * processDomainDecomposition[2]) !=
              nProcs) {
         if (rank == 0) {
            std::cerr << "(FSGRID) Domain decomposition failed, are you running on a prime number of tasks?"
                      << std::endl;
         }
         throw std::runtime_error("FSGrid computeDomainDecomposition failed");
      }
      if (rank == 0 && verbose) {
         std::cout << "(FSGRID) decomposition chosen as " << processDomainDecomposition[0] << " "
                   << processDomainDecomposition[1] << " " << processDomainDecomposition[2] << ", for processBox sizes "
                   << systemDim[0] / processDomainDecomposition[0] << " "
                   << systemDim[1] / processDomainDecomposition[1] << " "
                   << systemDim[2] / processDomainDecomposition[2] << " \n";
      }
   }
}
