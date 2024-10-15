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
#include "cartesian_grid.hpp"
#include "tools.hpp"
#include <array>
#include <cstdint>
#include <mpi.h>
#include <vector>

namespace fsgrid {
template <typename T> struct Neighbourhood {
private:
   //!< Lookup table from index to rank in the neighbour array (includes self)
   std::array<int32_t, 27> neighbourIndexToRank = {
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
       MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
   };
   //!< Lookup table from rank to index in the neighbour array
   std::vector<char> neighbourRankToIndex = {};
   //!< Datatype for sending data
   std::array<MPI_Datatype, 27> neighbourSendType = {
       MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL,
       MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL,
       MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL,
       MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL,
       MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL,
   };
   //!< Datatype for receiving data
   std::array<MPI_Datatype, 27> neighbourReceiveType = {
       MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL,
       MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL,
       MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL,
       MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL,
       MPI_DATATYPE_NULL, MPI_DATATYPE_NULL, MPI_DATATYPE_NULL,
   };

public:
   Neighbourhood(const CartesianGrid& cartesianGrid)
       : neighbourIndexToRank(FsGridTools::mapNeigbourIndexToRank(
             cartesianGrid.getTaskPosition(), cartesianGrid.getDecomposition(), cartesianGrid.getPeriodic(),
             cartesianGrid.getComm(), cartesianGrid.getRank())),
         neighbourRankToIndex(FsGridTools::mapNeighbourRankToIndex(neighbourIndexToRank, cartesianGrid.getSize())),
         neighbourSendType(FsGridTools::generateMPITypes<T>(
             cartesianGrid.getStorageSize(), cartesianGrid.getLocalSize(), cartesianGrid.getStencilSize(), true)),
         neighbourReceiveType(FsGridTools::generateMPITypes<T>(
             cartesianGrid.getStorageSize(), cartesianGrid.getLocalSize(), cartesianGrid.getStencilSize(), false)) {}

   void finalize() {
      for (size_t i = 0; i < neighbourReceiveType.size(); i++) {
         if (neighbourReceiveType[i] != MPI_DATATYPE_NULL) {
            FSGRID_MPI_CHECK(MPI_Type_free(&neighbourReceiveType[i]), "Failed to free MPI type");
         }
      }

      for (size_t i = 0; i < neighbourSendType.size(); i++) {
         if (neighbourSendType[i] != MPI_DATATYPE_NULL) {
            FSGRID_MPI_CHECK(MPI_Type_free(&neighbourSendType[i]), "Failed to free MPI type");
         }
      }
   }

   int32_t getNeighbourRank(int32_t neighbourIndex) { return neighbourIndexToRank[neighbourIndex]; }

   void updateGhostCells(int32_t rank, MPI_Comm comm, T* data) {
      std::array<MPI_Request, 27> receiveRequests;
      std::array<MPI_Request, 27> sendRequests;
      receiveRequests.fill(MPI_REQUEST_NULL);
      sendRequests.fill(MPI_REQUEST_NULL);

      for (int32_t x = -1; x <= 1; x++) {
         for (int32_t y = -1; y <= 1; y++) {
            for (int32_t z = -1; z <= 1; z++) {
               const int32_t shiftId = (x + 1) * 9 + (y + 1) * 3 + (z + 1);
               const int32_t receiveId = (1 - x) * 9 + (1 - y) * 3 + (1 - z);
               if (neighbourIndexToRank[receiveId] != MPI_PROC_NULL &&
                   neighbourSendType[shiftId] != MPI_DATATYPE_NULL) {
                  FSGRID_MPI_CHECK(MPI_Irecv(data, 1, neighbourReceiveType[shiftId], neighbourIndexToRank[receiveId],
                                             shiftId, comm, &(receiveRequests[shiftId])),
                                   "Rank ", rank, " failed to receive data from neighbor ", receiveId, " with rank ",
                                   neighbourIndexToRank[receiveId]);
               }
            }
         }
      }

      for (int32_t x = -1; x <= 1; x++) {
         for (int32_t y = -1; y <= 1; y++) {
            for (int32_t z = -1; z <= 1; z++) {
               const int32_t shiftId = (x + 1) * 9 + (y + 1) * 3 + (z + 1);
               const int32_t sendId = shiftId;
               if (neighbourIndexToRank[sendId] != MPI_PROC_NULL && neighbourSendType[shiftId] != MPI_DATATYPE_NULL) {
                  FSGRID_MPI_CHECK(MPI_Isend(data, 1, neighbourSendType[shiftId], neighbourIndexToRank[sendId], shiftId,
                                             comm, &(sendRequests[shiftId])),
                                   "Rank ", rank, " failed to send data to neighbor ", sendId, " with rank ",
                                   neighbourIndexToRank[sendId]);
               }
            }
         }
      }
      FSGRID_MPI_CHECK(MPI_Waitall(27, receiveRequests.data(), MPI_STATUSES_IGNORE),
                       "Synchronization at ghost cell update failed");
      FSGRID_MPI_CHECK(MPI_Waitall(27, sendRequests.data(), MPI_STATUSES_IGNORE),
                       "Synchronization at ghost cell update failed");
   }
};
}
