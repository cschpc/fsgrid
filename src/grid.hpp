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
#include "neighbourhood.hpp"
#include "tools.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <mpi.h>
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

public:
   /*! Constructor for this grid.
    * \param globalSize Cell size of the global simulation domain.
    * \param MPI_Comm The MPI communicator this grid should use.
    * \param periodic An array specifying, for each dimension, whether it is to be treated as periodic.
    */
   FsGrid(std::array<FsSize_t, 3> globalSize, MPI_Comm parentComm, std::array<bool, 3> periodic,
          const std::array<Task_t, 3>& decomposition = {0, 0, 0})
       : cartesianGrid(fsgrid::CartesianGrid(globalSize, parentComm, periodic, decomposition, stencil)),
         neighbourhood(fsgrid::Neighbourhood<T>(cartesianGrid)),
         data(cartesianGrid.getRank() == -1 ? 0 : cartesianGrid.getNumTotalStored()) {}

   //!< Finalize instead of destructor, as the MPI calls fail after the main program called MPI_Finalize().
   void finalize() {
      cartesianGrid.finalize();
      neighbourhood.finalize();
   }

   // ============================
   // Data access functions
   // ============================
   std::vector<T>& getData() { return data; }

   void copyData(FsGrid& other) {
      // Copy assignment
      data = other.getData();
   }

   // ------------------------------------------------------------
   // Neighbourhood stuff
   // ------------------------------------------------------------
   /*! Perform ghost cell communication.
    */
   void updateGhostCells() {
      const auto rank = cartesianGrid.getRank();
      const auto comm = cartesianGrid.getComm();

      if (rank == -1) {
         return;
      }

      neighbourhood.updateGhostCells(rank, comm, data.data());
   }

   /*! Get a reference to the field data in a cell
    * \param x x-Coordinate, in cells
    * \param y y-Coordinate, in cells
    * \param z z-Coordinate, in cells
    * \return A reference to cell data in the given cell
    */
   T* get(int32_t x, int32_t y, int32_t z) {

      // Keep track which neighbour this cell actually belongs to (13 = ourself)
      // TODO: compute in cartesian grid
      int32_t coord_shift[3] = {0, 0, 0};
      if (x < 0) {
         coord_shift[0] = 1;
      }
      if (x >= cartesianGrid.getLocalSize()[0]) {
         coord_shift[0] = -1;
      }
      if (y < 0) {
         coord_shift[1] = 1;
      }
      if (y >= cartesianGrid.getLocalSize()[1]) {
         coord_shift[1] = -1;
      }
      if (z < 0) {
         coord_shift[2] = 1;
      }
      if (z >= cartesianGrid.getLocalSize()[2]) {
         coord_shift[2] = -1;
      }

      const int32_t neighbourIndex = cartesianGrid.computeNeighbourIndex(x, y, z);
      const auto neigbourRank = neighbourhood.getNeighbourRank(neighbourIndex);
      if (neighbourIndex != 13) {
         // Check if the corresponding neighbour exists
         if (neigbourRank == MPI_PROC_NULL) {
            // TODO: change to assert?
            // Neighbour doesn't exist, we must be an outer boundary cell
            // (or something is quite wrong)
            return NULL;
         } else if (neigbourRank == cartesianGrid.getRank()) {
            // For periodic boundaries, where the neighbour is actually ourself,
            // return our own actual cell instead of the ghost
            x += coord_shift[0] * cartesianGrid.getLocalSize()[0];
            y += coord_shift[1] * cartesianGrid.getLocalSize()[1];
            z += coord_shift[2] * cartesianGrid.getLocalSize()[2];
         }
         // Otherwise we return the ghost cell
      }
      LocalID index = cartesianGrid.getLocalIDForCoords(x, y, z);

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
      return {physicalGlobalStart[0] + (cartesianGrid.getLocalStart()[0] + x) * DX,
              physicalGlobalStart[1] + (cartesianGrid.getLocalStart()[1] + y) * DY,
              physicalGlobalStart[2] + (cartesianGrid.getLocalStart()[2] + z) * DZ};
   }

   /*! Perform an MPI_Allreduce with this grid's internal communicator
    * Function syntax is identical to MPI_Allreduce, except the final (communicator
    * argument will not be needed) */
   int32_t Allreduce(void* sendbuf, void* recvbuf, int32_t count, MPI_Datatype datatype, MPI_Op op) {
      // If a normal FS-rank
      if (cartesianGrid.getRank() != -1) {
         return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, cartesianGrid.getComm());
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

   // ------------------------------------------------------------
   // Debug stuff
   // ------------------------------------------------------------
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
      if (cartesianGrid.getLocalSize()[0] > 1) {
         xmin = -stencil;
         xmax = cartesianGrid.getLocalSize()[0] + stencil;
      }
      if (cartesianGrid.getLocalSize()[1] > 1) {
         ymin = -stencil;
         ymax = cartesianGrid.getLocalSize()[1] + stencil;
      }
      if (cartesianGrid.getLocalSize()[2] > 1) {
         zmin = -stencil;
         zmax = cartesianGrid.getLocalSize()[2] + stencil;
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

   // TODO: move these somewhere?
   /*! Physical grid spacing and physical coordinate space start.
    */
   double DX = 0.0;
   double DY = 0.0;
   double DZ = 0.0;
   std::array<double, 3> physicalGlobalStart = {};

private:
   fsgrid::CartesianGrid cartesianGrid = {};
   fsgrid::Neighbourhood<T> neighbourhood = {};

   //! Actual storage of field data
   std::vector<T> data = {};
};
