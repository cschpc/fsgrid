#pragma once
/*
  Copyright (C) 2016 Finnish Meteorological Institute
  Copyright (C) 2016 CSC -IT Center for Science 

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
#include <array>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <limits>
#include <stdint.h>
#include <cassert>

/*! Simple cartesian, non-loadbalancing MPI Grid for use with the fieldsolver
 *
 * \param T datastructure containing the field in each cell which this grid manages
 * \param stencil ghost cell width of this grid
 */
template <typename T, int stencil> class FsGrid {

   public:

      typedef int64_t LocalID;
      typedef int64_t GlobalID;

      /*! Constructor for this grid.
       * \param globalSize Cell size of the global simulation domain.
       * \param MPI_Comm The MPI communicator this grid should use.
       * \param isPeriodic An array specifying, for each dimension, whether it is to be treated as periodic.
       */
      FsGrid(std::array<int32_t,3> globalSize, MPI_Comm parent_comm, std::array<bool,3> isPeriodic)
            : globalSize(globalSize) {
         int status;
         int size;

         status = MPI_Comm_size(parent_comm, &size);

         // Heuristically choose a good domain decomposition for our field size
         computeDomainDecomposition(globalSize, size, ntasks);
         
         //set private array
         periodic = isPeriodic;
         //set temporary int array for MPI_Cart_create
         std::array<int, 3> isPeriodicInt;
         for(int i=0; i < isPeriodic.size(); i++) {
            isPeriodicInt[i] = (int)isPeriodic[i];
         }  
         
         // Create cartesian communicator. Note, that reorder is false so
         // ranks match the ones in parent_comm
         status = MPI_Cart_create(parent_comm, 3, ntasks.data(), isPeriodicInt.data(), 0, &comm3d);
         if(status != MPI_SUCCESS) {
            std::cerr << "Creating cartesian communicatior failed when attempting to create FsGrid!" << std::endl;
            return;
         }

         status = MPI_Comm_rank(comm3d, &rank);
         if(status != MPI_SUCCESS) {
            std::cerr << "Getting rank failed when attempting to create FsGrid!" << std::endl;

            // Without a rank, there's really not much we can do. Just return an uninitialized grid
            // (I suppose we'll crash after this, anyway)
            return;
         }


         // Determine our position in the resulting task-grid
         status = MPI_Cart_coords(comm3d, rank, 3, taskPosition.data());
         if(status != MPI_SUCCESS) {
            std::cerr << "Rank " << rank
               << " unable to determine own position in cartesian communicator when attempting to create FsGrid!"
               << std::endl;
         }

         // Allocate the array of neighbours
         for(int i=0; i<size; i++) {
            neighbour_index.push_back(MPI_PROC_NULL);
         }
         for(int i=0; i<27; i++) {
            neighbour[i]=MPI_PROC_NULL;
         }

         // Get the IDs of the 26 direct neighbours
         for(int x=-1; x<=1;x++) {
            for(int y=-1; y<=1;y++) {
               for(int z=-1; z<=1; z++) {
                  std::array<int,3> neighPosition;

                  /*
                   * Figure out the coordinates of the neighbours in all three
                   * directions
                   */
                  neighPosition[0]=taskPosition[0]+x;
                  if(isPeriodic[0]) {
                     neighPosition[0] += ntasks[0];
                     neighPosition[0] %= ntasks[0];
                  }

                  neighPosition[1]=taskPosition[1]+y;
                  if(isPeriodic[1]) {
                     neighPosition[1] += ntasks[1];
                     neighPosition[1] %= ntasks[1];
                  }

                  neighPosition[2]=taskPosition[2]+z;
                  if(isPeriodic[2]) {
                     neighPosition[2] += ntasks[2];
                     neighPosition[2] %= ntasks[2];
                  }

                  /*
                   * If those coordinates exist, figure out the responsible CPU
                   * and store its rank
                   */
                  if(neighPosition[0]>=0 && neighPosition[0]<ntasks[0] && neighPosition[1]>=0
                        && neighPosition[1]<ntasks[1] && neighPosition[2]>=0 && neighPosition[2]<ntasks[2]) {

                     // Calculate the rank
                     int neighRank;
                     status = MPI_Cart_rank(comm3d, neighPosition.data(), &neighRank);
                     if(status != MPI_SUCCESS) {
                        std::cerr << "Rank " << rank << " can't determine neighbour rank at position [";
                        for(int i=0; i<3; i++) {
                           std::cerr << neighPosition[i] << ", ";
                        }
                        std::cerr << "]" << std::endl;
                     }

                     // Forward lookup table
                     neighbour[(x+1)*9+(y+1)*3+(z+1)]=neighRank;

                     // Reverse lookup table
                     if(neighRank >= 0 && neighRank < size) {
                        neighbour_index[neighRank]=(char) ((x+1)*9+(y+1)*3+(z+1));
                     }
                  } else {
                     neighbour[(x+1)*9+(y+1)*3+(z+1)]=MPI_PROC_NULL;
                  }
               }
            }
         }


         // Determine size of our local grid
         for(int i=0; i<3; i++) {
            localSize[i] = calcLocalSize(globalSize[i],ntasks[i], taskPosition[i]);
            localStart[i] = calcLocalStart(globalSize[i],ntasks[i], taskPosition[i]);
         }

         //std::cerr << "Rank " << rank << " here, my space starts at [" << localStart[0] << ", " << localStart[1] << ", "
         //   << localStart[2] << "] and ends at [" << (localStart[0] + localSize[0]) << ", " <<
         //   (localStart[1]+localSize[1]) << ", " << (localStart[2]+localSize[2]) << "]" << std::endl;

         // Allocate local storage array
         size_t totalStorageSize=1;
         for(int i=0; i<3; i++) {
            if(globalSize[i] <= 1) {
               // Collapsed dimension => only one cell thick
               storageSize[i] = 1;
            } else {
               // Size of the local domain + 2* size for the ghost cell stencil
               storageSize[i] = (localSize[i] + stencil*2);
            }
            totalStorageSize *= storageSize[i];
         }
         data.resize(totalStorageSize);


         MPI_Datatype mpiTypeT;
         MPI_Type_contiguous(sizeof(T), MPI_BYTE, &mpiTypeT);
         for(int x=-1; x<=1;x++) {
            for(int y=-1; y<=1;y++) {
               for(int z=-1; z<=1; z++) {
                  neighbourSendType[(x+1) * 9 + (y + 1) * 3 + (z + 1)] = MPI_DATATYPE_NULL;
                  neighbourReceiveType[(x+1) * 9 + (y + 1) * 3 + (z + 1)] = MPI_DATATYPE_NULL;
               }
            }
         }
         
         // Compute send and receive datatypes

         //loop through the shifts in the different directions
         for(int x=-1; x<=1;x++) {
            for(int y=-1; y<=1;y++) {
               for(int z=-1; z<=1; z++) {
                  std::array<int,3> subarraySize;
                  std::array<int,3> subarrayStart;                  
                  const int shiftId = (x+1) * 9 + (y + 1) * 3 + (z + 1);
                  
                  
                  if((storageSize[0] == 1 && x!= 0 ) ||
                     (storageSize[1] == 1 && y!= 0 ) ||
                     (storageSize[2] == 1 && z!= 0 ) ||
                     (x == 0 && y == 0 && z == 0)){
                     //skip flat dimension for 2 or 1D simulations, and self
                     neighbourSendType[shiftId] = MPI_DATATYPE_NULL;
                     neighbourReceiveType[shiftId] = MPI_DATATYPE_NULL;
                     continue;
                  }

                  subarraySize[0] = (x == 0) ? localSize[0] : stencil;
                  subarraySize[1] = (y == 0) ? localSize[1] : stencil;
                  subarraySize[2] = (z == 0) ? localSize[2] : stencil;

                  if( x == 0 || x == -1 )
                     subarrayStart[0] = stencil;
                  else if (x == 1)
                     subarrayStart[0] = storageSize[0] - 2 * stencil;
                  if( y == 0 || y == -1 )
                     subarrayStart[1] = stencil;
                  else if (y == 1)
                     subarrayStart[1] = storageSize[1] - 2 * stencil;
                  if( z == 0 || z == -1 )
                     subarrayStart[2] = stencil;
                  else if (z == 1)
                     subarrayStart[2] = storageSize[2] - 2 * stencil;
                  
                  for(int i = 0;i < 3; i++)
                     if(storageSize[i] == 1) 
                        subarrayStart[i] = 0;

                  std::array<int,3> swappedStorageSize = storageSize;
                  swapArray(swappedStorageSize);
                  swapArray(subarraySize);
                  swapArray(subarrayStart);                  
                  MPI_Type_create_subarray(3,
                                           swappedStorageSize.data(),
                                           subarraySize.data(),
                                           subarrayStart.data(),
                                           MPI_ORDER_C,
                                           mpiTypeT,
                                           &(neighbourSendType[shiftId]) );
                  
                  if(x == 1 )
                     subarrayStart[0] = 0;
                  else if (x == 0)
                     subarrayStart[0] = stencil;
                  else if (x == -1)
                     subarrayStart[0] = storageSize[0] -  stencil;
                  if(y == 1 )
                     subarrayStart[1] = 0;
                  else if (y == 0)
                     subarrayStart[1] = stencil;
                  else if (y == -1)
                     subarrayStart[1] = storageSize[1] -  stencil;
                  if(z == 1 )
                     subarrayStart[2] = 0;
                  else if (z == 0)
                     subarrayStart[2] = stencil;
                  else if (z == -1)
                     subarrayStart[2] = storageSize[2] -  stencil;
                  for(int i = 0;i < 3; i++)
                     if(storageSize[i] == 1) 
                        subarrayStart[i] = 0;
                  
                  swapArray(subarrayStart);                  
                  MPI_Type_create_subarray(3,
                                           swappedStorageSize.data(),
                                           subarraySize.data(),
                                           subarrayStart.data(),
                                           MPI_ORDER_C,
                                           mpiTypeT,
                                           &(neighbourReceiveType[shiftId]));
                  
               }
            }
         }
         
         for(int i=0;i<27;i++){
            if(neighbourReceiveType[i] != MPI_DATATYPE_NULL)
               MPI_Type_commit(&(neighbourReceiveType[i]));
            if(neighbourSendType[i] != MPI_DATATYPE_NULL)
               MPI_Type_commit(&(neighbourSendType[i]));
         }
         
         // Also set up coupling information to external grid (fill with MPI_PROC_NULL to begin with)
         // Before actual grid coupling can be done, this information has to be filled in.
         externalRank.resize(totalStorageSize);
         for(uint i=0; i<externalRank.size(); i++) {
            externalRank[i] = MPI_PROC_NULL;
         }


      }

      /*! Destructor, cleans up the cartesian communicator
       */
      ~FsGrid() {
         for(int i=0;i<27;i++){
            if(neighbourReceiveType[i] != MPI_DATATYPE_NULL)
               MPI_Type_free(&(neighbourReceiveType[i]));
            if(neighbourSendType[i] != MPI_DATATYPE_NULL)
               MPI_Type_free(&(neighbourSendType[i]));
         }
         MPI_Comm_free(&comm3d);
      }

      /*! Returns the task responsible, and its localID for handling the cell with the given GlobalID
       * \param id GlobalID of the cell for which task is to be determined
       * \return a task for the grid's cartesian communicator
       */
      std::pair<int,LocalID> getTaskForGlobalID(GlobalID id) {
         // Transform globalID to global cell coordinate
         std::array<int, 3> cell = globalIDtoCellCoord(id);

         // Find the index in the task grid this Cell belongs to
         std::array<int, 3> taskIndex;
         for(int i=0; i<3; i++) {
            int n_per_task = globalSize[i] / ntasks[i];
            int remainder = globalSize[i] % ntasks[i];

            if(cell[i] < remainder * (n_per_task+1)) {
               taskIndex[i] = cell[i] / (n_per_task + 1);
            } else {
               taskIndex[i] = remainder + (cell[i] - remainder*(n_per_task+1)) / n_per_task;
            }
         }

         // Get the task number from the communicator
         std::pair<int,LocalID> retVal;
         int status = MPI_Cart_rank(comm3d, taskIndex.data(), &retVal.first);
         if(status != MPI_SUCCESS) {
            std::cerr << "Unable to find FsGrid rank for global ID " << id << " (coordinates [";
            for(int i=0; i<3; i++) {
               std::cerr << cell[i] << ", ";
            }
            std::cerr << "]" << std::endl;
            return std::pair<int,LocalID>(MPI_PROC_NULL,0);
         }

         // Determine localID of that cell within the target task
         std::array<int, 3> thatTasksStart;
         std::array<int, 3> thatTaskStorageSize;
         for(int i=0; i<3; i++) {
            thatTasksStart[i] = calcLocalStart(globalSize[i], ntasks[i], taskIndex[i]);
            thatTaskStorageSize[i] = calcLocalSize(globalSize[i], ntasks[i], taskIndex[i]) + 2 * stencil;
         }

         retVal.second = 0;
         int stride = 1;
         for(int i=0; i<3; i++) {
            if(globalSize[i] <= 1) {
               // Collapsed dimension, doesn't contribute.
               retVal.second += 0;
            } else {
               retVal.second += stride*(cell[i] - thatTasksStart[i] + stencil);
               stride *= thatTaskStorageSize[i];
            }
         }
         
         return retVal;
      }

      /*! Determine the cell's GlobalID from its local x,y,z coordinates
       * \param x The cell's task-local x coordinate
       * \param y The cell's task-local y coordinate
       * \param z The cell's task-local z coordinate
       */
      GlobalID GlobalIDForCoords(int x, int y, int z) {
         return globalSize[0]*(x+localStart[0])
            + globalSize[1]*(y+localStart[1])
            + globalSize[2]*(z+localStart[2]);
      }
      /*! Determine the cell's LocalID from its local x,y,z coordinates
       * \param x The cell's task-local x coordinate
       * \param y The cell's task-local y coordinate
       * \param z The cell's task-local z coordinate
       */
      LocalID LocalIDForCoords(int x, int y, int z) {
         LocalID index=0;
         if(globalSize[2] > 1) {
            index += storageSize[0]*storageSize[1]*(stencil+z);
         }
         if(globalSize[1] > 1) {
            index += storageSize[0]*(stencil+y);
         }
         if(globalSize[0] > 1 ) {
            index += stencil+x;
         }

         return index;
      }

      /*! Prepare for transfer of Rank information from a coupled grid.
       * Setup MPI Irecv requests to recieve data for each grid cell, and prepare
       * buffer space to hold the MPI requests of the sends.
       *
       * \param cellsToCouple How many cells are going to be coupled
       */
      void setupForGridCoupling(int cellsToCouple) {
         int status;
         // Make sure we have sufficient buffer space to store our mpi requests
         requests.resize(localSize[0]*localSize[1]*localSize[2] + cellsToCouple);
         numRequests=0;

         for(int z=0; z<localSize[2]; z++) {
            for(int y=0; y<localSize[1]; y++) {
               for(int x=0; x<localSize[0]; x++) {
                  // Calculate LocalID for this cell
                  LocalID thisCell = LocalIDForCoords(x,y,z);
                  assert(numRequests < requests.size());
                  status = MPI_Irecv(&externalRank[thisCell], 1, MPI_INT, MPI_ANY_SOURCE, thisCell, comm3d,
                        &requests[numRequests++]);
                  if(status != MPI_SUCCESS) {
                     std::cerr << "Error setting up MPI Irecv in FsGrid::setupForGridCoupling" << std::endl;
                  }
               }
            }
         }
      }

      /*! Set the MPI rank responsible for external communication with the given cell.
       * If, for example, a separate grid is used in another part of the code, this would be the rank
       * in that grid, responsible for the information in this cell.
       *
       * This only needs to be called if the association changes (such as: on startup, and in a load
       * balancing step)
       *
       * \param id Global cell ID to be addressed
       * \iparam rank
       */
      void setGridCoupling(GlobalID id, int rank) {

         // Determine Task and localID that this cell belongs to
         std::pair<int,LocalID> TaskLid = getTaskForGlobalID(id);

         // Build the MPI Isend request for this cell
         int status;
         assert(numRequests < requests.size());
         status = MPI_Isend(&rank, 1, MPI_INT, TaskLid.first, TaskLid.second, comm3d,
               &requests[numRequests++]);
         if(status != MPI_SUCCESS) {
            std::cerr << "Error setting up MPI Isend in FsGrid::setGridCoupling" << std::endl;
         }
      }

      /*! Called after setting up the transfers into or out of this grid.
       * Basically only does a MPI_Waitall for all requests.
       */
      void finishGridCoupling() {
         assert(numRequests == requests.size());
         MPI_Waitall(numRequests, requests.data(), MPI_STATUSES_IGNORE);
         numRequests=0;
      }

      /*! Prepare for transfer of grid cell data into this grid. 
       * Setup MPI Irecv requests to recieve data for each grid cell, and prepare
       * buffer space to hold the MPI requests of the sends.
       *
       * \param cellsToSend How many cells are going to be sent by this task
       */
      void setupForTransferIn(int cellsToSend) {
         int status;
         // Make sure we have sufficient buffer space to store our mpi requests
         requests.resize(localSize[0]*localSize[1]*localSize[2] + cellsToSend);
         numRequests=0;

         for(int z=0; z<localSize[2]; z++) {
            for(int y=0; y<localSize[1]; y++) {
               for(int x=0; x<localSize[0]; x++) {
                  // Calculate LocalID for this cell
                  LocalID thisCell = LocalIDForCoords(x,y,z);
                  assert(numRequests < requests.size());
                  status = MPI_Irecv(get(thisCell), sizeof(T), MPI_BYTE, externalRank[thisCell], thisCell, comm3d,
                        &requests[numRequests++]);
                  if(status != MPI_SUCCESS) {
                     std::cerr << "Error setting up MPI Irecv in FsGrid::setupForTransferIn" << std::endl;
                  }
               }
            }
         }
      }

      /*! Set grid cell wtih the given ID to the value specified. Note that this doesn't need to be a local cell.
       * Note that this function should be concurrently called by all tasks, to allow for point-to-point communication.
       * \param id Global cell ID to be filled
       * \iparam value New value to be assigned to it
       */
      void transferDataIn(GlobalID id, T* value) {

         // Determine Task and localID that this cell belongs to
         std::pair<int,LocalID> TaskLid = getTaskForGlobalID(id);

         // Build the MPI Isend request for this cell
         int status;
         assert(numRequests < requests.size());
         status = MPI_Isend(value, sizeof(T), MPI_BYTE, TaskLid.first, TaskLid.second, comm3d,
               &requests[numRequests++]);
         if(status != MPI_SUCCESS) {
            std::cerr << "Error setting up MPI Isend in FsGrid::transferDataIn" << std::endl;
         }
      }
      /*! Called after setting up the transfers into or out of this grid.
       * Basically only does a MPI_Waitall for all requests.
       */
      void finishTransfersIn() {
         assert(numRequests == requests.size());
         MPI_Waitall(numRequests,requests.data(),MPI_STATUSES_IGNORE);
         numRequests = 0;
         
      }


      /*! Prepare for transfer of grid cell data out of this grid.
       * Setup MPI Isend requests to send data for each grid cell, and prepare
       * buffer space to hold the MPI requests of the sends.
       *
       * \param cellsToSend How many cells are going to be sent by this task
       */
      void setupForTransferOut(int cellsToRecieve) {
         // Make sure we have sufficient buffer space to store our mpi requests
         requests.resize(localSize[0]*localSize[1]*localSize[2] + cellsToRecieve);
         numRequests=0;
      }

      /*! Get the value of the grid cell wtih the given ID. Note that this doesn't need to be a local cell.
       * \param id Global cell ID to be read
       * \iparam target Location that the result is to be stored in
       */
      void transferDataOut(GlobalID id, T* target) {

         // Determine Task and localID that this cell belongs to
         std::pair<int,LocalID> TaskLid = getTaskForGlobalID(id);

         // Build the MPI Irecv request for this cell
         int status;
         assert(numRequests < requests.size());
         status = MPI_Irecv(&target, sizeof(T), MPI_BYTE, TaskLid.first, TaskLid.second, comm3d,
               &requests[numRequests++]);
         if(status != MPI_SUCCESS) {
            std::cerr << "Error setting up MPI Isend in FsGrid::transferDataOut" << std::endl;
         }
      }

      /*! Finalize the transfer of data out of this grid, re-integrating it back
       *! into the locations where it came from.
       * This method assumes that transferDataOut() has been called for each cell beforehand,
       * so that appropriate recieve buffers exist in the target ranks.
       */
      void finishTransfersOut() {
         int status;
         for(int z=0; z<localSize[2]; z++) {
            for(int y=0; y<localSize[1]; y++) {
               for(int x=0; x<localSize[0]; x++) {
                  // Calculate LocalID for this cell
                  LocalID thisCell = LocalIDForCoords(x,y,z);
                  assert(numRequests < requests.size());
                  status = MPI_Isend(get(thisCell), sizeof(T), MPI_BYTE, externalRank[thisCell], thisCell, comm3d,
                        &requests[numRequests++]);
                  if(status != MPI_SUCCESS) {
                     std::cerr << "Error setting up MPI Irecv in FsGrid::setupForTransferIn" << std::endl;
                  }
               }
            }
         }
         assert(numRequests == requests.size());
         MPI_Waitall(numRequests,requests.data(),MPI_STATUSES_IGNORE);
         numRequests = 0;
      }

      /*! Perform ghost cell communication.
       */
      void updateGhostCells() {
         //TODO, faster with simultaneous isends& ireceives?
         std::array<MPI_Request, 27> receiveRequests;
         std::array<MPI_Request, 27> sendRequests;
         
         for(int i = 0; i < 27; i++){
            receiveRequests[i] = MPI_REQUEST_NULL;
            sendRequests[i] = MPI_REQUEST_NULL;
         }
         
         
         for(int x=-1; x<=1;x++) {
            for(int y=-1; y<=1;y++) {
               for(int z=-1; z<=1; z++) {
                  int shiftId = (x+1) * 9 + (y + 1) * 3 + (z + 1);
                  int receiveId = (1 - x) * 9 + ( 1 - y) * 3 + ( 1 - z);
                  if(neighbour[receiveId] != MPI_PROC_NULL &&
                     neighbourSendType[shiftId] != MPI_DATATYPE_NULL) {
                     MPI_Irecv(data.data(), 1, neighbourReceiveType[shiftId], neighbour[receiveId], shiftId, comm3d, &(receiveRequests[shiftId]));
                  }
               }
            }
         }
         
         for(int x=-1; x<=1;x++) {
            for(int y=-1; y<=1;y++) {
               for(int z=-1; z<=1; z++) {
                  int shiftId = (x+1) * 9 + (y + 1) * 3 + (z + 1);
                  int sendId = shiftId;
                  if(neighbour[sendId] != MPI_PROC_NULL &&
                     neighbourSendType[shiftId] != MPI_DATATYPE_NULL) {
                     MPI_Isend(data.data(), 1, neighbourSendType[shiftId], neighbour[sendId], shiftId, comm3d, &(sendRequests[shiftId]));
                  }
               }
            }
         }
         MPI_Waitall(27, receiveRequests.data(), MPI_STATUSES_IGNORE);         
         MPI_Waitall(27, sendRequests.data(), MPI_STATUSES_IGNORE);
      }
   
   
   

      /*! Get the size of the local domain handled by this grid.
       */
      std::array<int32_t, 3>& getLocalSize() {
         return localSize;
      }

      /*! Calculate global cell position (XYZ in global cell space) from local cell coordinates.
       *
       * \param x x-Coordinate, in cells
       * \param y y-Coordinate, in cells
       * \param z z-Coordinate, in cells
       *
       * \return Global cell coordinates
       */
      std::array<int32_t, 3> getGlobalIndices(int x, int y, int z) {
         std::array<int32_t, 3> retval;
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
      T* get(int x, int y, int z) {

         // Santiy-Check that the requested cell is actually inside our domain
         // TODO: ugh, this is ugly.
         bool inside=true;
         if(localSize[0] <= 1 && !periodic[0]) {
            if(x != 0) {
               std::cerr << "x != 0 despite non-periodic x-axis with only one cell." << std::endl;
               inside = false;
            }
         } else {
            if(x < -stencil || x > localSize[0] + stencil) {
               std::cerr << "x = " << x << " is outside of [ " << -stencil << 
                  ", " << localSize[0] + stencil << "]!" << std::endl;
               inside = false;
            }
         }

         if(localSize[1] <= 1 && !periodic[1]) {
            if(y != 0) {
               std::cerr << "y != 0 despite non-periodic y-axis with only one cell." << std::endl;
               inside = false;
            }
         } else {
            if(y < -stencil || y > localSize[1] + stencil) {
               std::cerr << "y = " << y << " is outside of [ " << -stencil <<
                  ", " << localSize[1] + stencil << "]!" << std::endl;
               inside = false;
            }
         }

         if(localSize[2] <= 1 && !periodic[2]) {
            if(z != 0) {
               std::cerr << "z != 0 despite non-periodic z-axis with only one cell." << std::endl;
               inside = false;
            }
         } else {
            if(z < -stencil || z > localSize[2] + stencil) {
               inside = false;
               std::cerr << "z = " << z << " is outside of [ " << -stencil <<
                  ", " << localSize[2] + stencil << "]!" << std::endl;
            }
         }

         if(!inside) {
            std::cerr << "Out-of bounds access in FsGrid::get! Expect weirdness." << std::endl;
            return NULL;
         }
         
         LocalID index = LocalIDForCoords(x,y,z);

         return &data[index];
      }

      T* get(LocalID id) {
         if(id < 0 || id > data.size()) {
            std::cerr << "Out-of-bounds access in FsGrid::get!" << std::endl
               << "(LocalID = " << id << ", but storage space is " << data.size()
               << ". Expect weirdness." << std::endl;
            return NULL;
         }
         return &data[id];
      }

      /*! Physical grid spacing.
       * TODO: Should this be private and have accesor-functions?
       */
      double DX,DY,DZ;

   private:
      //! MPI Cartesian communicator used in this grid
      MPI_Comm comm3d;
      int rank; //!< This task's rank in the communicator
      std::vector<MPI_Request> requests;
      uint numRequests;

      std::vector<int> externalRank; //!< MPI rank that each cell is being communicated to externally

      std::array<int, 27> neighbour; //!< Tasks of the 26 neighbours (plus ourselves)
      std::vector<char> neighbour_index; //!< Lookup table from rank to index in the neighbour array

      // We have, fundamentally, two different coordinate systems we're dealing with:
      // 1) Task grid in the MPI_Cartcomm
      std::array<int, 3> ntasks; //!< Number of tasks in each direction
      std::array<int, 3> taskPosition; //!< This task's position in the 3d task grid
      // 2) Cell numbers in global and local view

      std::array<bool, 3> periodic; //!< Information about whether a given direction is periodic
      std::array<int32_t, 3> globalSize; //!< Global size of the simulation space, in cells
      std::array<int32_t, 3> localSize;  //!< Local size of simulation space handled by this task (without ghost cells)
      std::array<int32_t, 3> storageSize;  //!< Local size of simulation space handled by this task (including ghost cells)
      std::array<int32_t, 3> localStart; //!< Offset of the local
                                          //!coordinate system against
                                          //!the global one

      std::array<MPI_Datatype, 27> neighbourSendType; //!< Datatype for sending data
      std::array<MPI_Datatype, 27> neighbourReceiveType; //!< Datatype for receiving data



      //! Actual storage of field data
      std::vector<T> data;

      //! Helper function: given a global cellID, calculate the global cell coordinate from it.
      // This is then used do determine the task responsible for this cell, and the
      // local cell index in it.
      std::array<int, 3> globalIDtoCellCoord(GlobalID id) {

         // Transform globalID to global cell coordinate
         std::array<int, 3> cell;

         assert(id >= 0);
         assert(id < globalSize[0] * globalSize[1] * globalSize[2]);

         int stride=1;
         for(int i=0; i<3; i++) {
            cell[i] = (id / stride) % globalSize[i];
            stride *= globalSize[i];
         }

         return cell;
      }

      //! Helper function to optimize decomposition of this grid over the given number of tasks

   
      void computeDomainDecomposition(const std::array<int, 3>& GlobalSize, int nProcs, std::array<int,3>& processDomainDecomposition) {
         std::array<double, 3> systemDim;
         std::array<double, 3 > processBox;
         double optimValue = std::numeric_limits<double>::max();
         for(int i = 0; i < 3; i++) {
            systemDim[i] = (double)GlobalSize[i];
         }
         processDomainDecomposition = {1, 1, 1};
         for (int i = 1; i <= std::min(nProcs, GlobalSize[0]); i++) {
            processBox[0] = std::max(systemDim[0]/i, 1.0);
            for (int j = 1; j <= std::min(nProcs, GlobalSize[1]) ; j++) {
               if( i * j  > nProcs )
                  break;
               processBox[1] = std::max(systemDim[1]/j, 1.0);
               for (int k = 1; k <= std::min(nProcs, GlobalSize[2]); k++) {
                  if( i * j * k > nProcs )
                     break;
                  processBox[2] = std::max(systemDim[2]/k, 1.0);
                  double value = 
                     10 * processBox[0] * processBox[1] * processBox[2] + 
                     (i > 1 ? processBox[1] * processBox[2]: 0) +
                     (j > 1 ? processBox[0] * processBox[2]: 0) +
                     (k > 1 ? processBox[0] * processBox[1]: 0);
	
                  if(value < optimValue ){
                     optimValue = value;
                     processDomainDecomposition[0] = i;
                     processDomainDecomposition[1] = j;
                     processDomainDecomposition[2] = k;
                  }
               }
            }
         }
      }
   
      //! Helper function: calculate position of the local coordinate space for the given dimension
      // \param globalCells Number of cells in the global Simulation, in this dimension
      // \param ntasks Total number of tasks in this dimension
      // \param my_n This task's position in this dimension
      // \return Cell number at which this task's domains cells start (actual cells, not counting ghost cells)
      int32_t calcLocalStart(int32_t globalCells, int ntasks, int my_n) {
         int n_per_task = globalCells / ntasks;
         int remainder = globalCells % ntasks;

         if(my_n < remainder) {
            return my_n * (n_per_task+1);
         } else {
            return my_n * n_per_task + remainder;
         }
      }
      //! Helper function: calculate size of the local coordinate space for the given dimension
      // \param globalCells Number of cells in the global Simulation, in this dimension
      // \param ntasks Total number of tasks in this dimension
      // \param my_n This task's position in this dimension
      // \return Nmuber of cells for this task's local domain (actual cells, not counting ghost cells)
      int32_t calcLocalSize(int32_t globalCells, int ntasks, int my_n) {
         int n_per_task = globalCells/ntasks;
         int remainder = globalCells%ntasks;
         if(my_n < remainder) {
            return n_per_task+1;
         } else {
            return n_per_task;
         }
      }
   
      void swapArray(std::array<int, 3>& array) {
         int a = array[0];
         array[0] = array[2];
         array[2] = a;
      }
};
