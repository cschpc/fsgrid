#include "grid.hpp"
#include <gtest/gtest.h>
#include <mpi.h>
#include <tools.hpp>

TEST(FsGridTest, localToGlobalRoundtrip1) {
   const std::array<fsgrid::FsSize_t, 3> globalSize{1024, 666, 71};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, true, false};
   auto numProcs = 0;
   MPI_Comm_size(parentComm, &numProcs);

   const auto grid = fsgrid::FsGrid<std::array<double, 15>, 1>(globalSize, parentComm, numProcs, periodic,
                                                               {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   for (int32_t x = 0; x < localSize[0]; x++) {
      for (int32_t y = 0; y < localSize[1]; y++) {
         for (int32_t z = 0; z < localSize[2]; z++) {
            const auto global = grid.localToGlobal(x, y, z);
            const auto local = grid.globalToLocal(global[0], global[1], global[2]);
            ASSERT_EQ(local[0], x);
            ASSERT_EQ(local[1], y);
            ASSERT_EQ(local[2], z);
         }
      }
   }
}

TEST(FsGridTest, myGlobalIDCorrespondsToMyTask) {
   int rank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   const std::array<fsgrid::FsSize_t, 3> globalSize{6547, 16, 77};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, false, false};
   auto numProcs = 0;
   MPI_Comm_size(parentComm, &numProcs);

   const auto grid = fsgrid::FsGrid<std::array<double, 6>, 1>(globalSize, parentComm, numProcs, periodic,
                                                              {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   for (int32_t x = 0; x < localSize[0]; x++) {
      for (int32_t y = 0; y < localSize[1]; y++) {
         for (int32_t z = 0; z < localSize[2]; z++) {
            const auto gid = grid.globalIDFromLocalCoordinates(x, y, z);
            const auto task = grid.getTaskForGlobalID(gid);
            ASSERT_EQ(task, rank);
            ASSERT_EQ(task, rank);
            ASSERT_EQ(task, rank);
         }
      }
   }
}

TEST(FsGridTest, getTaskForGlobalID1) {
   const std::array<fsgrid::FsSize_t, 3> globalSize{11, 5, 1048};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, true, false};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 0;
   MPI_Comm_size(parentComm, &numProcs);

   auto grid = fsgrid::FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic,
                                                                    {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   constexpr auto id = 666;
   const auto task = grid.getTaskForGlobalID(id);
   printf("Task for id %d: %d\n", id, task);
   ASSERT_EQ(0, task);
}

TEST(FsGridTest, getTaskForGlobalID2) {
   const std::array<fsgrid::FsSize_t, 3> globalSize{11, 5, 1048};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, true, false};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 4;

   auto grid = fsgrid::FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic,
                                                                    {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   constexpr auto id = 666;
   const auto task = grid.getTaskForGlobalID(id);
   printf("Task for id %d: %d\n", id, task);
   ASSERT_EQ(0, task);
}
