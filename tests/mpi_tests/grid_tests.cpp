#include "grid.hpp"
#include <gtest/gtest.h>
#include <mpi.h>
#include <tools.hpp>

TEST(FsGridTest, localToGlobalRoundtrip1) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{1024, 666, 71};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, true, false};
   auto numProcs = 0;
   MPI_Comm_size(parentComm, &numProcs);

   const auto grid =
       FsGrid<std::array<double, 15>, 1>(globalSize, parentComm, numProcs, periodic, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
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

   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{6547, 16, 77};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, false, false};
   auto numProcs = 0;
   MPI_Comm_size(parentComm, &numProcs);

   const auto grid =
       FsGrid<std::array<double, 6>, 1>(globalSize, parentComm, numProcs, periodic, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
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

TEST(FsGridTest, localIdInBounds) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{647, 1, 666};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, false, true};
   auto numProcs = 0;
   MPI_Comm_size(parentComm, &numProcs);

   const auto grid =
       FsGrid<std::array<double, 50>, 1>(globalSize, parentComm, numProcs, periodic, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   for (int32_t x = 0; x < localSize[0]; x++) {
      for (int32_t y = 0; y < localSize[1]; y++) {
         for (int32_t z = 0; z < localSize[2]; z++) {
            const auto lid = grid.localIDFromLocalCoordinates(x, y, z);
            ASSERT_TRUE(grid.localIdInBounds(lid));
         }
      }
   }
}

TEST(FsGridTest, getNonPeriodic) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{12, 6, 2048};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{false, false, false};
   constexpr int32_t numGhostCells = 1;
   auto numProcs = 0;
   MPI_Comm_size(parentComm, &numProcs);

   auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic, {0.0, 0.0, 0.0},
                                                            {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   for (int32_t x = 0; x < localSize[0]; x++) {
      for (int32_t y = 0; y < localSize[1]; y++) {
         for (int32_t z = 0; z < localSize[2]; z++) {
            ASSERT_NE(grid.get(x, y, z), nullptr);
         }
      }
   }

   ASSERT_EQ(grid.get(-numGhostCells, 0, 0), nullptr);
   ASSERT_EQ(grid.get(-numGhostCells - 1, 0, 0), nullptr);
   ASSERT_EQ(grid.get(grid.getLocalSize()[0] + numGhostCells, 0, 0), nullptr);
   ASSERT_EQ(grid.get(grid.getLocalSize()[0] + numGhostCells - 1, 0, 0), nullptr);

   ASSERT_EQ(grid.get(0, -numGhostCells, 0), nullptr);
   ASSERT_EQ(grid.get(0, -numGhostCells - 1, 0), nullptr);
   ASSERT_EQ(grid.get(0, grid.getLocalSize()[1] + numGhostCells, 0), nullptr);
   ASSERT_EQ(grid.get(0, grid.getLocalSize()[1] + numGhostCells - 1, 0), nullptr);

   // This depends on the position on the grid
   if (grid.getLocalStart()[2] == 0) {
      ASSERT_EQ(grid.get(0, 0, -numGhostCells), nullptr);
   } else {
      ASSERT_NE(grid.get(0, 0, -numGhostCells), nullptr);
   }

   ASSERT_EQ(grid.get(0, 0, -numGhostCells - 1), nullptr);

   // This depends on the position on the grid
   if (grid.getLocalStart()[2] + grid.getLocalSize()[2] == static_cast<fsgrid_tools::FsIndex_t>(globalSize[2])) {
      ASSERT_EQ(grid.get(0, 0, grid.getLocalSize()[2] + numGhostCells - 1), nullptr);
   } else {
      ASSERT_NE(grid.get(0, 0, grid.getLocalSize()[2] + numGhostCells - 1), nullptr);
   }
   ASSERT_EQ(grid.get(0, 0, grid.getLocalSize()[2] + numGhostCells), nullptr);
}

TEST(FsGridTest, getPeriodic) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{120, 5, 1048};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, true, true};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 0;
   MPI_Comm_size(parentComm, &numProcs);

   auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic, {0.0, 0.0, 0.0},
                                                            {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   for (int32_t x = 0; x < localSize[0]; x++) {
      for (int32_t y = 0; y < localSize[1]; y++) {
         for (int32_t z = 0; z < localSize[2]; z++) {
            ASSERT_NE(grid.get(x, y, z), nullptr);
         }
      }
   }

   ASSERT_NE(grid.get(-numGhostCells, 0, 0), nullptr);
   ASSERT_EQ(grid.get(-numGhostCells - 1, 0, 0), nullptr);
   ASSERT_EQ(grid.get(grid.getLocalSize()[0] + numGhostCells, 0, 0), nullptr);
   ASSERT_NE(grid.get(grid.getLocalSize()[0] + numGhostCells - 1, 0, 0), nullptr);

   ASSERT_NE(grid.get(0, -numGhostCells, 0), nullptr);
   ASSERT_EQ(grid.get(0, -numGhostCells - 1, 0), nullptr);
   ASSERT_EQ(grid.get(0, grid.getLocalSize()[1] + numGhostCells, 0), nullptr);
   ASSERT_NE(grid.get(0, grid.getLocalSize()[1] + numGhostCells - 1, 0), nullptr);

   ASSERT_NE(grid.get(0, 0, -numGhostCells), nullptr);
   ASSERT_EQ(grid.get(0, 0, -numGhostCells - 1), nullptr);
   ASSERT_NE(grid.get(0, 0, grid.getLocalSize()[2] + numGhostCells - 1), nullptr);
   ASSERT_EQ(grid.get(0, 0, grid.getLocalSize()[2] + numGhostCells), nullptr);
}

TEST(FsGridTest, getTaskForGlobalID1) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{11, 5, 1048};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, true, false};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 0;
   MPI_Comm_size(parentComm, &numProcs);

   auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic, {0.0, 0.0, 0.0},
                                                            {0.0, 0.0, 0.0});
   constexpr auto id = 666;
   const auto task = grid.getTaskForGlobalID(id);
   printf("Task for id %d: %d\n", id, task);
   ASSERT_EQ(0, task);
}

TEST(FsGridTest, getTaskForGlobalID2) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{11, 5, 1048};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, true, false};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 4;

   auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic, {0.0, 0.0, 0.0},
                                                            {0.0, 0.0, 0.0});
   constexpr auto id = 666;
   const auto task = grid.getTaskForGlobalID(id);
   printf("Task for id %d: %d\n", id, task);
   ASSERT_EQ(0, task);
}

TEST(FsGridTest, shiftMultiplierNonperiodicXSplitOverX) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{1048, 5, 11};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{false, true, true};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 8;

   const auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic,
                                                                  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   const std::array xs{-numGhostCells, 0, localSize[0] + numGhostCells - 1};
   const std::array ys{-numGhostCells, 0, localSize[1] + numGhostCells - 1};
   const std::array zs{-numGhostCells, 0, localSize[2] + numGhostCells - 1};

   // clang-format off
   constexpr std::array values = {
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       1, 1, 1,
       1, 0, 1,
       1, 1, 1,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
   };
   // clang-format on

   auto i = 0ul;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            [[maybe_unused]] const auto a = grid.shiftMultiplier(x, y, z);
            if (grid.getRank() == 0) {
               printf("%d\n", a);
            }

            if (grid.getRank() != -1) {
               ASSERT_EQ(a, values[i++]);
            }
         }
      }
   }
}

TEST(FsGridTest, shiftMultiplierNonperiodicYSplitOverX) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{1048, 5, 11};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, false, true};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 8;

   const auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic,
                                                                  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   const std::array xs{-numGhostCells, 0, localSize[0] + numGhostCells - 1};
   const std::array ys{-numGhostCells, 0, localSize[1] + numGhostCells - 1};
   const std::array zs{-numGhostCells, 0, localSize[2] + numGhostCells - 1};

   // clang-format off
   constexpr std::array values = {
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       1, 0, 1,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
   };
   // clang-format on

   auto i = 0ul;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            [[maybe_unused]] const auto a = grid.shiftMultiplier(x, y, z);
            if (grid.getRank() == 0) {
               printf("%d\n", a);
            }

            if (grid.getRank() != -1) {
               ASSERT_EQ(a, values[i++]);
            }
         }
      }
   }
}

TEST(FsGridTest, shiftMultiplierNonperiodicZSplitOverX) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{1048, 5, 11};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, true, false};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 8;

   const auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic,
                                                                  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   const std::array xs{-numGhostCells, 0, localSize[0] + numGhostCells - 1};
   const std::array ys{-numGhostCells, 0, localSize[1] + numGhostCells - 1};
   const std::array zs{-numGhostCells, 0, localSize[2] + numGhostCells - 1};

   // clang-format off
   constexpr std::array values = {
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 1, 0,
       0, 0, 0,
       0, 1, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
   };
   // clang-format on

   auto i = 0ul;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            [[maybe_unused]] const auto a = grid.shiftMultiplier(x, y, z);
            if (grid.getRank() == 0) {
               printf("%d\n", a);
            }

            if (grid.getRank() != -1) {
               ASSERT_EQ(a, values[i++]);
            }
         }
      }
   }
}

TEST(FsGridTest, shiftMultiplierNonperiodicXSplitOverY) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{5, 1048, 11};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{false, true, true};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 8;

   const auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic,
                                                                  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   const std::array xs{-numGhostCells, 0, localSize[0] + numGhostCells - 1};
   const std::array ys{-numGhostCells, 0, localSize[1] + numGhostCells - 1};
   const std::array zs{-numGhostCells, 0, localSize[2] + numGhostCells - 1};

   // clang-format off
   constexpr std::array values = {
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       1, 0, 1,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
   };
   // clang-format on

   auto i = 0ul;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            [[maybe_unused]] const auto a = grid.shiftMultiplier(x, y, z);
            if (grid.getRank() == 0) {
               printf("%d\n", a);
            }

            if (grid.getRank() != -1) {
               ASSERT_EQ(a, values[i++]);
            }
         }
      }
   }
}

TEST(FsGridTest, shiftMultiplierNonperiodicYSplitOverY) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{5, 1048, 11};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, false, true};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 8;

   const auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic,
                                                                  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   const std::array xs{-numGhostCells, 0, localSize[0] + numGhostCells - 1};
   const std::array ys{-numGhostCells, 0, localSize[1] + numGhostCells - 1};
   const std::array zs{-numGhostCells, 0, localSize[2] + numGhostCells - 1};

   // clang-format off
   constexpr std::array values = {
       0, 0, 0,
       1, 1, 1,
       0, 0, 0,
       0, 0, 0,
       1, 0, 1,
       0, 0, 0,
       0, 0, 0,
       1, 1, 1,
       0, 0, 0,
   };
   // clang-format on

   auto i = 0ul;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            [[maybe_unused]] const auto a = grid.shiftMultiplier(x, y, z);
            if (grid.getRank() == 0) {
               printf("%d\n", a);
            }

            if (grid.getRank() != -1) {
               ASSERT_EQ(a, values[i++]);
            }
         }
      }
   }
}

TEST(FsGridTest, shiftMultiplierNonperiodicZSplitOverY) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{5, 1048, 11};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, true, false};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 8;

   const auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic,
                                                                  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   const std::array xs{-numGhostCells, 0, localSize[0] + numGhostCells - 1};
   const std::array ys{-numGhostCells, 0, localSize[1] + numGhostCells - 1};
   const std::array zs{-numGhostCells, 0, localSize[2] + numGhostCells - 1};

   // clang-format off
   constexpr std::array values = {
       0, 0, 0,
       0, 1, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 1, 0,
       0, 0, 0,
   };
   // clang-format on

   auto i = 0ul;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            [[maybe_unused]] const auto a = grid.shiftMultiplier(x, y, z);
            if (grid.getRank() == 0) {
               printf("%d\n", a);
            }

            if (grid.getRank() != -1) {
               ASSERT_EQ(a, values[i++]);
            }
         }
      }
   }
}

TEST(FsGridTest, shiftMultiplierNonperiodicXSplitOverZ) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{11, 5, 1048};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{false, true, true};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 8;

   const auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic,
                                                                  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   const std::array xs{-numGhostCells, 0, localSize[0] + numGhostCells - 1};
   const std::array ys{-numGhostCells, 0, localSize[1] + numGhostCells - 1};
   const std::array zs{-numGhostCells, 0, localSize[2] + numGhostCells - 1};

   // clang-format off
   constexpr std::array values = {
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 1, 0,
       0, 0, 0,
       0, 1, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
   };
   // clang-format on

   auto i = 0ul;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            [[maybe_unused]] const auto a = grid.shiftMultiplier(x, y, z);
            if (grid.getRank() == 0) {
               printf("%d\n", a);
            }

            if (grid.getRank() != -1) {
               ASSERT_EQ(a, values[i++]);
            }
         }
      }
   }
}

TEST(FsGridTest, shiftMultiplierNonperiodicYSplitOverZ) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{11, 5, 1048};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, false, true};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 8;

   const auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic,
                                                                  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   const std::array xs{-numGhostCells, 0, localSize[0] + numGhostCells - 1};
   const std::array ys{-numGhostCells, 0, localSize[1] + numGhostCells - 1};
   const std::array zs{-numGhostCells, 0, localSize[2] + numGhostCells - 1};

   // clang-format off
   constexpr std::array values = {
       0, 0, 0,
       0, 1, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 0, 0,
       0, 1, 0,
       0, 0, 0,
   };
   // clang-format on

   auto i = 0ul;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            [[maybe_unused]] const auto a = grid.shiftMultiplier(x, y, z);
            if (grid.getRank() == 0) {
               printf("%d\n", a);
            }

            if (grid.getRank() != -1) {
               ASSERT_EQ(a, values[i++]);
            }
         }
      }
   }
}

TEST(FsGridTest, shiftMultiplierNonperiodicZSplitOverZ) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{11, 5, 1048};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, true, false};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 8;

   const auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic,
                                                                  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
   const auto localSize = grid.getLocalSize();
   const std::array xs{-numGhostCells, 0, localSize[0] + numGhostCells - 1};
   const std::array ys{-numGhostCells, 0, localSize[1] + numGhostCells - 1};
   const std::array zs{-numGhostCells, 0, localSize[2] + numGhostCells - 1};

   // clang-format off
   constexpr std::array values = {
       0, 1, 0,
       0, 1, 0,
       0, 1, 0,
       0, 1, 0,
       0, 0, 0,
       0, 1, 0,
       0, 1, 0,
       0, 1, 0,
       0, 1, 0,
   };
   // clang-format on

   auto i = 0ul;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            [[maybe_unused]] const auto a = grid.shiftMultiplier(x, y, z);
            if (grid.getRank() == 0) {
               printf("%d\n", a);
            }

            if (grid.getRank() != -1) {
               ASSERT_EQ(a, values[i++]);
            }
         }
      }
   }
}

TEST(FsGridTest, localIDFromCellCoordinatesNonperiodicSplitOverX) {
   const std::array<fsgrid_tools::FsSize_t, 3> globalSize{1048, 11, 5};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{false, false, false};
   constexpr int32_t numGhostCells = 2;
   auto numProcs = 8;

   auto grid = FsGrid<std::array<double, 8>, numGhostCells>(globalSize, parentComm, numProcs, periodic, {0.0, 0.0, 0.0},
                                                            {0.0, 0.0, 0.0});

   constexpr auto value = std::numeric_limits<fsgrid_tools::LocalID>::min();

   const auto localSize = grid.getLocalSize();
   const std::array xs{-numGhostCells, 0, localSize[0] + numGhostCells - 1};
   const std::array ys{-numGhostCells, 0, localSize[1] + numGhostCells - 1};
   const std::array zs{-numGhostCells, 0, localSize[2] + numGhostCells - 1};
   const auto rank = grid.getRank();

   // clang-format off
   constexpr std::array valuesFirst = {
       value, value, value,
       value, value, value,
       value, value, value,
       value, value, value,
       value, 135l * 15l * 2l + 135l * 2l + 2l, value,
       value, value, value,
       value, value, value,
       value, 135l * 15l * 2l + 135l * 2l + 134l, value,
       value, value, value,
   };

   constexpr std::array valuesRest = {
       value, value, value,
       value, 135l * 15l * 2l + 135l * 2l, value,
       value, value, value,
       value, value, value,
       value, 135l * 15l * 2l + 135l * 2l + 2l, value,
       value, value, value,
       value, value, value,
       value, 135l * 15l * 2l + 135l * 2l + 134l, value,
       value, value, value,
   };

   constexpr std::array valuesLast = {
       value, value, value,
       value, 135l * 15l * 2l + 135l * 2l, value,
       value, value, value,
       value, value, value,
       value, 135l * 15l * 2l + 135l * 2l + 2l, value,
       value, value, value,
       value, value, value,
       value, value, value,
       value, value, value,
   };
   // clang-format on

   auto i = 0ul;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            [[maybe_unused]] const auto a = grid.localIDFromCellCoordinates(x, y, z);
            if (rank != -1) {
               if (rank == 0) {
                  ASSERT_EQ(a, valuesFirst[i++]);
               } else if (rank == numProcs - 1) {
                  ASSERT_EQ(a, valuesLast[i++]);
               } else {
                  ASSERT_EQ(a, valuesRest[i++]);
               }
            }
         }
      }
   }
}
