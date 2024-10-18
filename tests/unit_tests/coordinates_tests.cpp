#include "coordinates.hpp"

#include <gtest/gtest.h>

TEST(CoordinatesTest, singleRankCoordinates) {
   constexpr std::array<double, 3> physicalGridSpacing{0.1, 0.1, 0.1};
   constexpr std::array<double, 3> physicalGlobalStart{0.0, 0.0, 0.0};
   constexpr std::array<FsGridTools::FsSize_t, 3> globalSize{1024, 1, 512};
   constexpr std::array<bool, 3> periodic{false, false, false};
   constexpr std::array<FsGridTools::Task_t, 3> decomposition{0, 0, 0};
   constexpr std::array<FsGridTools::Task_t, 3> taskPosition{0, 0, 0};
   constexpr int32_t numRanks = 1;
   constexpr int32_t numGhostCells = 1;

   constexpr Coordinates coordinates(physicalGridSpacing, physicalGlobalStart, globalSize, periodic, decomposition,
                                     taskPosition, numRanks, numGhostCells);

   ASSERT_EQ(coordinates.numTasksPerDim[0], 1);
   ASSERT_EQ(coordinates.numTasksPerDim[1], 1);
   ASSERT_EQ(coordinates.numTasksPerDim[2], 1);

   ASSERT_EQ(coordinates.localStart[0], 0);
   ASSERT_EQ(coordinates.localStart[1], 0);
   ASSERT_EQ(coordinates.localStart[2], 0);

   ASSERT_EQ(coordinates.localSize[0], coordinates.globalSize[0]);
   ASSERT_EQ(coordinates.localSize[1], coordinates.globalSize[1]);
   ASSERT_EQ(coordinates.localSize[2], coordinates.globalSize[2]);

   auto i = 0;
   for (auto z = 0; z < coordinates.localSize[2]; z++) {
      for (auto y = 0; y < coordinates.localSize[1]; y++) {
         for (auto x = 0; x < coordinates.localSize[0]; x++) {
            const auto global = coordinates.localToGlobal(x, y, z);
            ASSERT_EQ(global[0], x);
            ASSERT_EQ(global[1], y);
            ASSERT_EQ(global[2], z);

            const auto local = coordinates.globalToLocal(global[0], global[1], global[2]);
            ASSERT_EQ(local[0], x);
            ASSERT_EQ(local[1], y);
            ASSERT_EQ(local[2], z);

            const auto physical = coordinates.getPhysicalCoords(x, y, z);
            ASSERT_EQ(physical[0], x * physicalGridSpacing[0]);
            ASSERT_EQ(physical[1], y * physicalGridSpacing[1]);
            ASSERT_EQ(physical[2], z * physicalGridSpacing[2]);

            ASSERT_EQ(coordinates.globalIDFromLocalCoordinates(x, y, z), i++);
         }
      }
   }
}

TEST(CoordinatesTest, neighbourIndex) {
   constexpr std::array<double, 3> physicalGridSpacing{0.1, 0.1, 0.1};
   constexpr std::array<double, 3> physicalGlobalStart{0.0, 0.0, 0.0};
   constexpr std::array<FsGridTools::FsSize_t, 3> globalSize{1024, 20, 512};
   constexpr std::array<bool, 3> periodic{false, false, false};
   constexpr std::array<FsGridTools::Task_t, 3> decomposition{0, 0, 0};
   constexpr std::array<FsGridTools::Task_t, 3> taskPosition{0, 0, 0};
   constexpr int32_t numRanks = 16;
   constexpr int32_t numGhostCells = 1;

   constexpr Coordinates coordinates(physicalGridSpacing, physicalGlobalStart, globalSize, periodic, decomposition,
                                     taskPosition, numRanks, numGhostCells);

   constexpr std::array xs{-1, 0, coordinates.localSize[0] + 1};
   constexpr std::array ys{-1, 0, coordinates.localSize[1] + 1};
   constexpr std::array zs{-1, 0, coordinates.localSize[2] + 1};

   size_t i = 0;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            ASSERT_EQ(coordinates.neighbourIndexFromCellCoordinates(x, y, z), i++);
         }
      }
   }
}

TEST(CoordinatesTest, shiftCellIndices) {
   constexpr std::array<double, 3> physicalGridSpacing{0.1, 0.1, 0.1};
   constexpr std::array<double, 3> physicalGlobalStart{0.0, 0.0, 0.0};
   constexpr std::array<FsGridTools::FsSize_t, 3> globalSize{1024, 20, 512};
   constexpr std::array<bool, 3> periodic{false, false, false};
   constexpr std::array<FsGridTools::Task_t, 3> decomposition{0, 0, 0};
   constexpr std::array<FsGridTools::Task_t, 3> taskPosition{0, 0, 0};
   constexpr int32_t numRanks = 16;
   constexpr int32_t numGhostCells = 1;

   constexpr Coordinates coordinates(physicalGridSpacing, physicalGlobalStart, globalSize, periodic, decomposition,
                                     taskPosition, numRanks, numGhostCells);

   constexpr std::array xs{-1, 0, coordinates.localSize[0] + 1};
   constexpr std::array ys{-1, 0, coordinates.localSize[1] + 1};
   constexpr std::array zs{-1, 0, coordinates.localSize[2] + 1};
   constexpr std::array<std::array<FsGridTools::FsIndex_t, 3>, 27> values{
       std::array{
           coordinates.localSize[0] - 1,
           coordinates.localSize[1] - 1,
           coordinates.localSize[2] - 1,
       },
       {
           coordinates.localSize[0] - 1,
           coordinates.localSize[1] - 1,
           0,
       },
       {
           coordinates.localSize[0] - 1,
           coordinates.localSize[1] - 1,
           1,
       },
       {
           coordinates.localSize[0] - 1,
           0,
           coordinates.localSize[2] - 1,
       },
       {
           coordinates.localSize[0] - 1,
           0,
           0,
       },
       {
           coordinates.localSize[0] - 1,
           0,
           1,
       },
       {
           coordinates.localSize[0] - 1,
           1,
           coordinates.localSize[2] - 1,
       },
       {
           coordinates.localSize[0] - 1,
           1,
           0,
       },
       {
           coordinates.localSize[0] - 1,
           1,
           1,
       },
       {
           0,
           coordinates.localSize[1] - 1,
           coordinates.localSize[2] - 1,
       },
       {
           0,
           coordinates.localSize[1] - 1,
           0,
       },
       {
           0,
           coordinates.localSize[1] - 1,
           1,
       },
       {
           0,
           0,
           coordinates.localSize[2] - 1,
       },
       {
           0,
           0,
           0,
       },
       {
           0,
           0,
           1,
       },
       {
           0,
           1,
           coordinates.localSize[2] - 1,
       },
       {
           0,
           1,
           0,
       },
       {
           0,
           1,
           1,
       },
       {
           1,
           coordinates.localSize[1] - 1,
           coordinates.localSize[2] - 1,
       },
       {
           1,
           coordinates.localSize[1] - 1,
           0,
       },
       {
           1,
           coordinates.localSize[1] - 1,
           1,
       },
       {
           1,
           0,
           coordinates.localSize[2] - 1,
       },
       {
           1,
           0,
           0,
       },
       {
           1,
           0,
           1,
       },
       {
           1,
           1,
           coordinates.localSize[2] - 1,
       },
       {
           1,
           1,
           0,
       },
       {
           1,
           1,
           1,
       },
   };

   size_t i = 0;
   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            const auto shifted = coordinates.shiftCellIndices(x, y, z);
            ASSERT_EQ(shifted[0], values[i][0]);
            ASSERT_EQ(shifted[1], values[i][1]);
            ASSERT_EQ(shifted[2], values[i][2]);
            i++;
         }
      }
   }
}

TEST(CoordinatesTest, indicesWithinDomain1) {
   constexpr std::array<double, 3> physicalGridSpacing{0.1, 0.1, 0.1};
   constexpr std::array<double, 3> physicalGlobalStart{0.0, 0.0, 0.0};
   constexpr std::array<FsGridTools::FsSize_t, 3> globalSize{1024, 20, 512};
   constexpr std::array<bool, 3> periodic{false, false, false};
   constexpr std::array<FsGridTools::Task_t, 3> decomposition{0, 0, 0};
   constexpr std::array<FsGridTools::Task_t, 3> taskPosition{0, 0, 0};
   constexpr int32_t numRanks = 16;
   constexpr int32_t numGhostCells = 1;

   constexpr Coordinates coordinates(physicalGridSpacing, physicalGlobalStart, globalSize, periodic, decomposition,
                                     taskPosition, numRanks, numGhostCells);

   constexpr std::array xs{-numGhostCells, 0, coordinates.localSize[0] + numGhostCells - 1};
   constexpr std::array ys{-numGhostCells, 0, coordinates.localSize[1] + numGhostCells - 1};
   constexpr std::array zs{-numGhostCells, 0, coordinates.localSize[2] + numGhostCells - 1};

   for (auto x : xs) {
      for (auto y : ys) {
         for (auto z : zs) {
            ASSERT_TRUE(coordinates.cellIndicesAreWithinBounds(x, y, z));
         }
      }
   }
}

TEST(CoordinatesTest, indicesWithinPeriodicDomain) {
   constexpr std::array<double, 3> physicalGridSpacing{0.1, 0.1, 0.1};
   constexpr std::array<double, 3> physicalGlobalStart{0.0, 0.0, 0.0};
   constexpr std::array<FsGridTools::FsSize_t, 3> globalSize{1024, 1, 512};
   constexpr std::array<bool, 3> periodic{false, true, false};
   constexpr std::array<FsGridTools::Task_t, 3> decomposition{0, 0, 0};
   constexpr std::array<FsGridTools::Task_t, 3> taskPosition{0, 0, 0};
   constexpr int32_t numRanks = 16;
   constexpr int32_t numGhostCells = 1;

   constexpr Coordinates coordinates(physicalGridSpacing, physicalGlobalStart, globalSize, periodic, decomposition,
                                     taskPosition, numRanks, numGhostCells);
   ASSERT_FALSE(coordinates.cellIndicesAreWithinBounds(0, -numGhostCells - 1, 0));
   ASSERT_TRUE(coordinates.cellIndicesAreWithinBounds(0, -numGhostCells, 0));
   ASSERT_TRUE(coordinates.cellIndicesAreWithinBounds(0, 0, 0));
   ASSERT_TRUE(coordinates.cellIndicesAreWithinBounds(0, numGhostCells, 0));
   ASSERT_FALSE(coordinates.cellIndicesAreWithinBounds(0, numGhostCells + 1, 0));
}

TEST(CoordinatesTest, indicesNotWithinDomain) {
   constexpr std::array<double, 3> physicalGridSpacing{0.1, 0.1, 0.1};
   constexpr std::array<double, 3> physicalGlobalStart{0.0, 0.0, 0.0};
   constexpr std::array<FsGridTools::FsSize_t, 3> globalSize{1024, 1, 512};
   constexpr std::array<bool, 3> periodic{false, false, false};
   constexpr std::array<FsGridTools::Task_t, 3> decomposition{0, 0, 0};
   constexpr std::array<FsGridTools::Task_t, 3> taskPosition{0, 0, 0};
   constexpr int32_t numRanks = 16;
   constexpr int32_t numGhostCells = 1;

   constexpr Coordinates coordinates(physicalGridSpacing, physicalGlobalStart, globalSize, periodic, decomposition,
                                     taskPosition, numRanks, numGhostCells);
   ASSERT_FALSE(coordinates.cellIndicesAreWithinBounds(0, -numGhostCells - 1, 0));
   ASSERT_FALSE(coordinates.cellIndicesAreWithinBounds(0, -numGhostCells, 0));
   ASSERT_TRUE(coordinates.cellIndicesAreWithinBounds(0, 0, 0));
   ASSERT_FALSE(coordinates.cellIndicesAreWithinBounds(0, numGhostCells, 0));
   ASSERT_FALSE(coordinates.cellIndicesAreWithinBounds(0, numGhostCells + 1, 0));
}

TEST(CoordinatesTest, indicesWithinNonPeriodicDomain) {
   constexpr std::array<double, 3> physicalGridSpacing{0.1, 0.1, 0.1};
   constexpr std::array<double, 3> physicalGlobalStart{0.0, 0.0, 0.0};
   constexpr std::array<FsGridTools::FsSize_t, 3> globalSize{1024, 2, 512};
   constexpr std::array<bool, 3> periodic{false, false, false};
   constexpr std::array<FsGridTools::Task_t, 3> decomposition{0, 0, 0};
   constexpr std::array<FsGridTools::Task_t, 3> taskPosition{0, 0, 0};
   constexpr int32_t numRanks = 16;
   constexpr int32_t numGhostCells = 1;

   constexpr Coordinates coordinates(physicalGridSpacing, physicalGlobalStart, globalSize, periodic, decomposition,
                                     taskPosition, numRanks, numGhostCells);
   ASSERT_FALSE(coordinates.cellIndicesAreWithinBounds(0, -numGhostCells - 1, 0));
   ASSERT_TRUE(coordinates.cellIndicesAreWithinBounds(0, -numGhostCells, 0));
   ASSERT_TRUE(coordinates.cellIndicesAreWithinBounds(0, 0, 0));
   ASSERT_TRUE(coordinates.cellIndicesAreWithinBounds(0, coordinates.localSize[1] + numGhostCells - 1, 0));
   ASSERT_FALSE(coordinates.cellIndicesAreWithinBounds(0, coordinates.localSize[1] + numGhostCells, 0));
}
