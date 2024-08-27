#include "tools.hpp"

#include <gtest/gtest.h>
#include <limits>

TEST(FsGridToolsTests, calcLocalStart1) {
    constexpr FsGridTools::FsSize_t numGlobalCells = 1024u;
    constexpr FsGridTools::Task_t numTasks = 32u;

    ASSERT_EQ(FsGridTools::calcLocalStart(numGlobalCells, numTasks, 0), 0);
    ASSERT_EQ(FsGridTools::calcLocalStart(numGlobalCells, numTasks, 1), 32);
    ASSERT_EQ(FsGridTools::calcLocalStart(numGlobalCells, numTasks, 2), 64);
    ASSERT_EQ(FsGridTools::calcLocalStart(numGlobalCells, numTasks, 3), 96);
}

TEST(FsGridToolsTests, calcLocalStart2) {
   constexpr FsGridTools::FsSize_t numGlobalCells = 666u;
   constexpr FsGridTools::Task_t numTasks = 64u;

   for (int i = 0; i < 26; i++) {
      ASSERT_EQ(FsGridTools::calcLocalStart(numGlobalCells, numTasks, i), i * 11);
   }
   for (int i = 26; i < numTasks; i++) {
      ASSERT_EQ(FsGridTools::calcLocalStart(numGlobalCells, numTasks, i), i * 10 + 26);
   }
}

TEST(FsGridToolsTests, calcLocalSize1) {
   constexpr FsGridTools::FsSize_t numGlobalCells = 1024u;
   constexpr FsGridTools::Task_t numTasks = 32u;

   for (int i = 0; i < numTasks; i++) {
      ASSERT_EQ(FsGridTools::calcLocalSize(numGlobalCells, numTasks, i), 32);
   }
}

TEST(FsGridToolsTests, calcLocalSize2) {
   constexpr FsGridTools::FsSize_t numGlobalCells = 666u;
   constexpr FsGridTools::Task_t numTasks = 64u;

   for (int i = 0; i < 26; i++) {
      ASSERT_EQ(FsGridTools::calcLocalSize(numGlobalCells, numTasks, i), 11);
   }

   for (int i = 26; i < numTasks; i++) {
      ASSERT_EQ(FsGridTools::calcLocalSize(numGlobalCells, numTasks, i), 10);
   }
}

TEST(FsGridToolsTests, globalIDtoCellCoord1) {
   constexpr std::array<FsGridTools::FsSize_t, 3> globalSize = {3, 7, 45};

   for (FsGridTools::FsSize_t i = 0; i < globalSize[0] * globalSize[1] * globalSize[2]; i++) {
      const std::array<FsGridTools::FsIndex_t, 3> result = {
          (FsGridTools::FsIndex_t)(i % globalSize[0]),
          (FsGridTools::FsIndex_t)((i / globalSize[0]) % globalSize[1]),
          (FsGridTools::FsIndex_t)((i / (globalSize[0] * globalSize[1])) % globalSize[2]),
      };
      ASSERT_EQ(FsGridTools::globalIDtoCellCoord(i, globalSize), result);
   }
}

TEST(FsGridToolsTests, globalIDtoCellCoord_globalSize_would_overflow) {
   constexpr uint32_t maxUint = std::numeric_limits<uint32_t>::max();
   constexpr std::array<FsGridTools::FsSize_t, 3> globalSize = {maxUint, 1, 1};
   const std::array<FsGridTools::FsIndex_t, 3> result = {
       -2147483648,
       -2147483648,
       -2147483648,
   };
   ASSERT_EQ(FsGridTools::globalIDtoCellCoord(globalSize[0] - 1, globalSize), result);
}

TEST(FsGridToolsTests, globalIDtoCellCoord_globalSize0_is_maximum_int) {
   constexpr int32_t maxInt = std::numeric_limits<int32_t>::max();
   constexpr std::array<FsGridTools::FsSize_t, 3> globalSize = {maxInt, maxInt, maxInt};

   std::array<FsGridTools::FsIndex_t, 3> result = {maxInt - 1, 0, 0};
   ASSERT_EQ(FsGridTools::globalIDtoCellCoord(globalSize[0] - 1, globalSize), result);

   result = {0, 1, 0};
   ASSERT_EQ(FsGridTools::globalIDtoCellCoord(globalSize[0], globalSize), result);

   result = {0, 0, 1};
   ASSERT_EQ(FsGridTools::globalIDtoCellCoord((int64_t)globalSize[0] * globalSize[0], globalSize), result);

   result = {maxInt - 1, maxInt - 1, 0};
   ASSERT_EQ(FsGridTools::globalIDtoCellCoord((int64_t)globalSize[0] * globalSize[0] - 1, globalSize), result);
}
