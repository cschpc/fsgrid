#include "tools.hpp"

#include <gtest/gtest.h>

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
