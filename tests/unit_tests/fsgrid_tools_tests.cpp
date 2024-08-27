#include "tools.hpp"

#include <gtest/gtest.h>

TEST(FsGridToolsTests, calcLocalStart1) {
    constexpr FsGridTools::FsSize_t numGlobalCells = 1024u;
    constexpr FsGridTools::Task_t numTasks = 32u;
    constexpr FsGridTools::Task_t taskIndex = 0;

    ASSERT_EQ(FsGridTools::calcLocalStart(numGlobalCells, numTasks, taskIndex), 0);
}

TEST(FsGridToolsTests, calcLocalStart2) {
   constexpr FsGridTools::FsSize_t numGlobalCells = 1024u;
   constexpr FsGridTools::Task_t numTasks = 32u;
   constexpr FsGridTools::Task_t taskIndex = 1;

   ASSERT_EQ(FsGridTools::calcLocalStart(numGlobalCells, numTasks, taskIndex), 32);
}

TEST(FsGridToolsTests, calcLocalStart3) {
   constexpr FsGridTools::FsSize_t numGlobalCells = 1024u;
   constexpr FsGridTools::Task_t numTasks = 32u;
   constexpr FsGridTools::Task_t taskIndex = 2;

   ASSERT_EQ(FsGridTools::calcLocalStart(numGlobalCells, numTasks, taskIndex), 64);
}

TEST(FsGridToolsTests, calcLocalStart2) {
   constexpr FsGridTools::FsSize_t numGlobalCells = 1024u;
   constexpr FsGridTools::Task_t numTasks = 32u;
   constexpr FsGridTools::Task_t taskIndex = 3;

   ASSERT_EQ(FsGridTools::calcLocalStart(numGlobalCells, numTasks, taskIndex), 96);
}
