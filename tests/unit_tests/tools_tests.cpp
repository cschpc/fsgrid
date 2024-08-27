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
