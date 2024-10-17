#include "grid.hpp"
#include <fstream>
#include <gtest/gtest.h>
#include <mpi.h>
#include <sstream>
#include <string>
#include <tools.hpp>
#include <utility>

std::pair<FsGrid<std::array<double, 15>, 1>, std::string>
makeFsGridAndFilename(std::array<FsGridTools::FsSize_t, 3> globalSize, MPI_Comm parentComm,
                      std::array<bool, 3> isPeriodic,
                      const std::array<FsGridTools::Task_t, 3>& decomposition = {0, 0, 0}) {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   std::stringstream ss;
   ss << "fsgrid_" << globalSize[0] << "_" << globalSize[1] << "_" << globalSize[2] << "_" << isPeriodic[0] << "_"
      << isPeriodic[1] << "_" << isPeriodic[2] << "_" << decomposition[0] << "_" << decomposition[1] << "_"
      << decomposition[2] << "_display_rank_" << rank << ".txt";

   return std::make_pair(FsGrid<std::array<double, 15>, 1>(globalSize, parentComm, isPeriodic, {0.0, 0.0, 0.0},
                                                           {0.0, 0.0, 0.0}, decomposition),
                         ss.str());
}

std::string getPrefix() { return "../testdata/"; }

std::pair<bool, std::string> readRefStr(const std::string& filename) {
   const auto path = getPrefix() + filename;
   std::stringstream contents;
   std::ifstream file(path);

   if (file.is_open()) {
      contents << file.rdbuf();
      return std::make_pair(true, contents.str());
   } else {
      contents << "Could not read reference string with filename " << path;
      return std::make_pair(false, contents.str());
   }
}

template <typename T> std::string makeTestStringFromGrid(const T& grid) {
   std::stringstream ss;
   ss << grid.display() << "\n";
   return ss.str();
}

template <typename T> void generateReferenceString(const T& grid, const std::string& filename) {
   const auto path = getPrefix() + filename;
   std::ofstream file(path);
   if (file.is_open()) {
      file << makeTestStringFromGrid(grid);
   }
}

TEST(FsGridTest, compareConstructedFsGridDisplayStringToReference1) {
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   const auto [grid, filename] =
       makeFsGridAndFilename({1024, 512, 64}, MPI_COMM_WORLD, {false, false, false}, {1, 1, size});
   const auto [success, refStr] = readRefStr(filename);
   ASSERT_TRUE(success) << refStr;

   const auto str = makeTestStringFromGrid(grid);
   if (rank == 0) {
      ASSERT_EQ(0, refStr.compare(str)) << "Generated string\n"
                                        << str << "\nand ref string\n"
                                        << refStr << "are not equal";
   } else {
      ASSERT_EQ(0, refStr.compare(str));
   }
}

TEST(FsGridTest, compareConstructedFsGridDisplayStringToReference2) {
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   const auto [grid, filename] =
       makeFsGridAndFilename({512, 128, 1024}, MPI_COMM_WORLD, {true, false, true}, {1, size, 1});
   const auto [success, refStr] = readRefStr(filename);
   ASSERT_TRUE(success) << refStr;

   const auto str = makeTestStringFromGrid(grid);
   if (rank == 0) {
      ASSERT_EQ(0, refStr.compare(str)) << "Generated string\n"
                                        << str << "\nand ref string\n"
                                        << refStr << "are not equal";
   } else {
      ASSERT_EQ(0, refStr.compare(str));
   }
}

TEST(FsGridTest, compareConstructedFsGridDisplayStringToReference3) {
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   const auto [grid, filename] =
       makeFsGridAndFilename({16, 512, 2048}, MPI_COMM_WORLD, {false, false, false}, {size, 1, 1});
   const auto [success, refStr] = readRefStr(filename);
   ASSERT_TRUE(success) << refStr;

   const auto str = makeTestStringFromGrid(grid);
   if (rank == 0) {
      ASSERT_EQ(0, refStr.compare(str)) << "Generated string\n"
                                        << str << "\nand ref string\n"
                                        << refStr << "are not equal";
   } else {
      ASSERT_EQ(0, refStr.compare(str));
   }
}

TEST(FsGridTest, compareConstructedFsGridDisplayStringToReference4) {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   const auto [grid, filename] = makeFsGridAndFilename({16, 512, 2048}, MPI_COMM_WORLD, {false, false, false});
   const auto [success, refStr] = readRefStr(filename);
   ASSERT_TRUE(success) << refStr;

   const auto str = makeTestStringFromGrid(grid);
   if (rank == 0) {
      ASSERT_EQ(0, refStr.compare(str)) << "Generated string\n"
                                        << str << "\nand ref string\n"
                                        << refStr << "are not equal";
   } else {
      ASSERT_EQ(0, refStr.compare(str));
   }
}

TEST(FsGridTest, compareConstructedFsGridDisplayStringToReference5) {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   const auto [grid, filename] = makeFsGridAndFilename({1024, 666, 71}, MPI_COMM_WORLD, {true, true, false});
   const auto [success, refStr] = readRefStr(filename);
   ASSERT_TRUE(success) << refStr;

   const auto str = makeTestStringFromGrid(grid);
   if (rank == 0) {
      ASSERT_EQ(0, refStr.compare(str)) << "Generated string\n"
                                        << str << "\nand ref string\n"
                                        << refStr << "are not equal";
   } else {
      ASSERT_EQ(0, refStr.compare(str));
   }
}

TEST(FsGridTest, localToGlobalRoundtrip1) {
   const std::array<FsGridTools::FsSize_t, 3> globalSize{1024, 666, 71};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, true, false};

   const auto grid =
       FsGrid<std::array<double, 15>, 1>(globalSize, parentComm, periodic, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
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

   const std::array<FsGridTools::FsSize_t, 3> globalSize{6547, 16, 77};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, false, false};

   const auto grid =
       FsGrid<std::array<double, 6>, 1>(globalSize, parentComm, periodic, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
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
   const std::array<FsGridTools::FsSize_t, 3> globalSize{647, 1, 666};
   const MPI_Comm parentComm = MPI_COMM_WORLD;
   const std::array<bool, 3> periodic{true, false, true};

   const auto grid =
       FsGrid<std::array<double, 50>, 1>(globalSize, parentComm, periodic, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
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
