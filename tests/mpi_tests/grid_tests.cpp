#include "grid.hpp"
#include <fstream>
#include <gtest/gtest.h>
#include <mpi.h>
#include <string>
#include <tools.hpp>
#include <utility>

std::pair<FsGrid<std::array<double, 15>, 1>, std::string>
makeFsGridAndFilename(std::array<FsGridTools::FsSize_t, 3> globalSize, MPI_Comm parentComm,
                      std::array<bool, 3> isPeriodic,
                      const std::array<FsGridTools::Task_t, 3>& decomposition = {0, 0, 0}, bool verbose = false) {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   std::stringstream ss;
   ss << "fsgrid_" << globalSize[0] << "_" << globalSize[1] << "_" << globalSize[2] << "_" << isPeriodic[0] << "_"
      << isPeriodic[1] << "_" << isPeriodic[2] << "_" << decomposition[0] << "_" << decomposition[1] << "_"
      << decomposition[2] << "_" << verbose << "_display_rank_" << rank << ".txt";

   return std::make_pair(FsGrid<std::array<double, 15>, 1>(globalSize, parentComm, isPeriodic, decomposition, verbose),
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

void generateReferenceString(const std::string& filename, const std::string& gridStr) {
   const auto path = getPrefix() + filename;
   std::ofstream file(path);
   if (file.is_open()) {
      file << gridStr;
   }
}

TEST(FsGridTest, compareConstructedFsGridDisplayStringToReference1) {
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   const auto [grid, filename] =
       makeFsGridAndFilename({1024, 512, 64}, MPI_COMM_WORLD, {false, false, false}, {1, 1, size}, false);
   const auto [success, refStr] = readRefStr(filename);
   ASSERT_TRUE(success) << refStr;

   const auto str = grid.display();
   ASSERT_EQ(0, refStr.compare(str)) << "Reference string\n"
                                     << refStr << "\nand generated string\n"
                                     << str << "are not equal";
}

TEST(FsGridTest, compareConstructedFsGridDisplayStringToReference2) {
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   const auto [grid, filename] =
       makeFsGridAndFilename({512, 128, 1024}, MPI_COMM_WORLD, {true, false, true}, {1, size, 1}, false);
   const auto [success, refStr] = readRefStr(filename);
   ASSERT_TRUE(success) << refStr;

   const auto str = grid.display();
   ASSERT_EQ(0, refStr.compare(str)) << "Reference string\n"
                                     << refStr << "\nand generated string\n"
                                     << str << "are not equal";
}

TEST(FsGridTest, compareConstructedFsGridDisplayStringToReference3) {
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   const auto [grid, filename] =
       makeFsGridAndFilename({16, 512, 2048}, MPI_COMM_WORLD, {false, false, false}, {size, 1, 1}, false);
   const auto [success, refStr] = readRefStr(filename);
   ASSERT_TRUE(success) << refStr;

   const auto str = grid.display();
   ASSERT_EQ(0, refStr.compare(str)) << "Reference string\n"
                                     << refStr << "\nand generated string\n"
                                     << str << "are not equal";
}
