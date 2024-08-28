#include "grid.hpp"
#include <cstdio>
#include <gtest/gtest.h>
#include <mpi.h>
#include <string>
#include <tools.hpp>

TEST(FsGridTest, displayGrid) {
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   const auto grid =
       FsGrid<std::array<double, 15>, 1>({1024, 512, 64}, MPI_COMM_WORLD, {false, false, false}, {1, 1, size});
   const auto str = grid.display();
   const auto filename = "grid_display_" + std::to_string(rank) + ".txt";
   if (const auto fstream = std::fopen(filename.c_str(), "w")) {
      std::fwrite(str.c_str(), sizeof str[0], str.size(), fstream);
      std::fclose(fstream);
   }
}
