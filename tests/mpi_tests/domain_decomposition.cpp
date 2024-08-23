#include "fsgrid.hpp"
#include <array>
#include <gtest/gtest.h>

struct Decomposition {
    int32_t x = 0u;
    int32_t y = 0u;
    int32_t z = 0u;
};

struct SystemSize {
    uint32_t x = 0u;
    uint32_t y = 0u;
    uint32_t z = 0u;
};

Decomposition computeDecomposition(const SystemSize systemSize,
                                   const uint32_t nProcs) {
    std::array<FsGridTools::Task_t, 3> dd;
    FsGridTools::computeDomainDecomposition(
        {
            systemSize.x,
            systemSize.y,
            systemSize.z,
        },
        nProcs, dd, 1, false);

    return Decomposition{dd[0], dd[1], dd[2]};
}

TEST(DomainDecomposition, size_256_256_256_nprocs_32) {
    const auto [x, y, z] = computeDecomposition(SystemSize{256, 256, 256}, 32);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 32);
}

TEST(DomainDecomposition, size_128_256_256_nprocs_32) {
    const auto [x, y, z] = computeDecomposition(SystemSize{128, 256, 256}, 32);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 32);
}

TEST(DomainDecomposition, size_256_128_256_nprocs_32) {
    const auto [x, y, z] = computeDecomposition(SystemSize{256, 128, 256}, 32);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 32);
}

TEST(DomainDecomposition, size_256_256_128_nprocs_32) {
    const auto [x, y, z] = computeDecomposition(SystemSize{256, 256, 128}, 32);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 32);
    ASSERT_EQ(z, 1);
}

TEST(DomainDecomposition, size_256_256_256_nprocs_1) {
    const auto [x, y, z] = computeDecomposition(SystemSize{256, 256, 256}, 1);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 1);
}

TEST(DomainDecomposition, size_128_256_256_nprocs_1) {
    const auto [x, y, z] = computeDecomposition(SystemSize{128, 256, 256}, 1);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 1);
}

TEST(DomainDecomposition, size_256_128_256_nprocs_1) {
    const auto [x, y, z] = computeDecomposition(SystemSize{256, 128, 256}, 1);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 1);
}

TEST(DomainDecomposition, size_256_256_128_nprocs_1) {
    const auto [x, y, z] = computeDecomposition(SystemSize{256, 256, 128}, 1);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 1);
}

TEST(DomainDecomposition, size_256_256_256_nprocs_64) {
    const auto [x, y, z] = computeDecomposition(SystemSize{256, 256, 256}, 64);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 64);
}

TEST(DomainDecomposition, size_128_256_256_nprocs_64) {
    const auto [x, y, z] = computeDecomposition(SystemSize{128, 256, 256}, 64);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 64);
}

TEST(DomainDecomposition, size_256_128_256_nprocs_64) {
    const auto [x, y, z] = computeDecomposition(SystemSize{256, 128, 256}, 64);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 64);
}

TEST(DomainDecomposition, size_256_256_128_nprocs_64) {
    const auto [x, y, z] = computeDecomposition(SystemSize{256, 256, 128}, 64);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 64);
    ASSERT_EQ(z, 1);
}

TEST(DomainDecomposition, size_1024_256_512_nprocs_64) {
    const auto [x, y, z] = computeDecomposition(SystemSize{1024, 256, 512}, 64);
    ASSERT_EQ(x, 64);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 1);
}

TEST(DomainDecomposition, size_256_512_128) {
    const auto [x, y, z] = computeDecomposition(SystemSize{256, 512, 128}, 64);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 64);
    ASSERT_EQ(z, 1);
}

TEST(DomainDecomposition, size_64_128_256_nprocs_64) {
    const auto [x, y, z] = computeDecomposition(SystemSize{64, 128, 256}, 64);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 64);
}

TEST(DomainDecomposition, size_64_256_1024_nprocs_64) {
    const auto [x, y, z] = computeDecomposition(SystemSize{64, 256, 1024}, 64);
    ASSERT_EQ(x, 1);
    ASSERT_EQ(y, 1);
    ASSERT_EQ(z, 64);
}
