#pragma once

/*
  Copyright (C) 2016 Finnish Meteorological Institute
  Copyright (C) 2016-2024 CSC -IT Center for Science

  This file is part of fsgrid

  fsgrid is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  fsgrid is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY;
  without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with fsgrid.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <array>
#include <cstddef>
#include <cstdint>

namespace fsgrid {
struct BitMask32 {
   constexpr BitMask32(uint32_t bits) : bits(bits) {}
   constexpr uint32_t operator[](uint32_t i) const {
      // Shifting by more than N - 1 is undefined behaviour for N bit values
      constexpr uint32_t n = sizeof(bits) * 8;
      const uint32_t mul = i < n;
      i &= n - 1;
      return mul * ((bits & (1u << i)) >> i);
   }

private:
   const uint32_t bits = 0u;
};

struct StencilConstants {
   const std::array<int32_t, 3> limits = {};
   const std::array<int32_t, 3> multipliers = {};
   const int32_t offset = 0;
   const BitMask32 shift = 0;
   const BitMask32 fallbackToCenter = 0;

   constexpr StencilConstants(const std::array<int32_t, 3>& limits, const std::array<int32_t, 3>& multipliers,
                              int32_t offset, BitMask32 shift, BitMask32 fallbackToCenter)
       : limits(limits), multipliers(multipliers), offset(offset), shift(shift), fallbackToCenter(fallbackToCenter) {}
   constexpr StencilConstants() {}
};

struct FsStencil {
   const int32_t i = 0;
   const int32_t j = 0;
   const int32_t k = 0;

private:
   const StencilConstants constants = {};

public:
   constexpr FsStencil(int32_t i, int32_t j, int32_t k, const StencilConstants& constants)
       : i(i), j(j), k(k), constants(constants) {}

   // clang-format off
   // These names come from the right hand rule, with
   // - x horizontal (left right)
   // - y vertical   (down up)
   // - z on the line of sight (far near)
   //
   //            up
   //             |
   //             |  far
   //             | /
   //             |/
   //  left ------o------ right
   //            /|
   //           / | 
   //       near  |
   //             |
   //            down

   constexpr size_t center() const          { return calculateIndex({i,     j    , k    }); }
   constexpr size_t near() const            { return calculateIndex({i,     j    , k + 1}); }
   constexpr size_t far()  const            { return calculateIndex({i,     j    , k - 1}); }

   constexpr size_t up()  const             { return calculateIndex({i,     j + 1, k    }); }
   constexpr size_t upnear() const          { return calculateIndex({i,     j + 1, k + 1}); }
   constexpr size_t upfar() const           { return calculateIndex({i,     j + 1, k - 1}); }

   constexpr size_t down() const            { return calculateIndex({i,     j - 1, k    }); }
   constexpr size_t downnear() const        { return calculateIndex({i,     j - 1, k + 1}); }
   constexpr size_t downfar() const         { return calculateIndex({i,     j - 1, k - 1}); }

   constexpr size_t right() const           { return calculateIndex({i + 1, j    , k    }); }
   constexpr size_t rightnear() const       { return calculateIndex({i + 1, j    , k + 1}); }
   constexpr size_t rightfar() const        { return calculateIndex({i + 1, j    , k - 1}); }

   constexpr size_t rightup() const         { return calculateIndex({i + 1, j + 1, k    }); }
   constexpr size_t rightupnear() const     { return calculateIndex({i + 1, j + 1, k + 1}); }
   constexpr size_t rightupfar() const      { return calculateIndex({i + 1, j + 1, k - 1}); }

   constexpr size_t rightdown() const       { return calculateIndex({i + 1, j - 1, k    }); }
   constexpr size_t rightdownnear() const   { return calculateIndex({i + 1, j - 1, k + 1}); }
   constexpr size_t rightdownfar() const    { return calculateIndex({i + 1, j - 1, k - 1}); }

   constexpr size_t left() const            { return calculateIndex({i - 1, j    , k    }); }
   constexpr size_t leftnear() const        { return calculateIndex({i - 1, j    , k + 1}); }
   constexpr size_t leftfar() const         { return calculateIndex({i - 1, j    , k - 1}); }

   constexpr size_t leftup() const          { return calculateIndex({i - 1, j + 1, k    }); }
   constexpr size_t leftupnear() const      { return calculateIndex({i - 1, j + 1, k + 1}); }
   constexpr size_t leftupfar() const       { return calculateIndex({i - 1, j + 1, k - 1}); }

   constexpr size_t leftdown() const        { return calculateIndex({i - 1, j - 1, k    }); }
   constexpr size_t leftdownnear() const    { return calculateIndex({i - 1, j - 1, k + 1}); }
   constexpr size_t leftdownfar() const     { return calculateIndex({i - 1, j - 1, k - 1}); }
   // clang-format on

   constexpr bool cellExists(int32_t io, int32_t jo, int32_t ko) const {
      const auto no = neighbourOffset({i + io, j + jo, k + ko});
      const auto ni = neighbourIndex(no);
      return static_cast<int32_t>(constants.fallbackToCenter[ni]) == 0;
   }

   constexpr std::array<size_t, 27> indices() const {
      // Return an array containing all indices
      // x changes fastest, then y, then z
      // clang-format off
       return {
           calculateIndex({i - 1, j - 1, k - 1}),
           calculateIndex({i    , j - 1, k - 1}),
           calculateIndex({i + 1, j - 1, k - 1}),
           calculateIndex({i - 1, j    , k - 1}),
           calculateIndex({i    , j    , k - 1}),
           calculateIndex({i + 1, j    , k - 1}),
           calculateIndex({i - 1, j + 1, k - 1}),
           calculateIndex({i    , j + 1, k - 1}),
           calculateIndex({i + 1, j + 1, k - 1}),
           calculateIndex({i - 1, j - 1, k    }),
           calculateIndex({i    , j - 1, k    }),
           calculateIndex({i + 1, j - 1, k    }),
           calculateIndex({i - 1, j    , k    }),
           calculateIndex({i    , j    , k    }),
           calculateIndex({i + 1, j    , k    }),
           calculateIndex({i - 1, j + 1, k    }),
           calculateIndex({i    , j + 1, k    }),
           calculateIndex({i + 1, j + 1, k    }),
           calculateIndex({i - 1, j - 1, k + 1}),
           calculateIndex({i    , j - 1, k + 1}),
           calculateIndex({i + 1, j - 1, k + 1}),
           calculateIndex({i - 1, j    , k + 1}),
           calculateIndex({i    , j    , k + 1}),
           calculateIndex({i + 1, j    , k + 1}),
           calculateIndex({i - 1, j + 1, k + 1}),
           calculateIndex({i    , j + 1, k + 1}),
           calculateIndex({i + 1, j + 1, k + 1}),
       };
      // clang-format on
   }

   constexpr size_t calculateIndex(std::array<int32_t, 3> cellIndex) const {
      const auto no = neighbourOffset(cellIndex);
      const auto ni = neighbourIndex(no);
      cellIndex = fallback(offsetValues(cellIndex, no, ni), ni);

      return applyMultipliersAndOffset(cellIndex);
   }

private:
   constexpr std::array<int32_t, 3> neighbourOffset(const std::array<int32_t, 3>& cellIndex) const {
      // clang-format off
      // A triplet of (-, 0, +) values, with 27 possibilities
      // The value for a coordinate is
      // - if the coordinate is below zero,
      // 0 if it's within bounds, or
      // + if it's at or above limit
      //
      // Visualized as 2D slices:
      // (values on the charts indicate (xyz) in order)
      //
      //     +Z plane
      // y
      // ^ -++  0++  +++
      // | -0+  00+  +0+
      // | --+  0-+  +-+
      // o-------------->x
      //
      //     0Z plane
      // y
      // ^ -+0  0+0  ++0
      // | -00  000  +00
      // | --0  0-0  +-0
      // o-------------->x
      //
      //     -Z plane
      // y
      // ^ -+-  0+-  ++-
      // | -0-  00-  +0-
      // | ---  0--  +--
      // o-------------->x
      //
      // clang-format on
      return {
          (cellIndex[0] >= constants.limits[0]) - (cellIndex[0] < 0),
          (cellIndex[1] >= constants.limits[1]) - (cellIndex[1] < 0),
          (cellIndex[2] >= constants.limits[2]) - (cellIndex[2] < 0),
      };
   }

   constexpr uint32_t neighbourIndex(const std::array<int32_t, 3>& no) const {
      // Translate a triplet of (-, 0, +) values to a single value between 0 and 26
      // - 0 is at (---) corner
      // - 13 is at (000) i.e. center
      // - 26 is at (+++) corner
      // - z changes fastest, then y, then x
      return static_cast<uint32_t>(13 + no[0] * 9 + no[1] * 3 + no[2]);
   }

   constexpr std::array<int32_t, 3> offsetValues(const std::array<int32_t, 3>& cellIndex,
                                                 const std::array<int32_t, 3>& no, uint32_t ni) const {
      // If the shift bit is 1 for neighbour index 'ni', add an offset to the given values
      const auto addOffset = static_cast<int32_t>(constants.shift[ni]);
      const auto offsets = shiftOffsets(no);
      return {
          cellIndex[0] + addOffset * offsets[0],
          cellIndex[1] + addOffset * offsets[1],
          cellIndex[2] + addOffset * offsets[2],
      };
   }

   constexpr std::array<int32_t, 3> shiftOffsets(const std::array<int32_t, 3>& no) const {
      // -limit, if the neihbour offset 'no' is +
      // 0, if the neihbour offset is 0
      // +limit, if neihbour offset is -
      return {
          -no[0] * constants.limits[0],
          -no[1] * constants.limits[1],
          -no[2] * constants.limits[2],
      };
   }

   constexpr std::array<int32_t, 3> fallback(const std::array<int32_t, 3>& cellIndex, uint32_t ni) const {
      // If the cellIndex is invalid, we'll use the center cell, i.e. (i, j, k)
      const auto invalid = static_cast<int32_t>(constants.fallbackToCenter[ni]);
      const auto valid = invalid ^ 1;
      return {
          valid * cellIndex[0] + invalid * i,
          valid * cellIndex[1] + invalid * j,
          valid * cellIndex[2] + invalid * k,
      };
   }

   constexpr size_t applyMultipliersAndOffset(const std::array<int32_t, 3>& cellIndex) const {
      // A dot product between cellIndex and constants.multipliers + an offset
      return static_cast<size_t>(constants.offset + constants.multipliers[0] * cellIndex[0] +
                                 constants.multipliers[1] * cellIndex[1] + constants.multipliers[2] * cellIndex[2]);
   }
};
} // namespace fsgrid
