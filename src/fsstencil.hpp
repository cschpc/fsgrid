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

#include "tools.hpp"

#include <array>

namespace fsgrid {

// This is constructed during FsGrid construction time
struct StencilConstants {
   const std::array<int32_t, 3> limits = {};
   const std::array<int32_t, 3> multipliers = {};
   const int32_t offset = 0;
   const fsgrid_tools::BitMask32 shift = 0;
   const fsgrid_tools::BitMask32 fallbackToCenter = 0;

   constexpr StencilConstants(const std::array<int32_t, 3>& limits, const std::array<int32_t, 3>& multipliers,
                              int32_t offset, fsgrid_tools::BitMask32 shift, fsgrid_tools::BitMask32 fallbackToCenter)
       : limits(limits), multipliers(multipliers), offset(offset), shift(shift), fallbackToCenter(fallbackToCenter) {}
   constexpr StencilConstants() {}
};

// TODO: Add center(), left() etc. functions
struct FsStencil {
private:
   const std::array<int32_t, 3> centerCell = {};
   const StencilConstants constants = {};

public:
   constexpr FsStencil(const std::array<int32_t, 3>& centerCell, const StencilConstants& constants)
       : centerCell(centerCell), constants(constants) {}

   constexpr std::array<int32_t, 3> localityMultipliers(const std::array<int32_t, 3>& values) const {
      return {
          (values[0] >= constants.limits[0]) - (values[0] < 0),
          (values[1] >= constants.limits[1]) - (values[1] < 0),
          (values[2] >= constants.limits[2]) - (values[2] < 0),
      };
   }

   constexpr uint32_t neighbourIndex(const std::array<int32_t, 3>& values) const {
      return static_cast<uint32_t>(13 + values[0] * 9 + values[1] * 3 + values[2]);
   }

   constexpr std::array<int32_t, 3> shiftOffsets(const std::array<int32_t, 3>& values) const {
      return {
          -values[0] * static_cast<int32_t>(constants.limits[0]),
          -values[1] * static_cast<int32_t>(constants.limits[1]),
          -values[2] * static_cast<int32_t>(constants.limits[2]),
      };
   }

   constexpr size_t applyMultipliersAndOffset(const std::array<int32_t, 3>& values) const {
      return static_cast<size_t>(constants.offset + constants.multipliers[0] * values[0] +
                                 constants.multipliers[1] * values[1] + constants.multipliers[2] * values[2]);
   }

   constexpr size_t calculateIndex(std::array<int32_t, 3> cellIndex) const {
      const auto locality = localityMultipliers(cellIndex);                       // -1, 0, 1
      const auto ni = neighbourIndex(locality);                                   // 0-26
      const auto fallback = static_cast<int32_t>(constants.fallbackToCenter[ni]); // 0, 1
      const auto valid = fallback ^ 1;                                            // Opposite of fallback
      const auto addOffset = static_cast<int32_t>(constants.shift[ni]);           // 0, 1
      const auto offsets = shiftOffsets(locality);                                // -limits, 0, limits

      // If the given cell indices are valid, we'll check if the indices need to be shifted, and if so, we'll add
      // the offsets. Else we'll fall back to the center cell of the stencil.
      cellIndex[0] = valid * (cellIndex[0] + addOffset * offsets[0]) + fallback * centerCell[0];
      cellIndex[1] = valid * (cellIndex[1] + addOffset * offsets[1]) + fallback * centerCell[1];
      cellIndex[2] = valid * (cellIndex[2] + addOffset * offsets[2]) + fallback * centerCell[2];

      return applyMultipliersAndOffset(cellIndex);
   }
};
} // namespace fsgrid
