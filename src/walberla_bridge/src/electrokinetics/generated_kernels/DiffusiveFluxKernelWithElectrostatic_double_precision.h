//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \\file DiffusiveFluxKernelWithElectrostatic_double_precision.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7, lbmpy v1.3.7, sympy v1.12.1,
// lbmpy_walberla/pystencils_walberla from waLBerla commit
// f36fa0a68bae59f0b516f6587ea8fa7c24a41141

#pragma once
#include "core/DataTypes.h"
#include "core/logging/Logging.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"

#include <functional>
#include <unordered_map>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreorder"
#endif

namespace walberla {
namespace pystencils {

class DiffusiveFluxKernelWithElectrostatic_double_precision {
public:
  DiffusiveFluxKernelWithElectrostatic_double_precision(
      BlockDataID jID_, BlockDataID phiID_, BlockDataID rhoID_, double D,
      double f_ext_0, double f_ext_1, double f_ext_2, double kT, double z)
      : jID(jID_), phiID(phiID_), rhoID(rhoID_), D_(D), f_ext_0_(f_ext_0),
        f_ext_1_(f_ext_1), f_ext_2_(f_ext_2), kT_(kT), z_(z) {}

  void run(IBlock *block);

  void runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers, IBlock *block);

  void operator()(IBlock *block) { run(block); }

  static std::function<void(IBlock *)> getSweep(
      const shared_ptr<DiffusiveFluxKernelWithElectrostatic_double_precision>
          &kernel) {
    return [kernel](IBlock *b) { kernel->run(b); };
  }

  static std::function<void(IBlock *)> getSweepOnCellInterval(
      const shared_ptr<DiffusiveFluxKernelWithElectrostatic_double_precision>
          &kernel,
      const shared_ptr<StructuredBlockStorage> &blocks,
      const CellInterval &globalCellInterval, cell_idx_t ghostLayers = 1) {
    return [kernel, blocks, globalCellInterval, ghostLayers](IBlock *b) {
      kernel->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b);
    };
  }

  std::function<void(IBlock *)> getSweep() {
    return [this](IBlock *b) { this->run(b); };
  }

  std::function<void(IBlock *)>
  getSweepOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks,
                         const CellInterval &globalCellInterval,
                         cell_idx_t ghostLayers = 1) {
    return [this, blocks, globalCellInterval, ghostLayers](IBlock *b) {
      this->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b);
    };
  }

  void configure(const shared_ptr<StructuredBlockStorage> & /*blocks*/,
                 IBlock * /*block*/) {}

  inline double getD() const { return D_; }
  inline double getF_ext_0() const { return f_ext_0_; }
  inline double getF_ext_1() const { return f_ext_1_; }
  inline double getF_ext_2() const { return f_ext_2_; }
  inline double getKt() const { return kT_; }
  inline double getZ() const { return z_; }
  inline void setD(const double value) { D_ = value; }
  inline void setF_ext_0(const double value) { f_ext_0_ = value; }
  inline void setF_ext_1(const double value) { f_ext_1_ = value; }
  inline void setF_ext_2(const double value) { f_ext_2_ = value; }
  inline void setKt(const double value) { kT_ = value; }
  inline void setZ(const double value) { z_ = value; }

private:
  BlockDataID jID;
  BlockDataID phiID;

public:
  inline void setPhiID(BlockDataID phiID_) { phiID = phiID_; }

private:
  BlockDataID rhoID;
  double D_;
  double f_ext_0_;
  double f_ext_1_;
  double f_ext_2_;
  double kT_;
  double z_;
};

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif
