/*
 * Copyright (C) 2010-2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "config/config.hpp"

#include "cell_system/CellStructure.hpp"

#ifdef CABANA

#include <Cabana_Core.hpp>

#endif

#include <cassert>

namespace detail {
/**
 * @brief Functor that returns true for
 *        any arguments.
 */

} // namespace detail

template <class BondKernel, class PairKernel,
          class VerletCriterion = detail::True>
void cabana_short_range(BondKernel bond_kernel, PairKernel pair_kernel,
                      CellStructure &cell_structure, double pair_cutoff,
                      double bond_cutoff,
                      ParticleRange particles, ParticleRange ghost_particles,
                      VerletCriterion const &verlet_criterion = {}) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  if (bond_cutoff >= 0.) {
    cell_structure.bond_loop(bond_kernel);
  }

  // Cabana short range loop
  if (pair_cutoff > 0.) {
    
    // save all particle positions in AoSoA

    // create CustomVerletList
    
    

    // (Maybe have to move this to CellStructure.hpp)
    // build verlet list
    // If rebuild: use Algorithm::link_cell to add all pair to our custom verlet list

    auto kernel = [&](Particle &p1, Particle &p2) {
      std::cout << "Pair: " << p1.id() << " " << p2.id() << std::endl;
    };

    cell_structure.cabana_verlet_list_loop(kernel, verlet_criterion);


    // loop over all pairs in the custom verlet list with the Cabana::NeighborList
    // use kernel

    // copy forces back to particles

    std::cout << "Iteration" << std::endl;

    // remove this
    cell_structure.non_bonded_loop(pair_kernel, verlet_criterion);

    //cell_structure.cabana_link_cell([&](Particle &p1, Particle &p2) {
    //  std::cout << "Pair: " << p1.id() << " " << p2.id() << std::endl;
    //    }
    //);

  }
}