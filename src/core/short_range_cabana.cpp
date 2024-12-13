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
#include <unordered_map>

// Dont know where to do this better
// ===================================================
using data_types = Cabana::MemberTypes<double[3], int, int, int, double[3]>;
using memory_space = Kokkos::SharedSpace;
using execution_space = Kokkos::DefaultExecutionSpace;

const int vector_length = 8;

// ===================================================

template <class BondKernel, class PairKernel,
          class VerletCriterion = detail::True>
void cabana_short_range(BondKernel bond_kernel, PairKernel pair_kernel,
                      CellStructure &cell_structure, double pair_cutoff,
                      double bond_cutoff, BoxGeometry &box_geo,
                      InteractionsNonBonded &nonbonded_ias,
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
    const int num_particles = particles.size() + ghost_particles.size();

    Cabana::AoSoA<data_types, memory_space, vector_length> particle_storage("particles", num_particles);
    auto slice_position = Cabana::slice<0>(particle_storage);
    auto slice_id = Cabana::slice<1>(particle_storage);
    auto slice_type = Cabana::slice<2>(particle_storage);
    auto slice_ghost = Cabana::slice<3>(particle_storage);
    auto slice_force = Cabana::slice<4>(particle_storage);

    std::unordered_map<int, int> id_to_index;

    int index = 0;
    // Write particles
    for (auto const& p : particles) {
      // do not use id here, we need to have a map from id to index
        auto const pos = p.pos();
        slice_position(index, 0) = pos[0];
        slice_position(index, 1) = pos[1];
        slice_position(index, 2) = pos[2];
        slice_id(index) = p.id();
        slice_type(index) = p.type();
        slice_ghost(index) = 0;

        id_to_index[p.id()] = index;
        index++;
    }

    // Write ghost particles
    for (auto const& p : ghost_particles) {
        // Check if the ID already exists in the map
        if (id_to_index.find(p.id()) == id_to_index.end()) {
            auto const pos = p.pos();
            slice_position(index, 0) = pos[0];
            slice_position(index, 1) = pos[1];
            slice_position(index, 2) = pos[2];
            slice_id(index) = p.id();
            slice_type(index) = p.type();
            slice_ghost(index) = 1;

            id_to_index[p.id()] = index;
            index++;
        }
    }

    std::cout << "Particles: " << particles.size() << " Ghosts: " << ghost_particles.size() << std::endl;
    std::cout << "index: " << index << std::endl;

    // create CustomVerletList
    using ListAlgorithm = Cabana::HalfNeighborTag;
    using ListType = Cabana::VerletList<memory_space, ListAlgorithm, Cabana::VerletLayout2D>;

    // gridmin etc
    double grid_min[3] = {0., 0., 0.};
    double grid_max[3] = {10, 10, 10};
    double cell_size_ratio = 1;
    double neighborhood_radius = 1;

    ListType verlet_list(slice_position, 0, slice_position.size(), neighborhood_radius, cell_size_ratio, grid_min, grid_max );

    // (Maybe have to move this to CellStructure.hpp)
    // build verlet list
    // If rebuild: use Algorithm::link_cell to add all pair to our custom verlet list

    auto kernel = [&](Particle &p1, Particle &p2) {
      std::cout << "Pair id: " << p1.id() << " " << p2.id() << " Index: " << id_to_index[p1.id()] << " " << id_to_index[p2.id()] << std::endl;
      //std::cout << "index" << id_to_index[p1.id()] << " " << id_to_index[p2.id()] << std::endl;
      verlet_list._data.addNeighbor(id_to_index[p1.id()], id_to_index[p2.id()]);
    };

    cell_structure.cabana_verlet_list_loop(kernel, verlet_criterion);

    auto first_neighbor_kernel = KOKKOS_LAMBDA(const int i, const int j) {
        
        std::cout << "Index: " << i << " " << j << "   Pair Id: " << slice_id(i) << " " << slice_id(j) << std::endl;

        
        if (slice_ghost(i) == 1 && slice_ghost(j) == 1) {
            return;
        }

        if (slice_id(i) == slice_id(j)) {
            return;
        }

        Utils::Vector3d pi = {slice_position(i, 0), slice_position(i, 1), slice_position(i, 2)};
        Utils::Vector3d pj = {slice_position(j, 0), slice_position(j, 1), slice_position(j, 2)};

        Utils::Vector3d dist_vec = box_geo.get_mi_vector(pj, pi);
        auto const dist = dist_vec.norm();

        IA_parameters const& ia_param = nonbonded_ias.get_ia_param(slice_type(i), slice_type(j));

        const ParticleForce kokkos_force = calc_central_radial_force(ia_param, dist_vec, dist);

        if (slice_ghost(i) == 0) {
            Kokkos::atomic_add(&slice_force(i, 0), kokkos_force.f[0]);
            Kokkos::atomic_add(&slice_force(i, 1), kokkos_force.f[1]);
            Kokkos::atomic_add(&slice_force(i, 2), kokkos_force.f[2]);
            ////std::cout << "Force i: " << slice_force(i,0) << " " << slice_force(i, 1) << " " << slice_force(i, 2) << " " << std::endl;

        }

        if (slice_ghost(j) == 0) {
            Kokkos::atomic_add(&slice_force(j, 0), -kokkos_force.f[0]);
            Kokkos::atomic_add(&slice_force(j, 1), -kokkos_force.f[1]);
            Kokkos::atomic_add(&slice_force(j, 2), -kokkos_force.f[2]);
            ////std::cout << "Force j: " << slice_force(j,0) << " " << slice_force(j, 1) << " " << slice_force(j, 2) << " " << std::endl;
        }
    };

    Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());

    Cabana::neighbor_parallel_for(policy, first_neighbor_kernel, verlet_list,
                                    Cabana::FirstNeighborsTag(),
                                    Cabana::SerialOpTag(), "verlet_list");
    
    Kokkos::fence(); 

    
    //std::size_t count Kokkos::atomic_fetch_add(&verlet_list._data.counts(3), 1);
    //verlet_list._data.neighbors(3, count) = 4;
    

    // print all neighbors in verlet_list._data
    for (int i = 0; i < verlet_list._data.counts.extent(0); i++) {
        for (int j = 0; j < verlet_list._data.counts(i); j++) {
            auto neighbor_index = verlet_list._data.neighbors(i, j);
            std::cout << "Particle " << i << " has neighbor " << neighbor_index << std::endl;
        }
    }

    for (auto & p : particles) {
        auto const id = id_to_index[p.id()];
        Utils::Vector3d f_vec{slice_force(id,0), slice_force(id, 1), slice_force(id, 2)};
        
        ParticleForce f(f_vec);
        p.force_and_torque() += f;
       
    }

    // loop over all pairs in the custom verlet list with the Cabana::NeighborList
    // use kernel

    // copy forces back to particles

    std::cout << "Iteration" << std::endl;

    // remove this
    //cell_structure.non_bonded_loop(pair_kernel, verlet_criterion);

    //cell_structure.cabana_link_cell([&](Particle &p1, Particle &p2) {
    //  std::cout << "Pair: " << p1.id() << " " << p2.id() << std::endl;
    //    }
    //);

  }
}