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
#include <unordered_set>
#include <utility>



// ===================================================

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2; // Combine the two hash values
    }
};

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

    // ===================================================
    // Setup Cabana Variables
    // ===================================================
    // Dont know where to do this better
    using data_types = Cabana::MemberTypes<double[3], int, int, int, double[3]>;
    using memory_space = Kokkos::SharedSpace;
    using execution_space = Kokkos::DefaultExecutionSpace;

    using ListAlgorithm = Cabana::HalfNeighborTag;
    using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm, Cabana::VerletLayout2D>;

    const int vector_length = 8;

    // ===================================================
    // Count unique particles and create Index map
    // ===================================================

    auto t1 = std::chrono::high_resolution_clock::now();

    std::unordered_map<int, int> id_to_index;
    int index = 0;

    for (auto const& p : particles) {
      id_to_index[p.id()] = index;
      index++;
    }

    for (auto const& p : ghost_particles) {
      if (id_to_index.find(p.id()) == id_to_index.end()) {
        id_to_index[p.id()] = index;
        index++;
      }
    }

    int number_of_unique_particles = index;

    auto t2 = std::chrono::high_resolution_clock::now();

    // ===================================================
    // Create and fill particle storage
    // ===================================================

    Cabana::AoSoA<data_types, memory_space, vector_length> particle_storage("particles", number_of_unique_particles);
    auto slice_position = Cabana::slice<0>(particle_storage);
    auto slice_id = Cabana::slice<1>(particle_storage);
    auto slice_type = Cabana::slice<2>(particle_storage);
    auto slice_ghost = Cabana::slice<3>(particle_storage);
    auto slice_force = Cabana::slice<4>(particle_storage);

    for (auto const& p : particles) {
      auto const pos = p.pos();
      slice_position(id_to_index[p.id()], 0) = pos[0];
      slice_position(id_to_index[p.id()], 1) = pos[1];
      slice_position(id_to_index[p.id()], 2) = pos[2];
      slice_id(id_to_index[p.id()]) = p.id();
      slice_type(id_to_index[p.id()]) = p.type();
      slice_ghost(id_to_index[p.id()]) = 0;
      slice_force(id_to_index[p.id()], 0) = 0.0;
      slice_force(id_to_index[p.id()], 1) = 0.0;
      slice_force(id_to_index[p.id()], 2) = 0.0;
    }

    for (auto const& p : ghost_particles) {
      if (id_to_index.find(p.id()) == id_to_index.end()) {
        auto const pos = p.pos();
        slice_position(id_to_index[p.id()], 0) = pos[0];
        slice_position(id_to_index[p.id()], 1) = pos[1];
        slice_position(id_to_index[p.id()], 2) = pos[2];
        slice_id(id_to_index[p.id()]) = p.id();
        slice_type(id_to_index[p.id()]) = p.type();
        slice_ghost(id_to_index[p.id()]) = 1;
        slice_force(id_to_index[p.id()], 0) = 0.0;
        slice_force(id_to_index[p.id()], 1) = 0.0;
        slice_force(id_to_index[p.id()], 2) = 0.0;
      }
    }

    auto t3 = std::chrono::high_resolution_clock::now();

    // ===================================================
    // Determine Verlet Pairs
    // ===================================================

    // TODO: this is not optimal
    // Either use dynamic memory allocation in custom verlet list
    // or
    // use a variant of this function to determine the maximum number of neighbors
    // and then allocate
    // Using a set to count particle pairs was slow


    // Neighbor counts, which particles has how many neighbors
    int neighbor_counts[number_of_unique_particles] = {0};
    //std::set<std::pair<int, int>> neighbor_map;

    // Counting kernel
    auto kernel_count = [&](Particle &p1, Particle &p2) {
      neighbor_counts[id_to_index[p1.id()]] = neighbor_counts[id_to_index[p1.id()]] + 1;
    };

    // Loop over all pairs to count neighbors
    cell_structure.cabana_verlet_list_loop(kernel_count, verlet_criterion);

    // Find max amount of neighbors from count for allocation in CustomVerletList
    auto max_neighbors = std::max_element(neighbor_counts, neighbor_counts + number_of_unique_particles);

    ListType verlet_list(slice_position, 0, slice_position.size(), *max_neighbors); 

    auto kernel = [&](Particle &p1, Particle &p2) {
      verlet_list.addNeighbor(id_to_index[p1.id()], id_to_index[p2.id()]);
    };



    cell_structure.cabana_verlet_list_loop(kernel, verlet_criterion);

    auto t4 = std::chrono::high_resolution_clock::now();

    // calculate max neighbors
    //std::cout << "Max neighbors: " << *max_neighbors << std::endl;

    // build empty customverletlist
    

    //for (auto const& pair : neighbor_map) {
    //  verlet_list.addNeighbor(pair.first, pair.second);
    //}

    auto t5 = std::chrono::high_resolution_clock::now();

    // fill customverletlist with pairs
    auto first_neighbor_kernel = KOKKOS_LAMBDA(const int i, const int j) {

        Utils::Vector3d pi = {slice_position(i, 0), slice_position(i, 1), slice_position(i, 2)};
        Utils::Vector3d pj = {slice_position(j, 0), slice_position(j, 1), slice_position(j, 2)};

        Utils::Vector3d const dist_vec = box_geo.get_mi_vector(pi, pj);
        //std::cout << "Dist vec: " << dist_vec[0] << " " << dist_vec[1] << " " << dist_vec[2] << std::endl;
        auto const dist = dist_vec.norm();

        IA_parameters const& ia_param = nonbonded_ias.get_ia_param(slice_type(i), slice_type(j));

        const ParticleForce kokkos_force = calc_central_radial_force(ia_param, dist_vec, dist);

        if (slice_ghost(i) == 0) {
            Kokkos::atomic_add(&slice_force(i, 0), kokkos_force.f[0]);
            Kokkos::atomic_add(&slice_force(i, 1), kokkos_force.f[1]);
            Kokkos::atomic_add(&slice_force(i, 2), kokkos_force.f[2]);
            //std::cout << "Force i: " << slice_force(i,0) << " " << slice_force(i, 1) << " " << slice_force(i, 2) << " " << std::endl;

        }

        if (slice_ghost(j) == 0) {
            Kokkos::atomic_add(&slice_force(j, 0), -kokkos_force.f[0]);
            Kokkos::atomic_add(&slice_force(j, 1), -kokkos_force.f[1]);
            Kokkos::atomic_add(&slice_force(j, 2), -kokkos_force.f[2]);
            //std::cout << "Force j: " << slice_force(j,0) << " " << slice_force(j, 1) << " " << slice_force(j, 2) << " " << std::endl;
        }
    };

    // loop over pairs with customverletlist
    Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());

    Cabana::neighbor_parallel_for(policy, first_neighbor_kernel, verlet_list,
                                    Cabana::FirstNeighborsTag(),
                                    Cabana::SerialOpTag(), "verlet_list");
    
    Kokkos::fence(); 

    auto t6 = std::chrono::high_resolution_clock::now();

    // copy forces back to particles
    for (auto & p : particles) {
        auto const id = id_to_index[p.id()];
        Utils::Vector3d f_vec{slice_force(id,0), slice_force(id, 1), slice_force(id, 2)};
        
        ParticleForce f(f_vec);
        p.force_and_torque() += f;
       
    }

    auto t7 = std::chrono::high_resolution_clock::now();

    std::cout << "Time to create particle storage: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() << " ns" << std::endl;
    std::cout << "Time to fill particle storage: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count() << " ns" << std::endl;
    std::cout << "Time to get verlet pairs: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count() << " ns" << std::endl;
    std::cout << "Time to create custom verlet list: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t4).count() << " ns" << std::endl;
    std::cout << "Time to run kernel: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t6 - t5).count() << " ns" << std::endl;
    std::cout << "Time to add forces: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t7 - t6).count() << " ns" << std::endl;
    std::cout << "Total time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t7 - t1).count() << " ns" << std::endl;
  }
}