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
// test hash function for custom pair
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2; // Combine the two hash values
    }
};

// ===================================================

template <class BondKernel,
          class VerletCriterion = detail::True>
void cabana_short_range(BondKernel bond_kernel,
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

    auto t1 = std::chrono::high_resolution_clock::now();

    // ===================================================
    // Setup Cabana Variables
    // ===================================================
    // Dont know where to do this better
    using data_types = Cabana::MemberTypes<double[3], int, int, int, double[3]>;
    using memory_space = Kokkos::SharedSpace;
    using execution_space = Kokkos::DefaultExecutionSpace;

    using ListAlgorithm = Cabana::HalfNeighborTag;
    using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm, Cabana::VerletLayout2D>;

    const int vector_length = 4;

    // ===================================================
    // Count unique particles and create Index map
    // ===================================================

    // TODO:
    // Move this to cell structure and only recalculate when verlet recalculates

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

    // print all particles pairs in map
    for (auto const& p : id_to_index) {
      //std::cout << p.first << " " << p.second << std::endl;
    }

    const int number_of_unique_particles = index;

    auto t2 = std::chrono::high_resolution_clock::now();

    // ===================================================
    // Create and fill particle storage
    // ===================================================

    // Maybe move to cell structure and only update positions?
    // Only recalculate id/type etc when verlet recalculates
    // Probably not worth it, writing is actually very fast.


    Cabana::AoSoA<data_types, memory_space, vector_length> particle_storage("particles", number_of_unique_particles);
    auto slice_position = Cabana::slice<0>(particle_storage);
    auto slice_id = Cabana::slice<1>(particle_storage);
    auto slice_type = Cabana::slice<2>(particle_storage);
    auto slice_ghost = Cabana::slice<3>(particle_storage);
    auto slice_force = Cabana::slice<4>(particle_storage);

    for (auto const& p : particles) {
      auto const pos = p.pos();
      auto const id = id_to_index.at(p.id());
      slice_position(id, 0) = pos[0];
      slice_position(id, 1) = pos[1];
      slice_position(id, 2) = pos[2];
      slice_id(id) = p.id();
      slice_type(id) = p.type();
      slice_ghost(id) = 0;
      slice_force(id, 0) = 0.0;
      slice_force(id, 1) = 0.0;
      slice_force(id, 2) = 0.0;
    }

    for (auto const& p : ghost_particles) {
      if (id_to_index.find(p.id()) == id_to_index.end()) {
        auto const pos = p.pos();
        auto const id = id_to_index.at(p.id());
        slice_position(id, 0) = pos[0];
        slice_position(id, 1) = pos[1];
        slice_position(id, 2) = pos[2];
        slice_id(id) = p.id();
        slice_type(id) = p.type();
        slice_ghost(id) = 1;
        slice_force(id, 0) = 0.0;
        slice_force(id, 1) = 0.0;
        slice_force(id, 2) = 0.0;
      }
    }

    auto t3 = std::chrono::high_resolution_clock::now();

    // ===================================================
    // Get Verlet Pairs and Fill list
    // ===================================================

    auto t4 = std::chrono::high_resolution_clock::now();
    // TODO: maybe save this until next verlet list rebuild
    // but could cause problems with particle indices changing

    // Potential speed up:
    // Instead of only saving counts up ahead, find a fast way to save all neighbors and then parallel iterate over them.
    
    ListType verlet_list;

    if (true) {
      verlet_list = ListType(slice_position, 0, slice_position.size(), 16); 

      auto kernel = [&](Particle &p1, Particle &p2) {
        verlet_list.addNeighbor(id_to_index.at(p1.id()), id_to_index.at(p2.id()));
      };
        
      cell_structure.cabana_verlet_list_loop(kernel, verlet_criterion);
    } else {

      bool rebuild = cell_structure.get_rebuild_verlet_list();
      std::cout << "Rebuild: " << rebuild << std::endl;
      
      
      if (rebuild) {

        verlet_list = ListType(slice_position, 0, slice_position.size(), 16);
        
        auto kernel = [&](Particle &p1, Particle &p2) {
          verlet_list.addNeighbor(id_to_index.at(p1.id()), id_to_index.at(p2.id()));
        };

        cell_structure.cabana_verlet_list_loop(kernel, verlet_criterion);
        cell_structure.saveObject(verlet_list);

      } else {
        verlet_list = cell_structure.getObject<ListType>();
      }
    }
    // fill customverletlist with pairs
    auto first_neighbor_kernel = KOKKOS_LAMBDA(const int i, const int j) {

        Utils::Vector3d const pi = {slice_position(i, 0), slice_position(i, 1), slice_position(i, 2)};
        Utils::Vector3d const pj = {slice_position(j, 0), slice_position(j, 1), slice_position(j, 2)};

        Utils::Vector3d const dist_vec = box_geo.get_mi_vector(pi, pj);
        auto const dist = dist_vec.norm();

        IA_parameters const& ia_param = nonbonded_ias.get_ia_param(slice_type(i), slice_type(j));

        const ParticleForce kokkos_force = calc_central_radial_force(ia_param, dist_vec, dist);

        Kokkos::atomic_add(&slice_force(i, 0), kokkos_force.f[0]);
        Kokkos::atomic_add(&slice_force(i, 1), kokkos_force.f[1]);
        Kokkos::atomic_add(&slice_force(i, 2), kokkos_force.f[2]);
        
        Kokkos::atomic_add(&slice_force(j, 0), -kokkos_force.f[0]);
        Kokkos::atomic_add(&slice_force(j, 1), -kokkos_force.f[1]);
        Kokkos::atomic_add(&slice_force(j, 2), -kokkos_force.f[2]);  
    };

    auto t5 = std::chrono::high_resolution_clock::now();

    // ===================================================
    // Execute Kernel
    // ===================================================
    
    Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());


    // TODO: Add option to switch "SerialOpTag" Between "TeamOpTag"
    Cabana::neighbor_parallel_for(policy, first_neighbor_kernel, verlet_list,
                                    Cabana::FirstNeighborsTag(),
                                    Cabana::TeamOpTag(), "verlet_list");
    
    Kokkos::fence(); 


    auto t6 = std::chrono::high_resolution_clock::now();

    // ===================================================
    // Add forces to particles
    // ===================================================

    for (auto & p : particles) {
        auto const id = id_to_index.at(p.id());
        Utils::Vector3d f_vec{slice_force(id,0), slice_force(id, 1), slice_force(id, 2)};
        
        ParticleForce f(f_vec);
        p.force_and_torque() += f;
       
    }

    // Gives wrong results
    /*
    for (auto & p : ghost_particles) {
        // TODO: CHECK ID IN INDEX
        auto const id = id_to_index.at(p.id());
        Utils::Vector3d f_vec{slice_force(id,0), slice_force(id, 1), slice_force(id, 2)};
        
        ParticleForce f(f_vec);
        p.force_and_torque() += f;
    }*/

    auto t7 = std::chrono::high_resolution_clock::now();

    std::cout << "Map creation: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us" << std::endl;
    std::cout << "Time to create particle storage: " << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << " us" << std::endl;
    std::cout << "Time to count neighbors: " << std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() << " us" << std::endl;
    std::cout << "Time to fill verlet list: " << std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count() << " us" << std::endl;
    std::cout << "Time to execute kernel: " << std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count() << " us" << std::endl;
    std::cout << "Time to add forces: " << std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count() << " us" << std::endl;
    std::cout << "Total time: " << std::chrono::duration_cast<std::chrono::microseconds>(t7 - t1).count() << " us" << std::endl;
    std::cout << "Number of unique particles: " << number_of_unique_particles << std::endl;
  }
}