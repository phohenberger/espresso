/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#include <iostream>

#ifdef SHARED_MEMORY_PARALLELISM

#include <Cabana_Core.hpp>
#include "custom_verlet_list.hpp"
#include <unordered_map>

using data_types = Cabana::MemberTypes<double[3], double[3], int, int>;
using memory_space = Kokkos::SharedSpace;
using execution_space = Kokkos::DefaultExecutionSpace;

int const vector_length = 8;

using ListAlgorithm = Cabana::HalfNeighborTag;
using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm, Cabana::VerletLayout2D>;
using AoSoA = Cabana::AoSoA<data_types, memory_space, vector_length>;

template <class SliceDouble3, class SliceInt>
inline void write_particle(Particle const &p, int id, SliceDouble3 &s_position, SliceDouble3 &s_force, SliceInt &s_id, SliceInt &s_type) {
  auto const pos = p.pos();
  s_position(id, 0) = pos[0];
  s_position(id, 1) = pos[1];
  s_position(id, 2) = pos[2];
  s_id(id) = p.id();
  s_type(id) = p.type();
  s_force(id, 0) = 0.0;
  s_force(id, 1) = 0.0;
  s_force(id, 2) = 0.0;
}

class CabanaData {

    std::unordered_map<int, int> id_to_index;
    AoSoA particle_storage;

public:
    CabanaData() {};
    
    void save_local_particles_to_aosoa(ParticleRange local_particles) {
        int n_local_particles = local_particles.size();

        particle_storage = AoSoA("particles", n_local_particles);
        auto slice_position = Cabana::slice<0>(particle_storage);
        auto slice_force = Cabana::slice<1>(particle_storage);
        auto slice_id = Cabana::slice<2>(particle_storage);
        auto slice_type = Cabana::slice<3>(particle_storage);

        int index = 0;

        for (auto const& p : local_particles) {
            id_to_index[p.id()] = index;
            
            write_particle(p, index, slice_position, slice_force, slice_id, slice_type);
            index++;
        }

    }

    void read_local_particles_from_aosoa() {
        auto slice_position = Cabana::slice<0>(particle_storage);
        auto slice_force = Cabana::slice<1>(particle_storage);
        auto slice_id = Cabana::slice<2>(particle_storage);
        auto slice_type = Cabana::slice<3>(particle_storage);

        for (int i = 0; i < 1; ++i) {
            std::cout << "Particle " << i << ": "
                      << "Position: (" << slice_position(i, 0) << ", "
                      << slice_position(i, 1) << ", "
                      << slice_position(i, 2) << "), "
                      << "Force: (" << slice_force(i, 0) << ", "
                      << slice_force(i, 1) << ", "
                      << slice_force(i, 2) << "), "
                      << "ID: " << slice_id(i) << ", "
                      << "Type: " << slice_type(i) << std::endl;
        }
    }

    AoSoA& get_aosoa() {
        return particle_storage;
    }

    ~CabanaData() {};
};
#endif 

