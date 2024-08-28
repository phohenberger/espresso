#pragma once

#include "config/config.hpp"

#ifdef CABANA
#include <Cabana_Core.hpp>
#include <Kokkos_Core.hpp> // Necessary?
#include <cmath> // Include cmath for sqrt
#endif

#include "cell_system/CellStructure.hpp"

using data_types = Cabana::MemberTypes<double[3], int, int, int>;
using memory_space = Kokkos::SharedSpace;
using execution_space = Kokkos::DefaultExecutionSpace;

const int vector_length = 64;

void verlet_loop_cabana(Cabana::AoSoA<data_types, memory_space, vector_length> aosoa_particles, CellStructure& cell_structure, BoxGeometry& box_geo, double cutoff, InteractionsNonBonded& nonbonded_ias) {
  // build verlet list
  auto slice_position = Cabana::slice<0>(aosoa_particles);
  auto slice_id = Cabana::slice<1>(aosoa_particles);
  auto slice_type = Cabana::slice<2>(aosoa_particles);
  auto slice_ghost = Cabana::slice<3>(aosoa_particles);


  // use HalfNeighbor tag so each pair only get saved once, and then apply forces to both
  using ListAlgorithm = Cabana::HalfNeighborTag;
  using ListType = Cabana::VerletList<memory_space, ListAlgorithm, Cabana::VerletLayout2D>;

  // Cell ratio?
  // neighborhood_radius + cutoff?
  // 
  double grid_min[3] = {0.0, 0.0, 0.0};
  auto box_l = box_geo.length();
  double grid_max[3] = {box_l[0], box_l[1], box_l[2]};
  double neighborhood_radius = cutoff; // + cell_structure.get_verlet_skin(); ?
  double cell_ratio = 1.0;
  
  ListType verlet_list(slice_position, 0, slice_position.size(), neighborhood_radius,
                        cell_ratio, grid_min, grid_max);

  auto first_neighbor_kernel = KOKKOS_LAMBDA( const int i, const int j ) {
        if (slice_ghost(i) == 1 && slice_ghost(j) == 1) {
          return;
        }

        IA_parameters const& ia_param = nonbonded_ias.get_ia_param(slice_type(i), slice_type(j));           
        // Calc Distance
        // TODO: Use minimum image distance?
        double dx = slice_position(i, 0) - slice_position(j, 0);
        double dy = slice_position(i, 1) - slice_position(j, 1);
        double dz = slice_position(i, 2) - slice_position(j, 2);
        double r2 = dx*dx + dy*dy + dz*dz;
        double dist = std::sqrt(r2);
        Utils::Vector3d dist_vec = {dx, dy, dz};

        ParticleForce force = calc_central_radial_force(ia_param, dist_vec,dist);
        // Kokkos::atomic_add(p1.f, force)
        // Kokkos::atomic_add(p2.f, -force)

        // TODO save particle force and then update them all at once?
  };

  Kokkos::RangePolicy<execution_space> policy( 0, aosoa_particles.size() );

  Cabana::neighbor_parallel_for( policy, first_neighbor_kernel, verlet_list,
                                   Cabana::FirstNeighborsTag(),
                                   Cabana::SerialOpTag(), "ex_1st_serial" );
  Kokkos::fence();
}

void short_range_cabana(ParticleRange& particles, ParticleRange& ghost_particles, CellStructure& cell_structure, BoxGeometry& box_geo, InteractionsNonBonded& nonbonded_ias, double cutoff) {
  std::cout << "=====" << std::endl;
  std::cout << "calling short_range_cabana" << std::endl;

  // Build the AoSoA
  // pos: double[3], particle id: int, type: int (ghost = 1, real = 0)
  const int particle_num = particles.size(); // + ghost_particles.size();

  Cabana::AoSoA<data_types, memory_space, vector_length> aosoa_particles( "all_particles" , particle_num);

  auto slice_position = Cabana::slice<0>(aosoa_particles);
  auto slice_id = Cabana::slice<1>(aosoa_particles);
  auto slice_type = Cabana::slice<2>(aosoa_particles);
  auto slice_ghost = Cabana::slice<3>(aosoa_particles);

  // Fill the AoSoA
  int id = 0;
  for (auto const& p : particles) {
    slice_position(id, 0) = p.pos()[0];
    slice_position(id, 1) = p.pos()[1];
    slice_position(id, 2) = p.pos()[2];
    slice_id(id) = p.id();
    slice_type(id) = p.type();
    slice_ghost(id) = 0;
    id++;
  }

  /*for (auto const& ghost_p : ghost_particles) {

    auto unfolded_position = box_geo.unfolded_position(ghost_p.pos(), ghost_p.image_box());
    slice_position(id, 0) = unfolded_position[0];
    slice_position(id, 1) = unfolded_position[1];
    slice_position(id, 2) = unfolded_position[2];
    slice_id(id) = ghost_p.id();
    slice_type(id) = 1;
    id++;
  }*/

  // TODO: use particle types to get ia params
  // TODO maybe only give slices to the verlet loop and not whole aosoa

  verlet_loop_cabana(aosoa_particles, cell_structure, box_geo, cutoff, nonbonded_ias);
  //cell_structure.cabana_verlet_loop(aosoa_particles);
}

