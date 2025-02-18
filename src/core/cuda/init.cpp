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

#include "config/config.hpp"

#ifdef CUDA

#include "init.hpp"
#include "utils.hpp"

#include "communication.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstring>
#include <iterator>
#include <set>
#include <vector>

/** Helper class for device sets.
 */
struct CompareDevices {
  bool operator()(const EspressoGpuDevice &a,
                  const EspressoGpuDevice &b) const {
    auto const name_comp = strncmp(a.proc_name, b.proc_name, 63);
    /* if both devices are from the same node, order by id */
    return (name_comp == 0) ? a.id < b.id : name_comp < 0;
  }
};

/** Gather list of CUDA devices on all nodes on the head node.
 *  It relies on <tt>MPI_Get_processor_name()</tt> to get a unique identifier
 *  of the physical node, as opposed to the logical rank of which there can
 *  be more than one per node.
 */
std::vector<EspressoGpuDevice> cuda_gather_gpus() {
  /* List of local devices */
  std::vector<EspressoGpuDevice> devices_local;
  /* Global unique device list (only relevant on the head node) */
  std::vector<EspressoGpuDevice> devices_global;

  int n_devices = 0;
  invoke_skip_cuda_exceptions(
      [&n_devices]() { n_devices = cuda_get_n_gpus(); });

  int proc_name_len;
  char proc_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(proc_name, &proc_name_len);
  proc_name[63] = '\0';

  invoke_skip_cuda_exceptions([&devices_local, n_devices, &proc_name]() {
    for (int i = 0; i < n_devices; ++i) {
      auto device = cuda_get_device_props(i);
      std::strncpy(device.proc_name, proc_name, 64);
      device.proc_name[63] = '\0';
      device.node = this_node;
      devices_local.emplace_back(device);
    }
  });

  auto const n_gpus = static_cast<int>(devices_local.size());
  auto const n_nodes = ::communicator.size;

  if (this_node == 0) {
    std::set<EspressoGpuDevice, CompareDevices> device_set;
    int *n_gpu_array = new int[static_cast<unsigned int>(n_nodes)];
    MPI_Gather(&n_gpus, 1, MPI_INT, n_gpu_array, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* insert local devices */
    std::ranges::copy(devices_local,
                      std::inserter(device_set, device_set.begin()));

    EspressoGpuDevice device;
    MPI_Status s;
    /* Get devices from other nodes */
    for (int i = 1; i < n_nodes; ++i) {
      for (int j = 0; j < n_gpu_array[i]; ++j) {
        MPI_Recv(&device, sizeof(EspressoGpuDevice), MPI_BYTE, i, 0,
                 MPI_COMM_WORLD, &s);
        device_set.insert(device);
      }
    }
    /* Copy unique devices to result, if any */
    std::ranges::copy(device_set, std::back_inserter(devices_global));
    delete[] n_gpu_array;
  } else {
    /* Send number of devices to head node */
    MPI_Gather(&n_gpus, 1, MPI_INT, nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /* Send devices to head node */
    for (auto const &device : devices_local) {
      MPI_Send(&device, sizeof(EspressoGpuDevice), MPI_BYTE, 0, 0,
               MPI_COMM_WORLD);
    }
  }
  return devices_global;
}

void cuda_on_program_start() {
  if (::communicator.this_node == 0) {
    invoke_skip_cuda_exceptions(cuda_init);
  }
}

#endif // CUDA
