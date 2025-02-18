#
# Copyright (C) 2010-2022 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Compare the diffusion coefficient of a single thermalized particle obtained
from the particle's mean square displacement and the auto correlation
function of its velocity to the expected value. Uses the
Observables/Correlators framework.
"""

import espressomd
import espressomd.accumulators
import espressomd.observables
import scipy.integrate
import numpy as np

gamma = 2.4
kT = 1.37
dt = 0.05

system = espressomd.System(box_l=[1.0, 1.0, 1.0])
np.random.seed(seed=42)

p = system.part.add(pos=(0, 0, 0))
system.time_step = dt
system.thermostat.set_langevin(kT=kT, gamma=gamma, seed=42)
system.cell_system.skin = 0.4
system.integrator.run(1000)

pos_obs = espressomd.observables.ParticlePositions(ids=(p.id,))
vel_obs = espressomd.observables.ParticleVelocities(ids=(p.id,))

c_pos = espressomd.accumulators.Correlator(
    obs1=pos_obs, tau_lin=16, tau_max=100., delta_N=10, compress1="discard1",
    corr_operation="square_distance_componentwise")
c_vel = espressomd.accumulators.Correlator(
    obs1=vel_obs, tau_lin=16, tau_max=20., delta_N=1, compress1="discard1",
    corr_operation="scalar_product")
system.auto_update_accumulators.add(c_pos)
system.auto_update_accumulators.add(c_vel)

system.integrator.run(1000000)

c_pos.finalize()
c_vel.finalize()

msd = np.column_stack((c_pos.lag_times(),
                       c_pos.sample_sizes(),
                       c_pos.result().reshape([-1, 3])))
vacf = np.column_stack((c_vel.lag_times(),
                        c_vel.sample_sizes(),
                        c_vel.result().reshape([-1, 1])))
np.savetxt("msd.dat", msd)
np.savetxt("vacf.dat", vacf)

# Integral of vacf via Green-Kubo
# D= 1/3 int_0^infty <v(t_0)v(t_0+t)> dt

# Integrate with trapezoidal rule
I = scipy.integrate.trapezoid(vacf[:, 2], vacf[:, 0])
ratio = 1. / 3. * I / (kT / gamma)
print("Ratio of measured and expected diffusion coefficients from Green-Kubo:",
      ratio)


def expected_msd(x):
    return 2. * kT / gamma * x


# Check MSD
print("Ratio of expected and measured msd")
print("#time ratio_x ratio_y ratio_z")
for i in range(4, msd.shape[0], 4):
    print(msd[i, 0], msd[i, 2:5] / expected_msd(msd[i, 0]))
