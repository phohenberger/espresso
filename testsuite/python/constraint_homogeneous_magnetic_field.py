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

import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd


class HomogeneousMagneticFieldTest(ut.TestCase):

    system = espressomd.System(box_l=[3.0, 3.0, 3.0])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    np.random.seed(seed=42)

    def tearDown(self):
        self.system.constraints.clear()

    def test_setter_and_getter(self):
        H_field1 = np.array([0.0, 1.0, 0.0])
        H_field2 = np.array([3.533, 5.842, 0.127])

        H_constraint = espressomd.constraints.HomogeneousMagneticField(
            H=H_field1)

        np.testing.assert_almost_equal(np.copy(H_constraint.H), H_field1)

        H_constraint.H = H_field2
        np.testing.assert_almost_equal(np.copy(H_constraint.H), H_field2)

    def test_default_value(self):
        H_field_default = np.array([1.0, 0.0, 0.0])
        H_constraint = espressomd.constraints.HomogeneousMagneticField()
        np.testing.assert_almost_equal(
            np.copy(H_constraint.H),
            H_field_default)

    @utx.skipIfMissingFeatures(["DIPOLES"])
    def test_add_energy_and_forces(self):
        H_field = [5.0, 3.0, 2.0]
        dip_mom0 = [2.0, 6.0, 1.]
        dip_mom1 = [-1.0, 0.5, -0.2]

        # check that the dipolar energy is zero initially, ...
        self.assertEqual(self.system.analysis.energy()["dipolar"], 0.0)

        H_constraint = espressomd.constraints.HomogeneousMagneticField(
            H=H_field)
        self.system.constraints.add(H_constraint)

        # ... and also after adding the constraint
        self.assertEqual(self.system.analysis.energy()["dipolar"], 0.0)

        # check dipolar energy when adding dipole moments
        p0 = self.system.part.add(
            pos=[0, 0, 0], dip=dip_mom0, rotation=(True, True, True))
        self.assertAlmostEqual(self.system.analysis.energy()["dipolar"],
                               -np.dot(H_field, dip_mom0), delta=1e-7)
        p1 = self.system.part.add(
            pos=[1, 1, 1], dip=dip_mom1, rotation=(True, True, True))
        self.assertAlmostEqual(self.system.analysis.energy()["dipolar"],
                               -(np.dot(H_field, dip_mom0) +
                                 np.dot(H_field, dip_mom1)),
                               delta=1e-7)

        # check that running the integrator leads to expected torques
        self.system.integrator.run(0, recalc_forces=True)
        torque_expected0 = np.cross(dip_mom0, H_field)
        torque_expected1 = np.cross(dip_mom1, H_field)
        np.testing.assert_allclose(np.copy(p0.torque_lab), torque_expected0,
                                   atol=1e-12, rtol=1e-10)
        np.testing.assert_allclose(np.copy(p1.torque_lab), torque_expected1,
                                   atol=1e-12, rtol=1e-10)


if __name__ == "__main__":
    ut.main()
