/*
 * Copyright (C) 2019-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "BoxGeometry.hpp"

#include <utils/Vector.hpp>

#include <stdexcept>

BOOST_AUTO_TEST_CASE(periodicity_test) {
  /* getter/default */
  {
    auto const box = BoxGeometry{};

    BOOST_CHECK(box.periodic(0u));
    BOOST_CHECK(box.periodic(1u));
    BOOST_CHECK(box.periodic(2u));
  }

  /* setter */
  {
    auto box = BoxGeometry{};

    box.set_periodic(0u, false);
    BOOST_CHECK(not box.periodic(0u));
    box.set_periodic(1u, false);
    BOOST_CHECK(not box.periodic(1u));
    box.set_periodic(2u, false);
    BOOST_CHECK(not box.periodic(2u));
    BOOST_CHECK_THROW(box.set_periodic(3u, false), std::out_of_range);
  }
}

BOOST_AUTO_TEST_CASE(length_test) {
  /* getter/default */
  {
    auto const box = BoxGeometry{};

    BOOST_CHECK(Utils::Vector3d::broadcast(1.) == box.length());
    BOOST_CHECK(Utils::Vector3d::broadcast(1.) == box.length_inv());
    BOOST_CHECK(Utils::Vector3d::broadcast(0.5) == box.length_half());
  }

  /* setter */
  {
    auto box = BoxGeometry{};
    auto const box_l = Utils::Vector3d{1., 2., 3.};
    auto const box_l_inv =
        Utils::Vector3d{1. / box_l[0], 1. / box_l[1], 1. / box_l[2]};

    box.set_length(box_l);

    BOOST_CHECK(box.length() == box_l);
    BOOST_CHECK(box.length_inv() == box_l_inv);
    BOOST_CHECK(box.length_half() == 0.5 * box_l);
  }
}
