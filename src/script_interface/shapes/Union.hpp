/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef SCRIPT_INTERFACE_SHAPES_SHAPE_UNION_HPP
#define SCRIPT_INTERFACE_SHAPES_SHAPE_UNION_HPP

#include "Shape.hpp"
#include "script_interface/ObjectList.hpp"

#include <shapes/Union.hpp>

#include <memory>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Shapes {

class Union : public ObjectList<Shape, Shape> {
public:
  using Base = ObjectList<Shape, Shape>;
  using value_type = typename Base::value_type;

  Union() : m_core_shape(std::make_shared<::Shapes::Union>()) {}

  ~Union() override { do_destruct(); }

private:
  bool has_in_core(value_type const &obj_ptr) const override {
    return m_core_shape->contains(obj_ptr->shape());
  }
  void add_in_core(value_type const &obj_ptr) override {
    m_core_shape->add(obj_ptr->shape());
  }
  void remove_in_core(value_type const &obj_ptr) final {
    m_core_shape->remove(obj_ptr->shape());
  }

public:
  std::shared_ptr<::Shapes::Shape> shape() const override {
    return m_core_shape;
  }

private:
  std::shared_ptr<::Shapes::Union> m_core_shape;
  std::vector<std::shared_ptr<Shapes::Shape>> m_shapes;
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
