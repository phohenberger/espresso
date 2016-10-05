#ifndef LBBOUNDARIES_LBBOUNDARY_HPP
#define LBBOUNDARIES_LBBOUNDARY_HPP

#include <memory>

#include "shapes/NoWhere.hpp"
#include "shapes/Shape.hpp"
#include "Vector.hpp"

namespace LBBoundaries {

class LBBoundary {
public:
  LBBoundary()
      : m_shape(std::make_shared<Shapes::NoWhere>()),
        m_velocity(Vector3d{0, 0, 0}),
        m_force(Vector3d{0, 0, 0}) {
#ifdef EK_BOUNDARIES
    m_charge_density = 0.0;
    m_net_charge = 0.0;
#endif
  }

  /* Calculate distance from the lbboundary */
  int calc_dist(const double *pos, double *dist, double *vec) const {
    return m_shape->calculate_dist(pos, dist, vec);
  }


  void set_shape(std::shared_ptr<Shapes::Shape> const &shape) {
    m_shape = shape;
  }

  void set_velocity(Vector3d velocity) { m_velocity = velocity; }

  Shapes::Shape const &shape() const { return *m_shape; }
  Vector3d &velocity() { return m_velocity; }
  Vector3d &force() { return m_force; }

#ifdef EK_BOUNDARIES //TODO: ugly. Better would be a class EKBoundaries, deriving from LBBoundaries, but that requires completely different initialization infrastructure.
  void set_charge_density(float charge_density) { m_charge_density = charge_density; }
  void set_net_charge(float net_charge) { m_net_charge = net_charge; }

  float charge_density() { return m_charge_density; }
  float net_charge() { return m_net_charge; }
#endif

private:
  /** Private methods */
  /* The actual boundary */
  std::shared_ptr<::LBBoundaries::LBBoundary> m_constraint;

  /** Private data members */
  std::shared_ptr<Shapes::Shape> m_shape; //TODO: I dont like this being a pointer just to get around the virtual limitations
  Vector3d m_velocity;
  Vector3d m_force;

#ifdef EK_BOUNDARIES //TODO: ugly. Better would be a class EKBoundaries, deriving from LBBoundaries, but that requires completely different initialization infrastructure.
  float m_charge_density;
  float m_net_charge;
#endif
};

} /* namespace LBBoundaries */

#endif
