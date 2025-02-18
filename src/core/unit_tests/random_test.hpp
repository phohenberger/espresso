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

#pragma once

#include <boost/test/unit_test.hpp>

/* Helper functions to compute random numbers covariance in a single pass */

#include <utils/Vector.hpp>
#include <utils/quaternion.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>
#include <boost/variant.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <numeric>
#include <tuple>
#include <vector>

namespace Utils {
using VariantVectorXd =
    boost::variant<double, Vector2d, Vector3d, Vector4d, Quaternion<double>>;
} // namespace Utils

using Utils::VariantVectorXd;

namespace {

using Utils::Vector;

class visitor_size : public boost::static_visitor<std::size_t> {
public:
  template <std::size_t N>
  std::size_t operator()(Vector<double, N> const &v) const {
    return v.size();
  }
  std::size_t operator()(Utils::Quaternion<double> const &) const { return 4u; }
  std::size_t operator()(double) const { return 1u; }
};

class visitor_get : public boost::static_visitor<double> {
public:
  template <std::size_t N>
  double operator()(Vector<double, N> const &v, std::size_t i) const {
    return v[i];
  }
  double operator()(Utils::Quaternion<double> const &q, std::size_t i) const {
    return q[i];
  }
  double operator()(double v, std::size_t i) const {
    assert(i == 0u);
    return v;
  }
};

std::size_t get_size(VariantVectorXd const &vec) {
  return boost::apply_visitor(visitor_size(), vec);
}

double get_value(VariantVectorXd const &vec, std::size_t i) {
  return boost::apply_visitor(
      std::bind(visitor_get(), std::placeholders::_1, i), vec);
}

template <typename T> auto square_matrix(std::size_t N) {
  return std::vector<std::vector<T>>(N, std::vector<T>(N));
}

} // namespace

/** Draw a large sample of 3D vectors from PRNGs and compute statistics.
 *  Parameter @p noise_function is a generator that returns @f$ N @f$ vectors
 *  of size @f$ M_i @f$. The following statistics are evaluated: @f$ N @f$ means
 *  and @f$ N @f$ variances (samples are uncorrelated across axes, so pooling
 *  them is fine), and a covariance and a correlation matrix of size
 *  @f$ \sum M_i @f$.
 */
template <typename NoiseKernel>
std::tuple<std::vector<double>, std::vector<double>,
           std::vector<std::vector<double>>, std::vector<std::vector<double>>>
noise_statistics(NoiseKernel &&noise_function, std::size_t sample_size) {

  // get size of the arrays and size of the triangular correlation matrix
  auto const first_value = noise_function();
  auto const n_vectors = first_value.size();
  std::vector<std::size_t> dimensions{};
  std::ranges::transform(first_value, std::back_inserter(dimensions), get_size);
  auto const matrix_dim = std::accumulate(dimensions.begin(), dimensions.end(),
                                          std::size_t{0u}, std::plus<>());

  // set up boost accumulators
  namespace ba = boost::accumulators;
  namespace bt = boost::accumulators::tag;
  using stat_variance = ba::stats<bt::mean, bt::variance(ba::lazy)>;
  using stat_covariance = ba::stats<bt::covariance<double, bt::covariate1>>;
  using boost_variance = ba::accumulator_set<double, stat_variance>;
  using boost_covariance = ba::accumulator_set<double, stat_covariance>;
  std::vector<boost_variance> acc_variance(n_vectors);
  auto acc_covariance = ::square_matrix<boost_covariance>(matrix_dim);

  // accumulate
  for (std::size_t step = 0u; step < sample_size; ++step) {
    auto const noise_tuple = noise_function();
    // for each vector, pool the random numbers of all columns
    for (std::size_t vec1 = 0u; vec1 < dimensions.size(); ++vec1) {
      for (std::size_t col1 = 0u; col1 < dimensions[vec1]; ++col1) {
        acc_variance[vec1](::get_value(noise_tuple[vec1], col1));
      }
    }
    // fill the covariance matrix (upper triangle)
    std::size_t index1 = 0u;
    for (std::size_t vec1 = 0u; vec1 < dimensions.size(); ++vec1) {
      for (std::size_t col1 = 0u; col1 < dimensions[vec1]; ++col1) {
        std::size_t index2 = index1;
        for (std::size_t vec2 = vec1; vec2 < dimensions.size(); ++vec2) {
          for (std::size_t col2 = (vec2 == vec1) ? col1 : std::size_t{0u};
               col2 < dimensions[vec2]; ++col2) {
            acc_covariance[index1][index2](
                ::get_value(noise_tuple[vec1], col1),
                ba::covariate1 = ::get_value(noise_tuple[vec2], col2));
            index2++;
          }
        }
        index1++;
      }
    }
  }

  // compute statistics
  std::vector<double> means(n_vectors);
  std::vector<double> variances(n_vectors);
  for (std::size_t i = 0u; i < n_vectors; ++i) {
    means[i] = ba::mean(acc_variance[i]);
    variances[i] = ba::variance(acc_variance[i]);
  }
  auto covariance = ::square_matrix<double>(matrix_dim);
  for (std::size_t i = 0u; i < matrix_dim; ++i) {
    for (std::size_t j = i; j < matrix_dim; ++j) {
      covariance[i][j] = covariance[j][i] =
          ba::covariance(acc_covariance[i][j]);
    }
  }
  auto correlation = ::square_matrix<double>(matrix_dim);
  for (std::size_t i = 0u; i < matrix_dim; ++i) {
    for (std::size_t j = i; j < matrix_dim; ++j) {
      correlation[i][j] = correlation[j][i] =
          covariance[i][j] / sqrt(covariance[i][i] * covariance[j][j]);
    }
  }

  return std::make_tuple(means, variances, covariance, correlation);
}

boost::test_tools::predicate_result correlation_almost_equal(
    std::vector<std::vector<double>> const &correlation_matrix, std::size_t i,
    std::size_t j, double reference, double threshold) {
  auto const value = correlation_matrix[i][j];
  auto const diff = std::abs(value - reference);
  if (diff > threshold) {
    boost::test_tools::predicate_result res(false);
    res.message() << "The correlation coefficient M[" << i << "][" << j << "]{"
                  << value << "} differs from " << reference << " by " << diff
                  << " (> " << threshold << ")";
    return res;
  }
  return true;
}
