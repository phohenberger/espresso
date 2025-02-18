/*
 * Copyright (C) 2014-2022 The ESPResSo project
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

#include "error_handling/RuntimeErrorCollector.hpp"

#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives.hpp>

#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace ErrorHandling {

RuntimeErrorCollector::RuntimeErrorCollector(boost::mpi::communicator comm)
    : m_comm(std::move(comm)) {}

RuntimeErrorCollector::~RuntimeErrorCollector() {
  if (!m_errors.empty()) {
    /* Print remaining error messages on destruction */
    std::cerr << "There were unhandled errors.\n";
    flush();
  }
}

void RuntimeErrorCollector::message(const RuntimeError &message) {
  m_errors.emplace_back(message);
}

void RuntimeErrorCollector::message(RuntimeError message) {
  m_errors.emplace_back(std::move(message));
}

void RuntimeErrorCollector::message(RuntimeError::ErrorLevel level,
                                    const std::string &msg,
                                    const char *function, const char *file,
                                    const int line) {
  m_errors.emplace_back(level, m_comm.rank(), msg, std::string(function),
                        std::string(file), line);
}

void RuntimeErrorCollector::warning(const std::string &msg,
                                    const char *function, const char *file,
                                    const int line) {
  m_errors.emplace_back(RuntimeError::ErrorLevel::WARNING, m_comm.rank(), msg,
                        std::string(function), std::string(file), line);
}

void RuntimeErrorCollector::warning(const char *msg, const char *function,
                                    const char *file, const int line) {
  warning(std::string(msg), function, file, line);
}

void RuntimeErrorCollector::warning(const std::ostringstream &mstr,
                                    const char *function, const char *file,
                                    const int line) {
  warning(mstr.str(), function, file, line);
}

void RuntimeErrorCollector::error(const std::string &msg, const char *function,
                                  const char *file, const int line) {
  m_errors.emplace_back(RuntimeError::ErrorLevel::ERROR, m_comm.rank(), msg,
                        std::string(function), std::string(file), line);
}

void RuntimeErrorCollector::error(const char *msg, const char *function,
                                  const char *file, const int line) {
  error(std::string(msg), function, file, line);
}

void RuntimeErrorCollector::error(const std::ostringstream &mstr,
                                  const char *function, const char *file,
                                  const int line) {
  error(mstr.str(), function, file, line);
}

int RuntimeErrorCollector::count() const {
  return boost::mpi::all_reduce(m_comm, static_cast<int>(m_errors.size()),
                                std::plus<>());
}

int RuntimeErrorCollector::count(RuntimeError::ErrorLevel level) {
  return static_cast<int>(std::ranges::count_if(
      m_errors, [level](auto const &e) { return e.level() >= level; }));
}

void RuntimeErrorCollector::clear() { m_errors.clear(); }

void RuntimeErrorCollector::flush() {
  for (auto const &e : m_errors) {
    std::cerr << e.format() << std::endl;
  }
  this->clear();
}

std::vector<RuntimeError> RuntimeErrorCollector::gather() {
  std::vector<RuntimeError> all_errors{};
  std::swap(all_errors, m_errors);

  Utils::Mpi::gather_buffer(all_errors, m_comm);

  return all_errors;
}

void RuntimeErrorCollector::gather_local() {
  Utils::Mpi::gather_buffer(m_errors, m_comm);

  this->clear();
}

} // namespace ErrorHandling
