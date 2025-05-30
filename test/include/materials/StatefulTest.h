//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 * Empty material for use in simple applications that don't need material properties.
 */
class StatefulTest : public Material
{
public:
  static InputParameters validParams();

  StatefulTest(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  // optional coupled variable
  const VariableValue * const _coupled_val;

  std::vector<std::string> _prop_names;
  std::vector<Real> _prop_values;

  unsigned int _num_props;

  std::vector<MaterialProperty<Real> *> _properties;
  std::vector<const MaterialProperty<Real> *> _properties_old;
  std::vector<const MaterialProperty<Real> *> _properties_older;
};
