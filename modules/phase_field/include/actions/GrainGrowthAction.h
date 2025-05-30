//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "Action.h"

#include "libmesh/fe_type.h"

class GrainGrowthAction : public Action
{
public:
  static InputParameters validParams();

  GrainGrowthAction(const InputParameters & params);

  virtual void act();

protected:
  void addVariables();
  void addBnds(const std::string & name_base);

  /// number of variables and variable name base for variable creation
  const unsigned int _op_num;
  const std::string _var_name_base;

  /// FEType for the variable being created
  const libMesh::FEType _fe_type;

  /// Take initial values from file?
  const bool _initial_from_file;

  /// use AD objects where possible
  const bool _use_ad;
};
