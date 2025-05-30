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

/**
 * An action for testing InputParameters::applyParameters
 */
class ApplyCoupledVariablesTestAction : public Action
{
public:
  /**
   * Class constructor
   */
  static InputParameters validParams();

  ApplyCoupledVariablesTestAction(const InputParameters & params);

  /**
   * Class destructor
   */
  virtual ~ApplyCoupledVariablesTestAction();

  /**
   * Creates an action for adding a Kernel
   */
  virtual void act();
};
