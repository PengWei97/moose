//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ReactorApp.h"

class ReactorUnitApp : public ReactorApp
{
public:
  ReactorUnitApp(const InputParameters & parameters);
  virtual ~ReactorUnitApp();

  static InputParameters validParams();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);

  /// So that we can construct MeshGenerators with metadata in unit apps
  virtual bool constructingMeshGenerators() const override { return true; }
};
