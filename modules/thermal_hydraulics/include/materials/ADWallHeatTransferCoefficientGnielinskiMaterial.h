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
 * Computes wall heat transfer coefficient for gases and water using the Gnielinski correlation
 */
class ADWallHeatTransferCoefficientGnielinskiMaterial : public Material
{
public:
  ADWallHeatTransferCoefficientGnielinskiMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Wall heat transfer coefficient
  ADMaterialProperty<Real> & _Hw;
  /// Density
  const ADMaterialProperty<Real> & _rho;
  /// Velocity
  const ADMaterialProperty<Real> & _vel;
  /// Hydraulic diameter
  const ADMaterialProperty<Real> & _D_h;
  /// Thermal conductivity
  const ADMaterialProperty<Real> & _k;
  /// Dynamic viscosity
  const ADMaterialProperty<Real> & _mu;
  /// Specific heat capacity
  const ADMaterialProperty<Real> & _cp;
  /// Fluid temperature
  const ADMaterialProperty<Real> & _T;
  /// Wall temperature
  const ADMaterialProperty<Real> & _T_wall;

public:
  static InputParameters validParams();
};
