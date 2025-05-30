//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IntegratedBC.h"

/**
 * Boundary condition for convective heat flux where temperature and heat transfer coefficient are
 * given by material properties.
 */
class ConvectiveHeatFluxBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  ConvectiveHeatFluxBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Far-field temperature variable
  const MaterialProperty<Real> & _T_infinity;

  /// Convective heat transfer coefficient
  const MaterialProperty<Real> & _htc;

  /// Derivative of convective heat transfer coefficient with respect to temperature
  const MaterialProperty<Real> & _htc_dT;
};
