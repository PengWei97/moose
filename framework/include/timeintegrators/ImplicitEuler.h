//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeIntegrator.h"

/**
 * Implicit Euler's method
 */
class ImplicitEuler : public TimeIntegrator
{
public:
  static InputParameters validParams();

  ImplicitEuler(const InputParameters & parameters);
  virtual ~ImplicitEuler();

  virtual int order() override { return 1; }
  virtual void computeTimeDerivatives() override;
  void computeADTimeDerivatives(ADReal & ad_u_dot,
                                const dof_id_type & dof,
                                ADReal & ad_u_dotdot) const override;
  virtual void postResidual(NumericVector<Number> & residual) override;
  virtual bool overridesSolve() const override { return false; }

  virtual Real timeDerivativeRHSContribution(const dof_id_type dof_id,
                                             const std::vector<Real> & factors) const override;
  virtual Real timeDerivativeMatrixContribution(const Real factor) const override;

protected:
  /**
   * Helper function that actually does the math for computing the time derivative
   */
  template <typename T, typename T2>
  void computeTimeDerivativeHelper(T & u_dot, const T2 & u_old) const;
};

template <typename T, typename T2>
void
ImplicitEuler::computeTimeDerivativeHelper(T & u_dot, const T2 & u_old) const
{
  u_dot -= u_old;
  u_dot *= 1. / _dt;
}
