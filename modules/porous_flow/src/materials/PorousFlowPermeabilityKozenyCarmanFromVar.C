//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPermeabilityKozenyCarmanFromVar.h"

registerMooseObject("PorousFlowApp", PorousFlowPermeabilityKozenyCarmanFromVar);
registerMooseObject("PorousFlowApp", ADPorousFlowPermeabilityKozenyCarmanFromVar);

template <bool is_ad>
InputParameters
PorousFlowPermeabilityKozenyCarmanFromVarTempl<is_ad>::validParams()
{
  InputParameters params = PorousFlowPermeabilityKozenyCarmanBase::validParams();
  params.addRequiredCoupledVar("A", "Variable used in permeability function.");
  params.addClassDescription("This Material calculates the permeability tensor from the "
                             "Kozeny-Carman equation for spatially varying initial properties.");
  return params;
}

template <bool is_ad>
PorousFlowPermeabilityKozenyCarmanFromVarTempl<
    is_ad>::PorousFlowPermeabilityKozenyCarmanFromVarTempl(const InputParameters & parameters)
  : PorousFlowPermeabilityKozenyCarmanBaseTempl<is_ad>(parameters), _A(coupledValue("A"))
{
}

template <bool is_ad>
Real
PorousFlowPermeabilityKozenyCarmanFromVarTempl<is_ad>::computeA() const
{
  if (_A[_qp] < 0)
    mooseError("The variable A must be greater than zero; A = ", _A[_qp], ".");
  return _A[_qp];
}

template class PorousFlowPermeabilityKozenyCarmanFromVarTempl<false>;
template class PorousFlowPermeabilityKozenyCarmanFromVarTempl<true>;
