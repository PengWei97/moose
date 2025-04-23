#ifdef MFEM_ENABLED

#pragma once
#include "libmesh/ignore_warnings.h"
#include "mfem/miniapps/common/pfem_extras.hpp"
#include "libmesh/restore_warnings.h"
#include "MFEMAuxKernel.h"

/*
Class to set an H(curl) auxvariable to be the gradient of a H1 scalar variable.
*/
class MFEMGradAux : public MFEMAuxKernel
{
public:
  static InputParameters validParams();

  MFEMGradAux(const InputParameters & parameters);

  virtual ~MFEMGradAux() = default;

  // Computes the auxvariable.
  virtual void execute() override;

protected:
  // Name of source MFEMVariable to take the gradient of.
  const VariableName _source_var_name;
  // Reference to source gridfunction.
  const mfem::ParGridFunction & _source_var;
  // Scalar factor to multiply the result by.
  const mfem::real_t _scale_factor;
  // Grad operator
  mfem::common::ParDiscreteGradOperator _grad;
};

#endif
