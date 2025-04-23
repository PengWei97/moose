#ifdef MFEM_ENABLED

#pragma once
#include "MFEMKernel.h"

/*
 * \f[
 * (\sigma u, \nabla V')
 * \f]
 */
class MFEMVectorFEWeakDivergenceKernel : public MFEMKernel<mfem::BilinearFormIntegrator>
{
public:
  static InputParameters validParams();

  MFEMVectorFEWeakDivergenceKernel(const InputParameters & parameters);

  virtual mfem::BilinearFormIntegrator * createIntegrator() override;

protected:
  const std::string _coef_name;
  mfem::Coefficient & _coef;
};

#endif
