#ifdef MFEM_ENABLED

#pragma once
#include "MFEMKernel.h"

/*
 * \f[
 * (\beta u, u')
 * \f]
 */
class MFEMMassKernel : public MFEMKernel<mfem::BilinearFormIntegrator>
{
public:
  static InputParameters validParams();

  MFEMMassKernel(const InputParameters & parameters);

  virtual mfem::BilinearFormIntegrator * createIntegrator() override;

protected:
  const std::string _coef_name;
  mfem::Coefficient & _coef;
};

#endif
