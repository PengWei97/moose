#ifdef MFEM_ENABLED

#include "MFEMMixedScalarCurlKernel.h"
#include "MFEMProblem.h"

registerMooseObject("MooseApp", MFEMMixedScalarCurlKernel);

InputParameters
MFEMMixedScalarCurlKernel::validParams()
{
  InputParameters params = MFEMMixedBilinearFormKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the bilinear form "
                             "$(k\\vec\\nabla \\times \\vec u, v)_\\Omega$ "
                             "arising from the weak form of the scalar curl operator "
                             "$k\\vec\\nabla \\times u$. The vector must be 2D.");
  params.addParam<MFEMScalarCoefficientName>(
      "coefficient", "Name of scalar property k to multiply the integrator by.");
  return params;
}

MFEMMixedScalarCurlKernel::MFEMMixedScalarCurlKernel(const InputParameters & parameters)
  : MFEMMixedBilinearFormKernel(parameters),
    _coef_name(getParam<MFEMScalarCoefficientName>("coefficient")),
    _coef(getScalarCoefficient(_coef_name))
{
}

mfem::BilinearFormIntegrator *
MFEMMixedScalarCurlKernel::createBFIntegrator()
{
  return new mfem::MixedScalarCurlIntegrator(_coef);
}

#endif
