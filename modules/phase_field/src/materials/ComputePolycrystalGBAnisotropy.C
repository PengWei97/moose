//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputePolycrystalGBAnisotropy.h"


registerMooseObject("PhaseFieldApp", ComputePolycrystalGBAnisotropy);

InputParameters
ComputePolycrystalGBAnisotropy::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Compute of grain boundary energy and grain boundary mobility based on misorientation");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides Grain ID according to element ID");  
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
    params.addParam<Real>(
      "GBenergy",
      1,
      "Grain boundary energy in J/m^2");
  params.addParam<Real>(
      "GBMobility",
      1,
      "GB mobility input in m^4/(J*s)");
  return params;   
}

ComputePolycrystalGBAnisotropy::ComputePolycrystalGBAnisotropy(const InputParameters & parameters)
  : Material(parameters),
    _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _GBEnergy(getParam<Real>("GBenergy")),
    _GBMobility(getParam<Real>("GBMobility")),
    _sigma_GB(declareProperty<Real>("sigma_GB")),
    _M_GB(declareProperty<Real>("M_GB")),
    _delta_theta(declareProperty<Real>("delta_theta"))
{
}

void
ComputePolycrystalGBAnisotropy::computeQpProperties()
{
  // Get list of active order parameters from grain tracker, std::vector<unsigned int>
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  unsigned int num_grain_id = 0; // the number of vaild graid IDs
  std::vector<std::vector<unsigned int>> grainID_varibaleIndex; // Create a vector of grain IDs to order parameter indices
  grainID_varibaleIndex.clear();

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];

    if (grain_id == FeatureFloodCount::invalid_id)
      continue;    
    
    grainID_varibaleIndex[num_grain_id][0] = grain_id;
    grainID_varibaleIndex[num_grain_id][1] = op_index; 
    num_grain_id++;
  }

  Real sum_val = 0;
  Real sum_delta_theta = 0;
  Real Val = 0.0;
  
  if (grainID_varibaleIndex.size() == 1) // inside the grain 
  {
    _sigma_GB[_qp] = 0; // grain boundary energy
    _M_GB[_qp] =0; // grain boundary mobility
    _delta_theta[_qp] =0;
  }
  else // on the grain boudary
  {
    for (unsigned int i = 0; i < grainID_varibaleIndex.size(); ++i)
    {
      for (unsigned int j = i+1; j < grainID_varibaleIndex.size(); ++j)
      {
        const RealVectorValue angles_i = _euler.getEulerAngles(grainID_varibaleIndex[i][0]); // Get the sequence of Euler angles for grain i
        const RealVectorValue angles_j = _euler.getEulerAngles(grainID_varibaleIndex[j][0]);

        unsigned int m = grainID_varibaleIndex[i][1]; // Get the order parameter index of grain i
        unsigned int n = grainID_varibaleIndex[j][1];

        Val = (100000.0 * ((*_vals[m])[_qp]) * ((*_vals[m])[_qp]) + 0.01) *
              (100000.0 * ((*_vals[n])[_qp]) * ((*_vals[n])[_qp]) + 0.01); // eta_i^2 * eta_j^2

        sum_val += Val;
        sum_delta_theta += std::abs(angles_i(0)-angles_j(0))*Val; // misorientation for only considering phi_1
      }
    }

    _delta_theta[_qp] = sum_delta_theta / sum_val;

    // Grain boundary energy according to the Read-Shockley law;
    const Real sigma_HGB = 1; // the grain boudary energy of a high angle grain boundary
    const Real delta_theta_HGB = 15; // the misorientation corresponding to the transition angle between low and high angle grain boudaries
    _sigma_GB[_qp] = sigma_HGB * _delta_theta[_qp] / delta_theta_HGB *std::log(1 - std::log(_delta_theta[_qp] / delta_theta_HGB));

    // Grain boundary mobility according to the sigmoidal law 
    const Real M_HGB = 1; // the grain boudary mobility of a high angle grain boundary
    const Real B = 5;
    const Real n = 4;
    _M_GB[_qp] = M_HGB * (1 - std::exp(-B * std::pow( _delta_theta[_qp] / delta_theta_HGB, n))); 
  } 
}

