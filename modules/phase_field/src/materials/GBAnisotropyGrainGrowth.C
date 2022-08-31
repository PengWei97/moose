//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBAnisotropyGrainGrowth.h"
#include "EulerAngleProvider.h"
#include "MooseMesh.h"

registerMooseObject("PhaseFieldApp", GBAnisotropyGrainGrowth);

#include <fstream>

InputParameters
GBAnisotropyGrainGrowth::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides Grain ID according to element ID");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");  
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addRequiredParam<Real>("wGB", "Diffuse GB width in nm");
  params.addCoupledVar("T", 450.0, "Temperature in Kelvin");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<Real>("time_scale", 1.0e-9, "Time scale in s, where default is ns");
  params.addParam<Real>("GBsigma_HAB",0.708, "the gb energy of a high angle grain boudary");
  params.addParam<Real>("GBmob_HAB", 2.5e-6, "the gb mobility of a high angle grain boudary");
  params.addParam<Real>("GBQ_HAB", 0.23, "the gb activate energy of a high angle grain boudary");
  params.addParam<Real>("rate1_HABvsLAB_mob", 0.5, "the ratio of low-angle gb to high-angle gb");
  params.addParam<Real>("rate2_HABvsLAB_mob", 0.5, "the initial ratio of low-angle gb to high-angle gb");
  params.addParam<Real>("rate1_HABvsLAB_sigma", 0.5, "the ratio of low-angle gb to high-angle gb");
  params.addParam<Real>("rate2_HABvsLAB_sigma", 0.5, "the ratio of low-angle gb to high-angle gb");
  params.addParam<Real>(
      "delta_sigma", 0.1, "factor determining inclination dependence of GB energy");
  params.addParam<Real>(
      "delta_mob", 0.1, "factor determining inclination dependence of GB mobility");
  params.addRequiredParam<bool>("inclination_anisotropy",
                                "The GB anisotropy inclination would be considered if true");
  params.addRequiredParam<bool>("gbEnergy_anisotropy",
                                "The GB energy anisotropy would be considered if true");
  params.addRequiredParam<bool>("gbMobility_anisotropy",
                                "The GB mobility anisotropy would be considered if true");                              
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

GBAnisotropyGrainGrowth::GBAnisotropyGrainGrowth(const InputParameters & parameters)
  : Material(parameters),
    _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider")),
    _wGB(getParam<Real>("wGB")),
    _mesh_dimension(_mesh.dimension()),
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _GBsigma_HAB(getParam<Real>("GBsigma_HAB")),
    _GBmob_HAB(getParam<Real>("GBmob_HAB")),
    _GBQ_HAB(getParam<Real>("GBQ_HAB")),
    _rate1_HABvsLAB_mob(getParam<Real>("rate1_HABvsLAB_mob")),
    _rate2_HABvsLAB_mob(getParam<Real>("rate2_HABvsLAB_mob")),
    _rate1_HABvsLAB_sigma(getParam<Real>("rate1_HABvsLAB_sigma")),
    _rate2_HABvsLAB_sigma(getParam<Real>("rate2_HABvsLAB_sigma")),
    _delta_sigma(getParam<Real>("delta_sigma")),
    _delta_mob(getParam<Real>("delta_mob")),
    _inclination_anisotropy(getParam<bool>("inclination_anisotropy")),
    _gbEnergy_anisotropy(getParam<bool>("gbEnergy_anisotropy")),
    _gbMobility_anisotropy(getParam<bool>("gbMobility_anisotropy")),
    _T(coupledValue("T")), //??
    _kappa(declareProperty<Real>("kappa_op")),
    _gamma(declareProperty<Real>("gamma_asymm")),
    _L(declareProperty<Real>("L")),
    _mu(declareProperty<Real>("mu")),
    _delta_theta(declareProperty<Real>("delta_theta")),
    _kb(8.617343e-5),      // Boltzmann constant in eV/K
    _JtoeV(6.24150974e18), // Joule to eV conversion
    _mu_qp(0.0),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _grad_vals(coupledGradients("v"))
{
  // reshape vectors
  _sigma.resize(_op_num);
  _mob.resize(_op_num);
  _Q.resize(_op_num);
  _kappa_gamma.resize(_op_num);
  _a_g2.resize(_op_num);
  // _grain_id.resize(_op_num);

  for (unsigned int op = 0; op < _op_num; ++op)
  {
    _sigma[op].resize(_op_num);
    _mob[op].resize(_op_num);
    _Q[op].resize(_op_num);
    _kappa_gamma[op].resize(_op_num);
    _a_g2[op].resize(_op_num);
  }
}

void
GBAnisotropyGrainGrowth::computeQpProperties()
{
  // calculate grain boundary energy, grain boundary mobility and activate energy based on the misorientation
  computeGBParamaterByMisorientaion();
   
  computeGBParamater();
  
  Real sum_kappa = 0.0;
  Real sum_gamma = 0.0;
  Real sum_L = 0.0;
  Real Val = 0.0;
  Real sum_val = 0.0;
  Real f_sigma = 1.0;
  Real f_mob = 1.0;
  Real gamma_value = 0.0;

  for (unsigned int m = 0; m < _op_num - 1; ++m)
  {
    for (unsigned int n = m + 1; n < _op_num; ++n) // m<n
    {
      gamma_value = _kappa_gamma[n][m];

      if (_inclination_anisotropy)
      {
        if (_mesh_dimension == 3)
          mooseError("This material doesn't support inclination dependence for 3D for now!");

        Real phi_ave = libMesh::pi * n / (2.0 * _op_num);
        Real sin_phi = std::sin(2.0 * phi_ave);
        Real cos_phi = std::cos(2.0 * phi_ave);

        Real a = (*_grad_vals[m])[_qp](0) - (*_grad_vals[n])[_qp](0);
        Real b = (*_grad_vals[m])[_qp](1) - (*_grad_vals[n])[_qp](1);
        Real ab = a * a + b * b + 1.0e-7; // for the sake of numerical convergence, the smaller the
                                          // more accurate, but more difficult to converge

        Real cos_2phi = cos_phi * (a * a - b * b) / ab + sin_phi * 2.0 * a * b / ab;
        Real cos_4phi = 2.0 * cos_2phi * cos_2phi - 1.0;

        f_sigma = 1.0 + _delta_sigma * cos_4phi;
        f_mob = 1.0 + _delta_mob * cos_4phi;

        Real g2 = _a_g2[n][m] * f_sigma;
        Real y = -5.288 * g2 * g2 * g2 * g2 - 0.09364 * g2 * g2 * g2 + 9.965 * g2 * g2 -
                 8.183 * g2 + 2.007;
        gamma_value = 1.0 / y;
      }

      Val = (100000.0 * ((*_vals[m])[_qp]) * ((*_vals[m])[_qp]) + 0.01) *
            (100000.0 * ((*_vals[n])[_qp]) * ((*_vals[n])[_qp]) + 0.01);

      sum_val += Val;
      sum_kappa += _kappa_gamma[m][n] * f_sigma * Val;
      sum_gamma += gamma_value * Val;
      // Following comes from substituting Eq. (36c) from the paper into (36b)
      sum_L += Val * _mob[m][n] * std::exp(-_Q[m][n] / (_kb * _T[_qp])) * f_mob * _mu_qp *
               _a_g2[n][m] / _sigma[m][n];
      // sum_L += Val * _mob[m][n] * f_mob * _mu_qp * _a_g2[n][m] / _sigma[m][n];
    }
  }

  _kappa[_qp] = sum_kappa / sum_val;
  _gamma[_qp] = sum_gamma / sum_val;
  _L[_qp] = sum_L / sum_val;
  _mu[_qp] = _mu_qp;
}

void
GBAnisotropyGrainGrowth::computeGBParamaterByMisorientaion()
{
  // Get list of active order parameters from grain tracker, std::vector<unsigned int>
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());  
  std::vector<unsigned int> grainID; // Create a vector of grain IDs
  std::vector<unsigned int> variableIndex; // Create a vector of order parameter indices
  Real bnd_val = 0;

  _delta_theta[_qp] = 0;
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index]; // grain id

    if (grain_id == FeatureFloodCount::invalid_id) // max 4294967295
      continue;    

    grainID.push_back(grain_id); // save the grain id
    variableIndex.push_back(op_index); // the id of the order paramater gr[0-num_op]
  }

  // initial the tensor of paramater
  Real data_sigma;
  Real data_mob;
  Real data_Q;

  Real GBsigma_LOW = _GBsigma_HAB * _rate2_HABvsLAB_sigma;
  Real GBmob_LOW = _GBmob_HAB * _rate2_HABvsLAB_mob;
  Real GBQ_LOW = _GBQ_HAB;

  for (unsigned int i = 0; i < _op_num; ++i) // column
  {
    std::vector<Real> row_sigma; // create an empty row of double values
    std::vector<Real> row_mob; 
    std::vector<Real> row_Q;
    for (unsigned int j = 0; j < _op_num; ++j) // line
    {
      data_sigma = GBsigma_LOW;
      data_mob = GBmob_LOW;   
      data_Q = GBQ_LOW;
      row_sigma.push_back(data_sigma);
      row_mob.push_back(data_mob);
      row_Q.push_back(data_Q);
    } // initialize three type paramaters
    
      _sigma[i] = row_sigma; // unit: J/m^2 GB energy

      _mob[i] = row_mob; // unit: m^4/(J*s) GB mobility

      _Q[i] = row_Q; // unit: eV Q
  }

  const Real & B = 5;
  const Real & n = 4;
  const Real & delta_theta_HAB = 15.0; // the misorientation corresponding to the transition angle between low and high angel grain boudaries.

  Real delta_euler = 0.0;
  // edit the tensor based on the delta_theta
  if (grainID.size() == 1 || grainID.size() == 0) // inside the grain   || grainID.size() == 0
  {
    _delta_theta[_qp] = 0.0;
  }
  else if (grainID.size() > 1)
  {
    for (unsigned int i = 0; i < grainID.size()-1; i++)
      for(unsigned int j = i+1; j < grainID.size(); j++)
      {
        const RealVectorValue angles_i = _euler.getEulerAngles(grainID[i]);
        const RealVectorValue angles_j = _euler.getEulerAngles(grainID[j]); // EulerAngles
        delta_euler = std::abs(angles_i(1)-angles_j(1))/35.0*15.0; // get the misorientation based grain id

        if (delta_euler > 0.0)
        {
          if (_gbMobility_anisotropy)
          {
            _mob[variableIndex[i]][variableIndex[j]] = _GBmob_HAB * ((1- std::exp(-B * std::pow( delta_euler / delta_theta_HAB, n))) * _rate1_HABvsLAB_mob + _rate2_HABvsLAB_mob );

            _mob[variableIndex[j]][variableIndex[i]] = _mob[variableIndex[i]][variableIndex[j]];
          }
          if (_gbEnergy_anisotropy)
          {
            if (delta_euler < delta_theta_HAB)
            {
              _sigma[variableIndex[i]][variableIndex[j]] = _GBsigma_HAB * ((delta_euler / delta_theta_HAB * (1 - std::log(delta_euler / delta_theta_HAB))) * _rate1_HABvsLAB_sigma + _rate2_HABvsLAB_sigma );
            }
            else
            {
              _sigma[variableIndex[i]][variableIndex[j]] = _GBsigma_HAB * ((15.0 / delta_theta_HAB * (1 - std::log(15.0 / delta_theta_HAB))) * _rate1_HABvsLAB_sigma + _rate2_HABvsLAB_sigma);
            }

            _sigma[variableIndex[j]][variableIndex[i]] = _sigma[variableIndex[i]][variableIndex[j]];  
          }
        }
      }
      _delta_theta[_qp] = 1.0;
  }
  
 
  // if (_delta_theta[_qp] > 80.0)
  // {
  //   for(unsigned int i = 0; i < variableIndex.size() - 1; i++)
  //     for(unsigned int j = i+1; j < variableIndex.size(); j++)
  //     {
  //       if (_gbMobility_anisotropy)
  //       {
  //         _mob[variableIndex[i]][variableIndex[j]] = _GBmob_HAB * ((1- std::exp(-B * std::pow(_delta_theta[_qp] / delta_theta_HAB, n))) * _rate1_HABvsLAB + _rate2_HABvsLAB );  // the sigmoidal law suggested by Humphreys
  //         _mob[variableIndex[j]][variableIndex[i]] = _mob[variableIndex[i]][variableIndex[j]];
  //       }

  //       if (_gbEnergy_anisotropy)
  //       {
  //         if (_delta_theta[_qp] < delta_theta_HAB)
  //         {
  //           _sigma[variableIndex[i]][variableIndex[j]] = _GBsigma_HAB * ((_delta_theta[_qp] / delta_theta_HAB * (1 - std::log(_delta_theta[_qp] / delta_theta_HAB))) * _rate1_HABvsLAB + _rate2_HABvsLAB ); 
  //         }
  //         else 
  //         {
  //           _sigma[variableIndex[i]][variableIndex[j]] = _GBsigma_HAB * ((15.0 / delta_theta_HAB * (1 - std::log(15.0 / delta_theta_HAB))) * _rate1_HABvsLAB + _rate2_HABvsLAB );
  //         }
            
  //         _sigma[variableIndex[j]][variableIndex[i]] = _sigma[variableIndex[i]][variableIndex[j]];
  //       }
  //     }
  // }
}

void
GBAnisotropyGrainGrowth::computeGBParamater()
{
  Real sigma_init;
  Real g2 = 0.0;
  Real f_interf = 0.0;
  Real a_0 = 0.75;
  Real a_star = 0.0;
  Real kappa_star = 0.0;
  Real gamma_star = 0.0;
  Real y = 0.0; // 1/gamma
  Real yyy = 0.0;

  Real sigma_big = 0.0;
  Real sigma_small = 0.0;

  for (unsigned int m = 0; m < _op_num - 1; ++m)
    for (unsigned int n = m + 1; n < _op_num; ++n)
    {
      // Convert units of mobility and energy
      _sigma[m][n] *= _JtoeV * (_length_scale * _length_scale); // eV/nm^2


      _mob[m][n] *= _time_scale / (_JtoeV * (_length_scale * _length_scale * _length_scale *
                                             _length_scale)); // Convert to nm^4/(eV*ns);

      if (m == 0 && n == 1)
      {
        sigma_big = _sigma[m][n];
        sigma_small = sigma_big;
      }

      else if (_sigma[m][n] > sigma_big)
        sigma_big = _sigma[m][n];

      else if (_sigma[m][n] < sigma_small)
        sigma_small = _sigma[m][n];
    }

  sigma_init = (sigma_big + sigma_small) / 2.0;
  _mu_qp = 6.0 * sigma_init / _wGB; // eV/Xm^2/Xm-- eV/Xm^3

  for (unsigned int m = 0; m < _op_num - 1; ++m)
    for (unsigned int n = m + 1; n < _op_num; ++n) // m<n
    {

      a_star = a_0;
      a_0 = 0.0;

      while (std::abs(a_0 - a_star) > 1.0e-9)
      {
        a_0 = a_star;
        kappa_star = a_0 * _wGB * _sigma[m][n];
        g2 = _sigma[m][n] * _sigma[m][n] / (kappa_star * _mu_qp);
        y = -5.288 * g2 * g2 * g2 * g2 - 0.09364 * g2 * g2 * g2 + 9.965 * g2 * g2 - 8.183 * g2 +
            2.007;
        gamma_star = 1 / y;
        yyy = y * y * y;
        f_interf = 0.05676 * yyy * yyy - 0.2924 * yyy * y * y + 0.6367 * yyy * y - 0.7749 * yyy +
                   0.6107 * y * y - 0.4324 * y + 0.2792;
        a_star = std::sqrt(f_interf / g2);
      }

      _kappa_gamma[m][n] = kappa_star; // upper triangle stores the discrete set of kappa values
      _kappa_gamma[n][m] = gamma_star; // lower triangle stores the discrete set of gamma values

      _a_g2[m][n] = a_star; // upper triangle stores "a" data.
      _a_g2[n][m] = g2;     // lower triangle stores "g2" data.
    }
}