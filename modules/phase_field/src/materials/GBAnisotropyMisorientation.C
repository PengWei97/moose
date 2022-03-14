//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBAnisotropyMisorientation.h"

registerMooseObject("PhaseFieldApp", GBAnisotropyMisorientation);

InputParameters
GBAnisotropyMisorientation::validParams()
{
  InputParameters params = Material::validParams();
  params.addCoupledVar("T", 300.0, "Temperature in Kelvin");
  params.addRequiredParam<Real>("wGB", "Diffuse GB width in nm");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<Real>("time_scale", 1.0e-9, "Time scale in s, where default is ns");
  params.addParam<Real>("molar_volume_value",
                        7.11e-6,
                        "molar volume of material in m^3/mol, by default it's the value of copper"); // 没用
  params.addRequiredParam<FileName>("Anisotropic_GB_file_name",
                                    "Name of the file containing: 1)GB mobility prefactor; 2) GB "
                                    "migration activation energy; 3)GB energy");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

GBAnisotropyMisorientation::GBAnisotropyMisorientation(const InputParameters & parameters)
  : Material(parameters), 
    _mesh_dimension(_mesh.dimension()),
    _wGB(getParam<Real>("wGB")),
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _M_V(getParam<Real>("molar_volume_value")),
    _Anisotropic_GB_file_name(getParam<FileName>("Anisotropic_GB_file_name")),
    _T(coupledValue("T")),
    _kappa(declareProperty<Real>("kappa_op")),
    _gamma(declareProperty<Real>("gamma_asymm")), // gamma
    _L(declareProperty<Real>("L")),
    _mu(declareProperty<Real>("mu")),
    _molar_volume(declareProperty<Real>("molar_volume")),
    _entropy_diff(declareProperty<Real>("entropy_diff")),
    _act_wGB(declareProperty<Real>("act_wGB")),
    _kb(8.617343e-5),      // Boltzmann constant in eV/K
    _JtoeV(6.24150974e18), // Joule to eV conversion
    _mu_qp(0.0),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _grad_vals(coupledGradients("v"))
{

  //  初始化晶界能、晶界迁移率、激活能、梯度自由能，梯度项-gamma，a_g^2
  // reshape vectors
  _sigma.resize(_op_num); // sigma, the gb energy, // unit: J/m^2
  _mob.resize(_op_num); // mu, the gb mobility, unit: m^4/(J*s)
  _Q.resize(_op_num); // Q, unit: eV
  _kappa_gamma.resize(_op_num); // kappa and gamma
  _a_g2.resize(_op_num); // ??

  for (unsigned int op = 0; op < _op_num; ++op)
  {
    _sigma[op].resize(_op_num);
    _mob[op].resize(_op_num);
    _Q[op].resize(_op_num);
    _kappa_gamma[op].resize(_op_num);
    _a_g2[op].resize(_op_num);
  }

  // 从文件中读取晶界能、晶界迁移率和激活能
  // Read in data from "Anisotropic_GB_file_name"
  std::ifstream inFile(_Anisotropic_GB_file_name.c_str());

  if (!inFile)
    paramError("Anisotropic_GB_file_name", "Can't open GB anisotropy input file");

  for (unsigned int i = 0; i < 2; ++i)
    inFile.ignore(255, '\n'); // ignore line

  Real data;
  for (unsigned int i = 0; i < 3 * _op_num; ++i) // 列
  {
    std::vector<Real> row; // create an empty row of double values
    for (unsigned int j = 0; j < _op_num; ++j) // 行
    {
      inFile >> data;
      row.push_back(data);
    }

    if (i < _op_num)
      _sigma[i] = row; // unit: J/m^2

    else if (i < 2 * _op_num)
      _mob[i - _op_num] = row; // unit: m^4/(J*s)

    else
      _Q[i - 2 * _op_num] = row; // unit: eV
  }

  inFile.close();


  Real sigma_init; //
  Real g2 = 0.0; 
  Real f_interf = 0.0; // f_{0, \text { interf }}(\gamma)}
  Real a_0 = 0.75; // a_init(gamma_init) or a_k
  Real a_star = 0.0; // a*(kappa*, gamma*)
  Real kappa_star = 0.0; 
  Real gamma_star = 0.0;
  Real y = 0.0; // 1/gamma
  Real yyy = 0.0; // 1/gamma^3

  // 用于设定sigma_init初始值
  Real sigma_big = 0.0;
  Real sigma_small = 0.0;

  // 对晶界能和晶界迁移率修改单位，并找到最大的晶界能和最小的晶界能
  for (unsigned int m = 0; m < _op_num - 1; ++m)
    for (unsigned int n = m + 1; n < _op_num; ++n) \\ 输出上三角的结果
    {
      // Convert units of mobility and energy
      _sigma[m][n] *= _JtoeV * (_length_scale * _length_scale); // eV/nm^2

      _mob[m][n] *= _time_scale / (_JtoeV * (_length_scale * _length_scale * _length_scale *
                                             _length_scale)); // Convert to nm^4/(eV*ns);

      if (m == 0 && n == 1) // 初始化sigma_big, sigma_small
      {
        sigma_big = _sigma[m][n];
        sigma_small = sigma_big;
      }

      else if (_sigma[m][n] > sigma_big)
        sigma_big = _sigma[m][n];

      else if (_sigma[m][n] < sigma_small)
        sigma_small = _sigma[m][n];
    }

  // 设置初始晶界能和局部自由能密度函数的前置因子
  sigma_init = (sigma_big + sigma_small) / 2.0;
  _mu_qp = 6.0 * sigma_init / _wGB; // 3/4 * 0.125 * sigma / l, model coefficient

  // 对于晶粒m-晶粒n的晶界
  for (unsigned int m = 0; m < _op_num - 1; ++m)
    for (unsigned int n = m + 1; n < _op_num; ++n) // 矩阵的上三角
    {
      a_star = a_0; // 0.75
      a_0 = 0.0; // a_0 = sqrt*(f_{0, interf}(gamma)) / g(gamma)

      while (std::abs(a_0 - a_star) > 1.0e-9) // whie (a_0 != a*)
      {
        a_0 = a_star;
        kappa_star = a_0 * _wGB * _sigma[m][n]; // (eq-36b)
        g2 = _sigma[m][n] * _sigma[m][n] / (kappa_star * _mu_qp); // (eq-12)
        y = -5.288 * g2 * g2 * g2 * g2 - 0.09364 * g2 * g2 * g2 + 9.965 * g2 * g2 - 8.183 * g2 +
            2.007; // g^{-1} ??
        gamma_star = 1 / y; // gamma* = g^{-1}
        yyy = y * y * y;
        f_interf = 0.05676 * yyy * yyy - 0.2924 * yyy * y * y + 0.6367 * yyy * y - 0.7749 * yyy +
                   0.6107 * y * y - 0.4324 * y + 0.2792; // f_interf(gamma*)
        a_star = std::sqrt(f_interf / g2);
      }

      _kappa_gamma[m][n] = kappa_star; // upper triangle stores the discrete set of kappa values
      _kappa_gamma[n][m] = gamma_star; // lower triangle stores the discrete set of gamma values

      _a_g2[m][n] = a_star; // upper triangle stores "a" data. a*
      _a_g2[n][m] = g2;     // lower triangle stores "g2" data. (eq-12)
    }
}

void
GBAnisotropyMisorientation::computeQpProperties()
{
  Real sum_kappa = 0.0;
  Real sum_gamma = 0.0;
  Real sum_L = 0.0;
  Real Val = 0.0;
  Real sum_val = 0.0;
  Real f_mob = 1.0;
  Real gamma_value = 0.0;

  for (unsigned int m = 0; m < _op_num - 1; ++m)
  {
    for (unsigned int n = m + 1; n < _op_num; ++n) // m<n 上三角
    {
      gamma_value = _kappa_gamma[n][m]; // gamma

      Val = (100000.0 * ((*_vals[m])[_qp]) * ((*_vals[m])[_qp]) + 0.01) *
            (100000.0 * ((*_vals[n])[_qp]) * ((*_vals[n])[_qp]) + 0.01);

      sum_val += Val;
      sum_kappa += _kappa_gamma[m][n] * Val; // kappa
      sum_gamma += gamma_value * Val;
      // Following comes from substituting Eq. (36c) from the paper into (36b)
      sum_L += Val * _mob[m][n] * std::exp(-_Q[m][n] / (_kb * _T[_qp])) * f_mob * _mu_qp *
               _a_g2[n][m] / _sigma[m][n]; // the sum of L
    }
  }

  _kappa[_qp] = sum_kappa / sum_val;
  _gamma[_qp] = sum_gamma / sum_val;
  _L[_qp] = sum_L / sum_val;
  _mu[_qp] = _mu_qp;

  _molar_volume[_qp] =
      _M_V / (_length_scale * _length_scale * _length_scale); // m^3/mol converted to ls^3/mol
  _entropy_diff[_qp] = 9.5 * _JtoeV;                          // J/(K mol) converted to eV(K mol)
  _act_wGB[_qp] = 0.5e-9 / _length_scale;                     // 0.5 nm
}
