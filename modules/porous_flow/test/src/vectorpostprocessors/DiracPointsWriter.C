//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiracPointsWriter.h"

#include "libmesh/parallel_algebra.h"

registerMooseObject("PorousFlowTestApp", DiracPointsWriter);

InputParameters
DiracPointsWriter::validParams()
{
  return GeneralVectorPostprocessor::validParams();
}

DiracPointsWriter::DiracPointsWriter(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters),
    _xs(declareVector("x")),
    _ys(declareVector("y")),
    _zs(declareVector("z"))
{
}

void
DiracPointsWriter::execute()
{
  _xs.clear();
  _ys.clear();
  _zs.clear();
  std::set<Point> points;
  for (auto & entry : _subproblem.diracKernelInfo().getPoints())
    if (entry.first->active())
      points.insert(entry.second.first.begin(), entry.second.first.end());

  // Not every processor might know about every point
  this->comm().set_union(points);

  for (auto & p : points)
  {
    _xs.push_back(p(0));
    _ys.push_back(p(1));
    _zs.push_back(p(2));
  }
}
