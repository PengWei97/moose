//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"

// libmesh
#include "libmesh/elem_range.h"
#include "libmesh/threads.h"

class FEProblemBase;

class ComputeInitialConditionThread
{
public:
  ComputeInitialConditionThread(FEProblemBase & fe_problem);
  // Splitting Constructor
  ComputeInitialConditionThread(ComputeInitialConditionThread & x, Threads::split split);

  void operator()(const libMesh::ConstElemRange & range);
  void join(const ComputeInitialConditionThread & /*y*/);

protected:
  /// Print information about the loop, mostly order of execution of objects
  void printGeneralExecutionInformation() const;

  FEProblemBase & _fe_problem;
  THREAD_ID _tid;
};
