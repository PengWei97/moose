//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterWrapperSolutionTransferBase.h"
#include "MultiApp.h"
#include "FEProblemBase.h"
#include "DisplacedProblem.h"
#include "InterWrapperMesh.h"

InputParameters
InterWrapperSolutionTransferBase::validParams()
{
  InputParameters params = MultiAppTransfer::validParams();
  params.addRequiredParam<std::vector<AuxVariableName>>("variable",
                                                        "The auxiliary variables to transfer.");
  return params;
}

InterWrapperSolutionTransferBase::InterWrapperSolutionTransferBase(
    const InputParameters & parameters)
  : MultiAppTransfer(parameters), _var_names(getParam<std::vector<AuxVariableName>>("variable"))
{
  if (_directions.contains(Transfer::FROM_MULTIAPP))
    paramError("from_multiapp", "This transfer works only into multi-app.");
}

void
InterWrapperSolutionTransferBase::initialSetup()
{
  MultiAppTransfer::initialSetup();
  for (std::size_t var_index = 0; var_index < _var_names.size(); ++var_index)
  {
    // Check source variable on regular subchannel problem
    MooseVariableFieldBase & from_var = _subproblem.getVariable(
        0, _var_names[var_index], Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_ANY);
    System & from_sys = from_var.sys().system();
    const auto & fe_type = from_sys.variable_type(from_var.number());

    if (fe_type.family != LAGRANGE || fe_type.order != FIRST)
      paramError("variable",
                 "This transfer requires a first order Lagrange variable for the source variable");

    // Check target variable in visualization mesh
    MooseVariableFieldBase & to_var = _to_problems[0]->getVariable(
        0, _var_names[var_index], Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_ANY);

    System & to_sys = to_var.sys().system();
    const auto & fe_type_target = to_sys.variable_type(to_var.number());

    if (fe_type_target.family != LAGRANGE || fe_type_target.order != FIRST)
      paramError("variable",
                 "This transfer requires a first order Lagrange variable for the source variable");
  }
}

void
InterWrapperSolutionTransferBase::execute()
{
  getAppInfo();

  switch (_current_direction)
  {
    case TO_MULTIAPP:
      transferToMultiApps();
      break;

    default:
      break;
  }
}

void
InterWrapperSolutionTransferBase::transferToMultiApps()
{
  mooseAssert(_from_meshes.size() == 1, "Only one source mesh can be active in this transfer.");
  if (dynamic_cast<InterWrapperMesh *>(_from_meshes[0]) == nullptr)
    mooseError("This transfer works only with InterWrapperMesh classes.");

  for (unsigned int i = 0; i < _multi_app->numGlobalApps(); i++)
    if (_multi_app->hasLocalApp(i))
      transferVarsToApp(i);
}

void
InterWrapperSolutionTransferBase::transferVarsToApp(unsigned int app_idx)
{
  transferNodalVars(app_idx);
}

void
InterWrapperSolutionTransferBase::transferNodalVars(unsigned int app_idx)
{
  Moose::ScopedCommSwapper swapper(_multi_app->comm());

  FEProblemBase & to_problem = _multi_app->appProblemBase(app_idx);
  MooseMesh * mesh = NULL;
  if (_displaced_target_mesh && to_problem.getDisplacedProblem())
    mesh = &to_problem.getDisplacedProblem()->mesh();
  else
    mesh = &to_problem.mesh();

  const InterWrapperMesh & from_mesh = dynamic_cast<InterWrapperMesh &>(*_from_meshes[0]);
  FEProblemBase & from_problem = *_from_problems[0];

  for (auto & node : mesh->getMesh().local_node_ptr_range())
  {
    Node * from_node = getFromNode(from_mesh, *node);

    for (auto & var_name : _var_names)
    {
      System * to_sys = find_sys(to_problem.es(), var_name);
      unsigned int to_sys_num = to_sys->number();
      unsigned int to_var_num = to_sys->variable_number(var_name);

      if (node->n_dofs(to_sys_num, to_var_num) > 0)
      {
        System * from_sys = find_sys(from_problem.es(), var_name);
        unsigned int from_sys_num = from_sys->number();
        unsigned int from_var_num = from_sys->variable_number(var_name);

        swapper.forceSwap();
        NumericVector<Real> * from_solution = from_sys->solution.get();
        dof_id_type from_dof = from_node->dof_number(from_sys_num, from_var_num, 0);
        Real from_value = (*from_solution)(from_dof);
        swapper.forceSwap();

        NumericVector<Real> & to_solution = _multi_app->appTransferVector(app_idx, var_name);
        dof_id_type to_dof = node->dof_number(to_sys_num, to_var_num, 0);
        to_solution.set(to_dof, from_value);
      }
    }
  }

  for (auto & var_name : _var_names)
  {
    _multi_app->appTransferVector(app_idx, var_name).close();
    find_sys(to_problem.es(), var_name)->update();
  }
}
