#ifdef MFEM_ENABLED

#include "EquationSystemProblemOperator.h"

namespace Moose::MFEM
{
void
EquationSystemProblemOperator::SetGridFunctions()
{
  _trial_var_names = GetEquationSystem()->TrialVarNames();
  ProblemOperator::SetGridFunctions();
}

void
EquationSystemProblemOperator::Init(mfem::BlockVector & X)
{
  ProblemOperator::Init(X);

  GetEquationSystem()->BuildEquationSystem();
}

void
EquationSystemProblemOperator::Solve(mfem::Vector &)
{
  GetEquationSystem()->BuildJacobian(_true_x, _true_rhs);

  if ((_problem.jacobian_solver->isLOR() || _problem.jacobian_preconditioner->isLOR()) &&
      _equation_system->_test_var_names.size() > 1)
    mooseError("LOR solve is only supported for single-variable systems");

  _problem.jacobian_solver->updateSolver(
      *_equation_system->_blfs.Get(_equation_system->_test_var_names.at(0)),
      _equation_system->_ess_tdof_lists.at(0));

  _problem.nonlinear_solver->SetSolver(*_problem.jacobian_solver->getSolver());
  _problem.nonlinear_solver->SetOperator(*GetEquationSystem());
  _problem.nonlinear_solver->Mult(_true_rhs, _true_x);

  GetEquationSystem()->RecoverFEMSolution(_true_x, _problem.gridfunctions);
}

} // namespace Moose::MFEM

#endif
