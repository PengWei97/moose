[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Executioner]
  type = Steady

  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
  perf_graph = true
[]

[MultiApps]
  [full_solve]
    type = FullSolveMultiApp
    # not setting app_type to use the same app type of parent, i.e. MooseTestApp
    execute_on = TIMESTEP_BEGIN
    positions = '0 0 0'
    input_files = sub.i
  []
[]

[AuxVariables]
  [pid]
    order = CONSTANT
    family = MONOMIAL
    [AuxKernel]
      type = ProcessorIDAux
      execute_on = 'initial'
    []
  []
[]

[Transfers]
  [pid]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = 'full_solve'
    source_variable = pid
    variable = parent_pid
    search_value_conflicts = false
  []
[]
