[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 2
  xmin = 0
  xmax = 10
  ymin = 0
  ymax = 2
[]

[Variables]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = 'x0:=5.0;thk:=0.5;m:=2;r:=abs(x-x0);v:=exp(-(r/thk)^m);0.1+0.1*v'
    [../]
  [../]
  [./mu]
  [../]
  [./jx]
  [../]
  [./jy]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[AuxVariables]
  [./gb]
    family = LAGRANGE
    order  = FIRST
  [../]
  [./creep_strain_xx]
    family = MONOMIAL
    order  = CONSTANT
  [../]
  [./creep_strain_yy]
    family = MONOMIAL
    order  = CONSTANT
  [../]
  [./creep_strain_xy]
    family = MONOMIAL
    order  = CONSTANT
  [../]
  [./stress_xx]
    family = MONOMIAL
    order  = CONSTANT
  [../]
  [./stress_yy]
    family = MONOMIAL
    order  = CONSTANT
  [../]
  [./stress_xy]
    family = MONOMIAL
    order  = CONSTANT
  [../]
[]

[Kernels]
  [./conc]
    type = CHSplitConcentration
    variable = c
    mobility = mobility_prop
    chemical_potential_var = mu
  [../]
  [./chempot]
    type = CHSplitChemicalPotential
    variable = mu
    chemical_potential_prop = mu_prop
    c = c
  [../]
  [./flux_x]
    type = CHSplitFlux
    variable = jx
    component = 0
    mobility_name = mobility_prop
    mu = mu
    c = c
  [../]
  [./flux_y]
    type = CHSplitFlux
    variable = jy
    component = 1
    mobility_name = mobility_prop
    mu = mu
    c = c
  [../]
  [./time]
    type = TimeDerivative
    variable = c
  [../]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
[]

[AuxKernels]
  [./gb]
    type = FunctionAux
    variable = gb
    function = 'x0:=5.0;thk:=0.5;m:=2;r:=abs(x-x0);v:=exp(-(r/thk)^m);v'
  [../]
  [./creep_strain_xx]
    type = RankTwoAux
    variable = creep_strain_xx
    rank_two_tensor = creep_strain
    index_i = 0
    index_j = 0
  [../]
  [./creep_strain_yy]
    type = RankTwoAux
    variable = creep_strain_yy
    rank_two_tensor = creep_strain
    index_i = 1
    index_j = 1
  [../]
  [./creep_strain_xy]
    type = RankTwoAux
    variable = creep_strain_xy
    rank_two_tensor = creep_strain
    index_i = 0
    index_j = 1
  [../]
  [./stress_xx]
    type = RankTwoAux
    variable = stress_xx
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
  [../]
  [./stress_xy]
    type = RankTwoAux
    variable = stress_xy
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
  [../]
[]

[Materials]
  [./chemical_potential]
    type = DerivativeParsedMaterial
    block = 0
    property_name = mu_prop
    coupled_variables = c
    expression = 'c'
    derivative_order = 1
  [../]
  [./var_dependence]
    type = DerivativeParsedMaterial
    block = 0
    expression = 'c*(1.0-c)'
    coupled_variables = c
    property_name = var_dep
    derivative_order = 1
  [../]
  [./mobility]
    type = CompositeMobilityTensor
    block = 0
    M_name = mobility_prop
    tensors = diffusivity
    weights = var_dep
    args = c
  [../]
  [./phase_normal]
    type = PhaseNormalTensor
    phase = gb
    normal_tensor_name = gb_normal
  [../]
  [./aniso_tensor]
    type = GBDependentAnisotropicTensor
    gb = gb
    bulk_parameter = 0.1
    gb_parameter = 1
    gb_normal_tensor_name = gb_normal
    gb_tensor_prop_name = aniso_tensor
  [../]
  [./diffusivity]
    type = GBDependentDiffusivity
    gb = gb
    bulk_parameter = 0.1
    gb_parameter = 1
    gb_normal_tensor_name = gb_normal
    gb_tensor_prop_name = diffusivity
  [../]
  [./diffuse_strain_increment]
    type = FluxBasedStrainIncrement
    xflux = jx
    yflux = jy
    gb = gb
    property_name = diffuse
  [../]
  [./diffuse_creep_strain]
    type = SumTensorIncrements
    tensor_name = creep_strain
    coupled_tensor_increment_names = diffuse
  [../]
  [./strain]
   type = ComputeIncrementalStrain
    displacements = 'disp_x disp_y'
  [../]
  [./stress]
    type = ComputeStrainIncrementBasedStress
    inelastic_strain_names = creep_strain
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '120.0 80.0'
    fill_method = symmetric_isotropic
  [../]
[]

[BCs]
  [./Periodic]
    [./cbc]
      auto_direction = 'x y'
      variable = c
    [../]
  [../]
  [./fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK

  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      31                  preonly       lu           1'

  nl_rel_tol = 1e-10
  nl_max_its = 5

  l_tol = 1e-4
  l_max_its = 20

  dt = 1
  num_steps = 5
[]

[Preconditioning]
  [./smp]
     type = SMP
     full = true
  [../]
[]

[Outputs]
  exodus = true
[]
