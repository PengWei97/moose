# Testing permeability from porosity
# Trivial test, checking calculated permeability is correct
# k = k_anisotropic * k0 * (1-phi0)^m/phi0^n * phi^n/(1-phi)^m
# Block 1 k0 twice that of block 0 so permeability is twice has high in block 1

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 3
    xmin = 0
    xmax = 3
  []
  [top_two_elements]
    type = SubdomainBoundingBoxGenerator
    input = gmg
    bottom_left = '1.1 0 0'
    top_right = '3.1 0 0'
    block_id = 1
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [pp]
    [InitialCondition]
      type = ConstantIC
      value = 0
    []
  []
[]

[Kernels]
  [flux]
    type = PorousFlowAdvectiveFlux
    gravity = '0 0 0'
    variable = pp
  []
[]

[BCs]
  [ptop]
    type = DirichletBC
    variable = pp
    boundary = right
    value = 0
  []
  [pbase]
    type = DirichletBC
    variable = pp
    boundary = left
    value = 1
  []
[]

[AuxVariables]
  [poro]
    order = CONSTANT
    family = MONOMIAL
  []
  [perm_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [perm_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [perm_z]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [poro]
    type = PorousFlowPropertyAux
    property = porosity
    variable = poro
  []
  [perm_x]
    type = PorousFlowPropertyAux
    property = permeability
    variable = perm_x
    row = 0
    column = 0
  []
  [perm_y]
    type = PorousFlowPropertyAux
    property = permeability
    variable = perm_y
    row = 1
    column = 1
  []
  [perm_z]
    type = PorousFlowPropertyAux
    property = permeability
    variable = perm_z
    row = 2
    column = 2
  []
[]

[Postprocessors]
  [perm_x_bottom]
    type = PointValue
    variable = perm_x
    point = '0 0 0'
  []
  [perm_y_bottom]
    type = PointValue
    variable = perm_y
    point = '0 0 0'
  []
  [perm_z_bottom]
    type = PointValue
    variable = perm_z
    point = '0 0 0'
  []
  [perm_x_top]
    type = PointValue
    variable = perm_x
    point = '3 0 0'
  []
  [perm_y_top]
    type = PointValue
    variable = perm_y
    point = '3 0 0'
  []
  [perm_z_top]
    type = PointValue
    variable = perm_z
    point = '3 0 0'
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    # unimportant in this fully-saturated test
    m = 0.8
    alpha = 1e-4
  []
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 2.2e9
    viscosity = 1e-3
    density0 = 1000
    thermal_expansion = 0
  []
[]

[AuxVariables]
  [A_var]
    order = CONSTANT
    family = MONOMIAL
  []
  [A_var_bad]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [A]
    type = ParsedAux
    variable = A_var
    expression = 'if(x<1.1,0.11552,0.23104)'
    use_xyzt = true
  []
  [A_bad]
    type = ParsedAux
    variable = A_var_bad
    expression = 'if(x<1.1,0.11552,-.01)'
    use_xyzt = true
  []
[]

[Materials]
  inactive = 'permeability_all permeability_0A permeability_1A var_error param_error'
  [permeability_0]
    type = PorousFlowPermeabilityKozenyCarman
    k_anisotropy = '1 0 0  0 2 0  0 0 0.1'
    poroperm_function = kozeny_carman_phi0
    k0 = 1e-10
    phi0 = 0.05
    m = 2
    n = 7
    block = 0
  []
  [permeability_1]
    type = PorousFlowPermeabilityKozenyCarman
    k_anisotropy = '1 0 0  0 2 0  0 0 0.1'
    poroperm_function = kozeny_carman_phi0
    k0 = 2e-10
    phi0 = 0.05
    m = 2
    n = 7
    block = 1
  []
  [permeability_0A]
    type = PorousFlowPermeabilityKozenyCarman
    k_anisotropy = '1 0 0  0 2 0  0 0 0.1'
    poroperm_function = kozeny_carman_A
    A = 0.11552
    m = 2
    n = 7
    block = 0
  []
  [permeability_1A]
    type = PorousFlowPermeabilityKozenyCarman
    k_anisotropy = '1 0 0  0 2 0  0 0 0.1'
    poroperm_function = kozeny_carman_A
    A = 0.23104
    m = 2
    n = 7
    block = 1
  []
  [permeability_all]
    type = PorousFlowPermeabilityKozenyCarmanFromVar
    k_anisotropy = '1 0 0  0 2 0  0 0 0.1'
    m = 2
    n = 7
    A = A_var
  []
  [var_error]
    type = PorousFlowPermeabilityKozenyCarmanFromVar
    k_anisotropy = '1 0 0  0 2 0  0 0 0.1'
    m = 2
    n = 7
    A = A_var_bad
  []
  [param_error]
    type = PorousFlowPermeabilityKozenyCarman
    k_anisotropy = '1 0 0  0 2 0  0 0 0.1'
    poroperm_function = kozeny_carman_A
    A = 0.23104
    phi0 = .01
    m = 2
    n = 7
  []

  [temperature]
    type = PorousFlowTemperature
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
  []
  [ppss]
    type = PorousFlow1PhaseP
    porepressure = pp
    capillary_pressure = pc
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  []
  [porosity]
    type = PorousFlowPorosity
    porosity_zero = 0.1
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 0 # unimportant in this fully-saturated situation
    phase = 0
  []
[]

[Preconditioning]
  [andy]
    type = SMP
    full = true
  []
[]

[Executioner]
  solve_type = Newton
  type = Steady
  l_tol = 1E-5
  nl_abs_tol = 1E-3
  nl_rel_tol = 1E-8
  l_max_its = 200
  nl_max_its = 400

  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              lu            gmres     200'
[]

[Outputs]
  csv = true
  execute_on = 'timestep_end'
[]
