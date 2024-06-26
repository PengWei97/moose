##########################################################
# Simulation of Gallium Melting Experiment
# Ref: Gau, C., & Viskanta, R. (1986). Melting and solidification of a pure metal on a vertical wall.
# Key physics: melting/solidification, convective heat transfer, natural convection
##########################################################

mu = 1.81e-3
rho_solid = 6093
rho_liquid = 6093
k_solid = 32
k_liquid = 32
cp_solid = 381.5
cp_liquid = 381.5
L = 80160
alpha_b = 1.2e-4
T_solidus = 302.93
T_liquidus = '${fparse T_solidus + 0.1}'
advected_interp_method = 'upwind'
velocity_interp_method = 'rc'
T_cold = 301.15
T_hot = 311.15
Nx = 100
Ny = 50

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = vel_x
    v = vel_y
    pressure = pressure
  []
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 88.9e-3
    ymin = 0
    ymax = 63.5e-3
    nx = ${Nx}
    ny = ${Ny}
  []
[]

[AuxVariables]
  [U]
    type = MooseVariableFVReal
  []
  [fl]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
  [density]
    type = MooseVariableFVReal
  []
  [th_cond]
    type = MooseVariableFVReal
  []
  [cp_var]
    type = MooseVariableFVReal
  []
  [darcy_coef]
    type = MooseVariableFVReal
  []
  [fch_coef]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [mag]
    type = VectorMagnitudeAux
    variable = U
    x = vel_x
    y = vel_y
  []
  [compute_fl]
    type = NSLiquidFractionAux
    variable = fl
    temperature = T
    T_liquidus = '${T_liquidus}'
    T_solidus = '${T_solidus}'
    execute_on = 'TIMESTEP_END'
  []
  [rho_out]
    type = FunctorAux
    functor = 'rho_mixture'
    variable = 'density'
  []
  [th_cond_out]
    type = FunctorAux
    functor = 'k_mixture'
    variable = 'th_cond'
  []
  [cp_out]
    type = FunctorAux
    functor = 'cp_mixture'
    variable = 'cp_var'
  []
  [darcy_out]
    type = FunctorAux
    functor = 'Darcy_coefficient'
    variable = 'darcy_coef'
  []
  [fch_out]
    type = FunctorAux
    functor = 'Forchheimer_coefficient'
    variable = 'fch_coef'
  []
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
    initial_condition = 0.0
  []
  [vel_y]
    type = INSFVVelocityVariable
    initial_condition = 0.0
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [lambda]
    family = SCALAR
    order = FIRST
  []
  [T]
    type = INSFVEnergyVariable
    initial_condition = '${T_cold}'
    scaling = 1e-4
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = rho_mixture
  []
  [mean_zero_pressure]
    type = FVIntegralValueConstraint
    variable = pressure
    lambda = lambda
    phi0 = 0.0
  []

  [u_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_x
    rho = rho_mixture
    momentum_component = 'x'
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = rho_mixture
    momentum_component = 'x'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = vel_x
    momentum_component = 'x'
    pressure = pressure
  []
  [u_friction]
    type = PINSFVMomentumFriction
    variable = vel_x
    momentum_component = 'x'
    u = vel_x
    v = vel_y
    Darcy_name = 'Darcy_coeff'
    Forchheimer_name = 'Forchheimer_coeff'
    rho = ${rho_liquid}
    mu = ${mu}
    standard_friction_formulation = false
  []
  [u_buoyancy]
    type = INSFVMomentumBoussinesq
    variable = vel_x
    T_fluid = T
    gravity = '0 -9.81 0'
    rho = '${rho_liquid}'
    ref_temperature = ${T_cold}
    momentum_component = 'x'
  []
  [u_gravity]
    type = INSFVMomentumGravity
    variable = vel_x
    gravity = '0 -9.81 0'
    rho = '${rho_liquid}'
    momentum_component = 'x'
  []

  [v_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_y
    rho = rho_mixture
    momentum_component = 'y'
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = rho_mixture
    momentum_component = 'y'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = vel_y
    momentum_component = 'y'
    pressure = pressure
  []
  [v_friction]
    type = PINSFVMomentumFriction
    variable = vel_y
    momentum_component = 'y'
    u = vel_x
    v = vel_y
    Darcy_name = 'Darcy_coeff'
    Forchheimer_name = 'Forchheimer_coeff'
    rho = ${rho_liquid}
    mu = ${mu}
    standard_friction_formulation = false
  []
  [v_buoyancy]
    type = INSFVMomentumBoussinesq
    variable = vel_y
    T_fluid = T
    gravity = '0 -9.81 0'
    rho = '${rho_liquid}'
    ref_temperature = ${T_cold}
    momentum_component = 'y'
  []
  [v_gravity]
    type = INSFVMomentumGravity
    variable = vel_y
    gravity = '0 -9.81 0'
    rho = '${rho_liquid}'
    momentum_component = 'y'
  []

  [T_time]
    type = INSFVEnergyTimeDerivative
    variable = T
    rho = rho_mixture
    dh_dt = dh_dt
  []
  [energy_advection]
    type = INSFVEnergyAdvection
    variable = T
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
  []
  [energy_diffusion]
    type = FVDiffusion
    coeff = k_mixture
    variable = T
  []
  [energy_source]
    type = NSFVPhaseChangeSource
    variable = T
    L = ${L}
    liquid_fraction = fl
    T_liquidus = ${T_liquidus}
    T_solidus = ${T_solidus}
    rho = 'rho_mixture'
  []
[]

[FVBCs]
  [walls-u]
    type = INSFVNoSlipWallBC
    boundary = 'left right top bottom'
    variable = vel_x
    function = 0
  []
  [walls-v]
    type = INSFVNoSlipWallBC
    boundary = 'left right top bottom'
    variable = vel_y
    function = 0
  []
  [hot_wall]
    type = FVDirichletBC
    variable = T
    value = '${T_hot}'
    boundary = 'left'
  []
  [cold_wall]
    type = FVDirichletBC
    variable = T
    value = '${T_cold}'
    boundary = 'right'
  []
[]

[FunctorMaterials]
  [ins_fv]
    type = INSFVEnthalpyFunctorMaterial
    rho = rho_mixture
    cp = cp_mixture
    temperature = 'T'
  []
  [eff_cp]
    type = NSFVMixtureFunctorMaterial
    phase_2_names = '${cp_solid} ${k_solid} ${rho_solid}'
    phase_1_names = '${cp_liquid} ${k_liquid} ${rho_liquid}'
    prop_names = 'cp_mixture k_mixture rho_mixture'
    phase_1_fraction = fl
  []
  [mushy_zone_resistance]
    type = INSFVMushyPorousFrictionFunctorMaterial
    liquid_fraction = 'fl'
    mu = '${mu}'
    rho_l = '${rho_liquid}'
    dendrite_spacing_scaling = 1e-1
  []
  [friction]
    type = ADGenericVectorFunctorMaterial
    prop_names = 'Darcy_coeff Forchheimer_coeff'
    prop_values = 'darcy_coef darcy_coef darcy_coef fch_coef fch_coef fch_coef'
  []
  [const_functor]
    type = ADGenericFunctorMaterial
    prop_names = 'alpha_b'
    prop_values = '${alpha_b}'
  []
[]

[Executioner]
  type = Transient

  # Time-stepping parameters
  start_time = 0.0
  end_time = 200.0
  num_steps = 2

  [TimeStepper]
    type = IterationAdaptiveDT
    # Raise time step often but not by as much
    # There's a rough spot for convergence near 10% fluid fraction
    optimal_iterations = 15
    growth_factor = 1.5
    dt = 0.1
  []

  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  nl_rel_tol = 1e-6
  nl_max_its = 30
  line_search = 'none'
[]

[Postprocessors]
  [ave_p]
    type = ElementAverageValue
    variable = 'pressure'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [ave_fl]
    type = ElementAverageValue
    variable = 'fl'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [ave_T]
    type = ElementAverageValue
    variable = 'T'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[VectorPostprocessors]
  [vel_x]
    type = ElementValueSampler
    variable = 'vel_x fl'
    sort_by = 'x'
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
