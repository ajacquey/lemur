# Terzaghi's problem of consolodation of a drained medium
#
# See Arnold Verruijt "Theory and Problems of Poroelasticity" 2015
# Section 2.2 Terzaghi's problem

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  # nz = 100
  nz = 10
  xmin = -0.1
  xmax = 0.1
  ymin = -0.1
  ymax = 1
  zmin = 0
  zmax = 1
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
  [./pf]
  [../]
[]

[Kernels]
  [./grad_stress_x]
    type = LMStressDivergence
    variable = disp_x
    fluid_pressure = pf
    component = 0
  [../]
  [./grad_stress_y]
    type = LMStressDivergence
    variable = disp_y
    fluid_pressure = pf
    component = 1
  [../]
  [./grad_stress_z]
    type = LMStressDivergence
    variable = disp_z
    fluid_pressure = pf
    component = 2
  [../]
  [./pf_time_derivative]
    type = LMFluidFlowTimeDerivative
    variable = pf
  [../]
  [./darcy]
    type = LMFluidFlowDarcy
    variable = pf
  [../]
[]

[AuxVariables]
  [./phi]
    initial_condition = 0.1
  [../]
[]

[AuxKernels]
  [./phi_aux]
    type = ConstantAux
    variable = phi
    value = 0.1
  [../]
[]

[BCs]
  [./confinex]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'left right'
    preset = true
  [../]
  [./confiney]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'bottom top'
    preset = true
  [../]
  [./basefixed]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = back
    preset = true
  [../]
  [./topdrained]
    type = DirichletBC
    variable = pf
    value = 0
    boundary = front
  [../]
  [./topload]
    type = NeumannBC
    variable = disp_z
    value = -1
    boundary = front
  [../]
[]

[Materials]
  [./mechanical]
    type = LMMechMaterial
    displacements = 'disp_x disp_y disp_z'
    bulk_modulus = 4
    shear_modulus = 3
  [../]
  [./hydraulic]
    type = LMPoroMaterial
    porosity = phi
    permeability = 1.5e-02
    fluid_viscosity = 1.395348837e-01
    fluid_modulus = 8
    solid_modulus = 10
  [../]
[]

[VectorPostprocessors]
  [./line_pf]
    type = LineValueSampler
    variable = pf
    start_point = '0.0 0.0 0.0'
    end_point = '0.0 0.0 1.0'
    num_points = 10
    sort_by = 'z'
    outputs = 'csv'
  [../]
[]

[Preconditioning]
  [./hypre]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew -snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-pc_type -pc_hypre_type
                           -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_max_levels
                           -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type
                           -pc_hypre_boomeramg_P_max -pc_hypre_boomeramg_truncfacto -snes_atol'
    petsc_options_value = 'hypre boomeramg
                           0.7 4 5 25
                           HMIS ext+i
                           2 0.3 1.0e-14'
  [../]
[]

[Functions]
  # [./time_stepper_fct]
  #   type = PiecewiseConstant
  #   x = '0      0.01  0.1  1.0'
  #   y = '0.0001 0.001 0.01 0.1'
  # [../]
  [./time_stepper_fct]
    type = PiecewiseConstant
    x = '0      0.01  0.1'
    y = '0.001 0.01 0.1'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  # automatic_scaling = true
  start_time = 0
  end_time = 10 # ~10 s
  [./TimeStepper]
    type = FunctionDT
    function = time_stepper_fct
  [../]
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = true
  execute_on = 'TIMESTEP_END'
  exodus = true
  [./csv]
    type = CSV
    sync_only = true
    sync_times = '0.001 0.01 0.05 0.1 0.5 1.0'
  [../]
[]
