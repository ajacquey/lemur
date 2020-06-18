[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 4
  nz = 1
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  zmin = 0
  zmax = 0.25
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[Kernels]
  [./mech_x]
    type = LMStressDivergence
    variable = disp_x
    component = 0
  [../]
  [./mech_y]
    type = LMStressDivergence
    variable = disp_y
    component = 1
  [../]
  [./mech_z]
    type = LMStressDivergence
    variable = disp_z
    component = 2
  [../]
[]

[AuxVariables]
  [./Se]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ed]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ed_v]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eta_e]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./Se_aux]
    type = LMVonMisesStressAux
    variable = Se
  [../]
  [./Ed_aux]
    type = LMEqvStrainAux
    variable = Ed
  [../]
  [./Ed_v_aux]
    type = LMEqvStrainAux
    variable = Ed_v
    strain_type = viscous
  [../]
  [./eta_e_aux]
    type = ADMaterialRealAux
    variable = eta_e
    property = effective_viscosity
  [../]
[]

[BCs]
  [./no_ux]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
    preset = true
  [../]
  [./ux_right]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = '-1.0e-14*t'
  [../]
  [./no_uy]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = 0.0
    preset = true
  [../]
  [./uy_bottom]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = bottom
    function = '-1.0e-14*t'
  [../]
  [./no_uz]
    type = DirichletBC
    variable = disp_z
    boundary = 'front back'
    value = 0.0
    preset = true
  [../]
[]

[Materials]
  [./elastic_mat]
    type = LMMechMaterial
    displacements = 'disp_x disp_y disp_z'
    bulk_modulus = 1.0e+10
    shear_modulus = 1.0e+10
    viscoelastic_model = 'maxwell'
  [../]
  [./maxwell]
    type = LMNonLinearViscosity
    viscosity = 1.0e+22
    exponent = 1.9
  [../]
[]

[Preconditioning]
  [./precond]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -ksp_max_it -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres asm 1E-15 1E-10 20 50 ilu NONZERO'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true
  start_time = 0.0
  end_time = 3.1536e+13
  dt = 3.1536e+11
[]

[Outputs]
  execute_on = 'TIMESTEP_END'
  print_linear_residuals = false
  perf_graph = true
  exodus = true
[]
