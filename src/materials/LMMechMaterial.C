/******************************************************************************/
/*                            This file is part of                            */
/*                       LEMUR, a MOOSE-based application                     */
/*          muLtiphysics of gEomaterials using MUltiscale Rheologies          */
/*                                                                            */
/*                  Copyright (C) 2019 by Antoine B. Jacquey                  */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*            Licensed under GNU Lesser General Public License v2.1           */
/*                       please see LICENSE for details                       */
/*                 or http://www.gnu.org/licenses/lgpl.html                   */
/******************************************************************************/

#include "LMMechMaterial.h"
#include "LMViscoPlasticUpdate.h"
#include "Function.h"

registerADMooseObject("LemurApp", LMMechMaterial);

defineADValidParams(
    LMMechMaterial,
    ADMaterial,
    params.addClassDescription("Base class calculating the strain and stress of a material.");
    // Coupled variables
    params.addRequiredCoupledVar(
        "displacements",
        "The displacements appropriate for the simulation geometry and coordinate system.");
    // Strain parameters
    MooseEnum strain_model("small=0 finite=1", "small");
    params.addParam<MooseEnum>("strain_model",
                               strain_model,
                               "The model to use to calculate the strain rate tensor.");
    // Elastic moduli parameters
    params.addRequiredRangeCheckedParam<Real>("bulk_modulus",
                                              "bulk_modulus > 0.0",
                                              "The bulk modulus of the material.");
    params.addRequiredRangeCheckedParam<Real>("shear_modulus",
                                              "shear_modulus > 0.0",
                                              "The shear modulus of the material.");
    // Initial stress
    params.addParam<std::vector<FunctionName>>(
        "initial_stress",
        "The initial stress principal components (negative in compression).");
    // Visco-Plastic model
    params.addParam<MaterialName>("viscoplastic_model",
                                  "The material object to use for the viscoplastic correction.");
    params.suppressParameter<bool>("use_displaced_mesh"););

template <ComputeStage compute_stage>
LMMechMaterial<compute_stage>::LMMechMaterial(const InputParameters & parameters)
  : ADMaterial<compute_stage>(parameters),
    // Coupled variables
    _ndisp(coupledComponents("displacements")),
    _grad_disp(3),
    _grad_disp_old(3),
    // Strain parameters
    _strain_model(getParam<MooseEnum>("strain_model")),
    // Elastic moduli parameters
    _bulk_modulus(getParam<Real>("bulk_modulus")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    // Initial stress
    _initial_stress_fct(getParam<std::vector<FunctionName>>("initial_stress")),
    _num_ini_stress(_initial_stress_fct.size()),
    // Visco-Plastic model
    _has_vp(isParamValid("viscoplastic_model")),
    // Strain properties
    _strain_increment(declareADProperty<RankTwoTensor>("strain_increment")),
    _spin_increment(declareADProperty<RankTwoTensor>("spin_increment")),
    _elastic_strain_incr(declareADProperty<RankTwoTensor>("elastic_strain_increment")),
    // Stress properties
    _stress(declareADProperty<RankTwoTensor>("stress")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>("stress"))
{
  if (getParam<bool>("use_displaced_mesh"))
    paramError("use_displaced_mesh",
               "The strain and stress calculator needs to run on the undisplaced mesh.");

  if (_num_ini_stress != 3 && _num_ini_stress != 0)
    paramError("initial_stress", "You need to provide 3 components for the initial stress.");
  
  _initial_stress.resize(_num_ini_stress);

  for (unsigned int i = 0; i < _num_ini_stress; i++)
    _initial_stress[i] = &getFunctionByName(_initial_stress_fct[i]);
}

template <ComputeStage compute_stage>
void
LMMechMaterial<compute_stage>::initialSetup()
{
  displacementIntegrityCheck();
  // Fetch coupled variables and gradients
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _grad_disp[i] = &adCoupledGradient("displacements", i);
    if (_fe_problem.isTransient())
      _grad_disp_old[i] = &coupledGradientOld("displacements", i);
    else
      _grad_disp_old[i] = &_grad_zero;
  }

  // Set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _grad_disp[i] = &adZeroGradient();
    _grad_disp_old[i] = &_grad_zero;
  }

  // Fetch viscoplastic model object
  if (_has_vp)
  {
    MaterialName vp_model = getParam<MaterialName>("viscoplastic_model");

    LMViscoPlasticUpdate<compute_stage> * vp_r =
        dynamic_cast<LMViscoPlasticUpdate<compute_stage> *>(
            &this->template getMaterialByName<compute_stage>(vp_model));

    _vp_model = vp_r;
  }
  else
    _vp_model = nullptr;
}

template <ComputeStage compute_stage>
void
LMMechMaterial<compute_stage>::displacementIntegrityCheck()
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    paramError(
        "displacements",
        "The number of variables supplied in 'displacements' must match the mesh dimension.");
}

template <ComputeStage compute_stage>
void
LMMechMaterial<compute_stage>::initQpStatefulProperties()
{
  _stress[_qp].zero();
  RankTwoTensor init_stress_tensor = RankTwoTensor();
  if (_num_ini_stress > 0)
  {
    std::vector<Real> init_stress(_num_ini_stress, 0.0);
    for (unsigned int i = 0; i < _num_ini_stress; i++)
      init_stress[i] = (*_initial_stress[i]).value(_t, _q_point[_qp]);
    init_stress_tensor.fillFromInputVector(init_stress);
  }
  _stress[_qp] += init_stress_tensor;
}

template <ComputeStage compute_stage>
void
LMMechMaterial<compute_stage>::computeQpProperties()
{
  computeQpStrainIncrement();
  computeQpElasticityTensor();
  computeQpStress();
}

template <ComputeStage compute_stage>
void
LMMechMaterial<compute_stage>::computeQpStrainIncrement()
{
  ADRankTwoTensor grad_tensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);
  RankTwoTensor grad_tensor_old(
      (*_grad_disp_old[0])[_qp], (*_grad_disp_old[1])[_qp], (*_grad_disp_old[2])[_qp]);

  switch (_strain_model)
  {
    case 0: // SMALL STRAIN
      computeQpSmallStrain(grad_tensor, grad_tensor_old);
      break;
    case 1: // FINITE STRAIN
      computeQpFiniteStrain(grad_tensor, grad_tensor_old);
      break;
    default:
      mooseError("Unknown strain model. Specify 'small' or 'finite'!");
  }
}

template <ComputeStage compute_stage>
void
LMMechMaterial<compute_stage>::computeQpSmallStrain(const ADRankTwoTensor & grad_tensor,
                                                    const RankTwoTensor & grad_tensor_old)
{
  ADRankTwoTensor A = grad_tensor - grad_tensor_old;

  _strain_increment[_qp] = 0.5 * (A + A.transpose());
  _spin_increment[_qp] = 0.5 * (A - A.transpose());
}

template <ComputeStage compute_stage>
void
LMMechMaterial<compute_stage>::computeQpFiniteStrain(const ADRankTwoTensor & grad_tensor,
                                                     const RankTwoTensor & grad_tensor_old)
{
  ADRankTwoTensor F = grad_tensor;
  RankTwoTensor F_old = grad_tensor_old;
  F.addIa(1.0);
  F_old.addIa(1.0);

  // Increment gradient
  ADRankTwoTensor L = -F_old * F.inverse();
  L.addIa(1.0);

  _strain_increment[_qp] = 0.5 * (L + L.transpose());
  _spin_increment[_qp] = 0.5 * (L - L.transpose());
}

template <ComputeStage compute_stage>
void
LMMechMaterial<compute_stage>::computeQpElasticityTensor()
{
  _Cijkl.fillGeneralIsotropic(_bulk_modulus - 2.0 / 3.0 * _shear_modulus, _shear_modulus, 0.0);
}

template <ComputeStage compute_stage>
void
LMMechMaterial<compute_stage>::computeQpStress()
{
  // Elastic guess
  computeQpElasticGuess();

  // Viscoplastic correction
  if (_has_vp)
  {
    _vp_model->setQp(_qp);
    _vp_model->viscoPlasticUpdate(_stress[_qp], _Cijkl, _elastic_strain_incr[_qp]);
  }
}

template <ComputeStage compute_stage>
void
LMMechMaterial<compute_stage>::computeQpElasticGuess()
{
  _elastic_strain_incr[_qp] = _strain_increment[_qp];
  _stress[_qp] = spinRotation(_stress_old[_qp]) + _Cijkl * _strain_increment[_qp];
}

template <ComputeStage compute_stage>
ADRankTwoTensor
LMMechMaterial<compute_stage>::spinRotation(const ADRankTwoTensor & tensor)
{
  return tensor + _spin_increment[_qp] * tensor.deviatoric() -
         tensor.deviatoric() * _spin_increment[_qp];
}

adBaseClass(LMMechMaterial);