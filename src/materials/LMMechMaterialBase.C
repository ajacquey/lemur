/******************************************************************************/
/*                            This file is part of                            */
/*                       LEMUR, a MOOSE-based application                     */
/*          muLtiphysics of gEomaterials using MUltiscale Rheologies          */
/*                                                                            */
/*                  Copyright (C) 2020 by Antoine B. Jacquey                  */
/*                    Massachusetts Institute of Technology                   */
/*                                                                            */
/*            Licensed under GNU Lesser General Public License v2.1           */
/*                       please see LICENSE for details                       */
/*                 or http://www.gnu.org/licenses/lgpl.html                   */
/******************************************************************************/

#include "LMMechMaterialBase.h"
#include "LMViscoElasticUpdate.h"
#include "LMViscoPlasticUpdate.h"
#include "Function.h"

InputParameters
LMMechMaterialBase::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("Base class calculating the strain and stress of a material.");
  // Coupled variables
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system.");
  // Strain parameters
  MooseEnum strain_model("small=0 finite=1", "small");
  params.addParam<MooseEnum>(
      "strain_model", strain_model, "The model to use to calculate the strain rate tensor.");
  // Initial stress
  params.addParam<std::vector<FunctionName>>(
      "initial_stress", "The initial stress principal components (negative in compression).");
  // Visco-Elastic model
  params.addParam<MaterialName>("viscoelastic_model",
                                "The material object to use for the viscoelastic correction.");
  // Visco-Plastic model
  params.addParam<MaterialName>("viscoplastic_model",
                                "The material object to use for the viscoplastic correction.");
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

LMMechMaterialBase::LMMechMaterialBase(const InputParameters & parameters)
  : ADMaterial(parameters),
    // Coupled variables
    _ndisp(coupledComponents("displacements")),
    _grad_disp(3),
    _grad_disp_old(3),
    // Strain parameters
    _strain_model(getParam<MooseEnum>("strain_model")),
    // Initial stress
    _initial_stress_fct(getParam<std::vector<FunctionName>>("initial_stress")),
    _num_ini_stress(_initial_stress_fct.size()),
    // Visco-Elastic model
    _has_ve(isParamValid("viscoelastic_model")),
    // Visco-Plastic model
    _has_vp(isParamValid("viscoplastic_model")),
    // Strain properties
    _strain_increment(declareADProperty<RankTwoTensor>("strain_increment")),
    _spin_increment(declareADProperty<RankTwoTensor>("spin_increment")),
    _elastic_strain_incr(declareADProperty<RankTwoTensor>("elastic_strain_increment")),
    // Stress properties
    _K(declareADProperty<Real>("bulk_modulus")),
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

void
LMMechMaterialBase::initialSetup()
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

  // Fetch viscoelastic model object
  if (_has_ve)
  {
    MaterialName ve_model = getParam<MaterialName>("viscoelastic_model");

    LMViscoElasticUpdate * ve_r =
        dynamic_cast<LMViscoElasticUpdate *>(&this->getMaterialByName(ve_model));

    _ve_model = ve_r;
  }
  else
    _ve_model = nullptr;

  // Fetch viscoplastic model object
  if (_has_vp)
  {
    MaterialName vp_model = getParam<MaterialName>("viscoplastic_model");

    LMViscoPlasticUpdate * vp_r =
        dynamic_cast<LMViscoPlasticUpdate *>(&this->getMaterialByName(vp_model));

    _vp_model = vp_r;
  }
  else
    _vp_model = nullptr;
}

void
LMMechMaterialBase::displacementIntegrityCheck()
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    paramError(
        "displacements",
        "The number of variables supplied in 'displacements' must match the mesh dimension.");
}

void
LMMechMaterialBase::initQpStatefulProperties()
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

void
LMMechMaterialBase::computeQpProperties()
{
  computeQpStrainIncrement();
  computeQpElasticityTensor();
  computeQpStress();
}

void
LMMechMaterialBase::computeQpStrainIncrement()
{
  ADRankTwoTensor grad_tensor = ADRankTwoTensor::initializeFromRows(
      (*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);
  RankTwoTensor grad_tensor_old = RankTwoTensor::initializeFromRows(
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

void
LMMechMaterialBase::computeQpSmallStrain(const ADRankTwoTensor & grad_tensor,
                                         const RankTwoTensor & grad_tensor_old)
{
  ADRankTwoTensor A = grad_tensor - grad_tensor_old;

  _strain_increment[_qp] = 0.5 * (A + A.transpose());
  _spin_increment[_qp] = 0.5 * (A - A.transpose());
}

void
LMMechMaterialBase::computeQpFiniteStrain(const ADRankTwoTensor & grad_tensor,
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

void
LMMechMaterialBase::computeQpStress()
{
  // Elastic guess
  computeQpElasticGuess();

  // Viscoelastic correction
  if (_has_ve)
  {
    _ve_model->setQp(_qp);
    _ve_model->viscoElasticUpdate(_stress[_qp], _Cijkl, _elastic_strain_incr[_qp]);
  }

  // Viscoplastic correction
  if (_has_vp)
  {
    _vp_model->setQp(_qp);
    _vp_model->viscoPlasticUpdate(_stress[_qp], _Cijkl, _elastic_strain_incr[_qp]);
  }
}

void
LMMechMaterialBase::computeQpElasticGuess()
{
  _elastic_strain_incr[_qp] = _strain_increment[_qp];
  _stress[_qp] = spinRotation(_stress_old[_qp]) + _Cijkl * _strain_increment[_qp];
}

ADRankTwoTensor
LMMechMaterialBase::spinRotation(const ADRankTwoTensor & tensor)
{
  return tensor + _spin_increment[_qp] * tensor.deviatoric() -
         tensor.deviatoric() * _spin_increment[_qp];
}