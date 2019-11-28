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

#pragma once

#include "ADMaterial.h"

template <ComputeStage>
class LMMechMaterial;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;
template <ComputeStage>
class LMViscoPlasticUpdate;

declareADValidParams(LMMechMaterial);

template <ComputeStage compute_stage>
class LMMechMaterial : public ADMaterial<compute_stage>
{
public:
  LMMechMaterial(const InputParameters & parameters);
  void initialSetup() override;
  void displacementIntegrityCheck();

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  virtual void computeQpStrainIncrement();
  virtual void computeQpSmallStrain(const ADRankTwoTensor & grad_tensor,
                                    const RankTwoTensor & grad_tensor_old);
  virtual void computeQpFiniteStrain(const ADRankTwoTensor & grad_tensor,
                                     const RankTwoTensor & grad_tensor_old);
  virtual void computeQpElasticityTensor();
  virtual void computeQpStress();
  virtual ADRankTwoTensor spinRotation(const ADRankTwoTensor & tensor);

  // Coupled variables
  const unsigned int _ndisp;
  std::vector<const ADVariableGradient *> _grad_disp;
  std::vector<const VariableGradient *> _grad_disp_old;

  // Strain parameters
  const unsigned int _strain_model;

  // Elastic parameters
  const Real _bulk_modulus;
  const Real _shear_modulus;

  // Initial stress
  const std::vector<FunctionName> _initial_stress_fct;
  const unsigned int _num_ini_stress;

  // Viscoplastoc model
  const bool _has_vp;

  // Strain properties
  ADMaterialProperty(RankTwoTensor) & _strain_increment;
  ADMaterialProperty(RankTwoTensor) & _spin_increment;
  ADMaterialProperty(RankTwoTensor) & _elastic_strain_incr;

  // Stress properties
  ADMaterialProperty(RankTwoTensor) & _stress;
  const MaterialProperty<RankTwoTensor> & _stress_old;

  // Initial stresses
  std::vector<const Function *> _initial_stress;
  
  // Viscoplastic model
  LMViscoPlasticUpdate<compute_stage> * _vp_model;

  // Elasticity tensor
  RankFourTensor _Cijkl;

  usingMaterialMembers;
};