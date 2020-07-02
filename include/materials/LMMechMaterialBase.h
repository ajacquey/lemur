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

#pragma once

#include "ADMaterial.h"

class LMViscoElasticUpdate;
class LMViscoPlasticUpdate;

class LMMechMaterialBase : public ADMaterial
{
public:
  static InputParameters validParams();
  LMMechMaterialBase(const InputParameters & parameters);
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
  virtual void computeQpElasticityTensor() = 0;
  virtual void computeQpStress();
  virtual void computeQpElasticGuess();
  virtual ADRankTwoTensor spinRotation(const ADRankTwoTensor & tensor);

  // Coupled variables
  const unsigned int _ndisp;
  std::vector<const ADVariableGradient *> _grad_disp;
  std::vector<const VariableGradient *> _grad_disp_old;

  // Strain parameters
  const unsigned int _strain_model;

  // Initial stress
  const std::vector<FunctionName> _initial_stress_fct;
  const unsigned int _num_ini_stress;

  // Viscoelastic model
  const bool _has_ve;

  // Viscoplastic model
  const bool _has_vp;

  // Strain properties
  ADMaterialProperty<RankTwoTensor> & _strain_increment;
  ADMaterialProperty<RankTwoTensor> & _spin_increment;
  ADMaterialProperty<RankTwoTensor> & _elastic_strain_incr;

  // Stress properties
  ADMaterialProperty<Real> & _K;
  ADMaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _stress_old;

  // Initial stresses
  std::vector<const Function *> _initial_stress;

  // Viscoelastic model
  LMViscoElasticUpdate * _ve_model;

  // Viscoplastic model
  LMViscoPlasticUpdate * _vp_model;

  // Elasticity tensor
  ADRankFourTensor _Cijkl;
};