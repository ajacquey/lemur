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

#define usingViscoPlasticUpdateMembers                                                             \
  usingMaterialMembers;                                                                            \
  using LMViscoPlasticUpdate<compute_stage>::_abs_tol;                                             \
  using LMViscoPlasticUpdate<compute_stage>::_rel_tol;                                             \
  using LMViscoPlasticUpdate<compute_stage>::_max_its;                                             \
  using LMViscoPlasticUpdate<compute_stage>::_eta_p;                                               \
  using LMViscoPlasticUpdate<compute_stage>::_n;                                                   \
  using LMViscoPlasticUpdate<compute_stage>::_yield_function;                                      \
  using LMViscoPlasticUpdate<compute_stage>::_plastic_strain_incr

template <ComputeStage>
class LMViscoPlasticUpdate;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;

declareADValidParams(LMViscoPlasticUpdate);

template <ComputeStage compute_stage>
class LMViscoPlasticUpdate : public ADMaterial<compute_stage>
{
public:
  LMViscoPlasticUpdate(const InputParameters & parameters);
  void setQp(unsigned int qp);
  virtual void viscoPlasticUpdate(ADRankTwoTensor & stress,
                                  const RankFourTensor & Cijkl,
                                  ADRankTwoTensor & elastic_strain_incr) = 0;
  void resetQpProperties() final {}
  void resetProperties() final {}

protected:
  Real _abs_tol;
  Real _rel_tol;
  unsigned int _max_its;
  Real _eta_p;
  Real _n;

  ADMaterialProperty(Real) & _yield_function;
  ADMaterialProperty(RankTwoTensor) & _plastic_strain_incr;

  usingMaterialMembers;
};