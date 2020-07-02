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

#include "LMPoroMechMaterial.h"

registerMooseObject("LemurApp", LMPoroMechMaterial);

InputParameters
LMPoroMechMaterial::validParams()
{
  InputParameters params = LMMechMaterialBase::validParams();
  params.addClassDescription("Base class calculating the strain and stress of a porous material "
                             "based on non-linear elasticity.");
  // Coupled variables
  params.addRequiredCoupledVar("porosity", "The porosity variable.");
  // Elastic moduli parameters
  params.addRequiredRangeCheckedParam<Real>("compression_index",
                                            "compression_index > 0.0",
                                            "The compression index of the porous material.");
  params.addRequiredRangeCheckedParam<Real>(
      "poisson_ratio", "poisson_ratio > 0.0", "The poisson ratio of the porous material.");
  return params;
}

LMPoroMechMaterial::LMPoroMechMaterial(const InputParameters & parameters)
  : LMMechMaterialBase(parameters),
    // Coupled variables
    _porosity(coupledValue("porosity")),
    // Elastic moduli parameters
    _compression_idx(getParam<Real>("compression_index")),
    _poisson_ratio(getParam<Real>("poisson_ratio")),
    // Strain quantities
    _elastic_strain_incr_old(getMaterialPropertyOld<RankTwoTensor>("elastic_strain_increment"))
{
}

void
LMPoroMechMaterial::initQpStatefulProperties()
{
  LMMechMaterialBase::initQpStatefulProperties();

  if (_stress[_qp].trace() == 0.0)
    mooseWarning("LMPoroMechMaterial: Initial pressure is zero!");
}

void
LMPoroMechMaterial::computeQpElasticityTensor()
{
  RankTwoTensor I = RankTwoTensor(RankTwoTensor::initIdentity);
  RankFourTensor I4 = RankFourTensor(RankFourTensor::initIdentitySymmetricFour);

  Real r = 3.0 * (1.0 - 2.0 * _poisson_ratio) / (2.0 * (1.0 + _poisson_ratio));
  RankFourTensor c_hat = I.outerProduct(I) + 2.0 * r * (I4 - I.outerProduct(I) / 3.0);

  ADReal p_old = -_stress_old[_qp].trace() / 3.0;
  ADReal strain_vol_incr = -_strain_increment[_qp].trace();
  ADReal strain_el_vol_incr = -_elastic_strain_incr_old[_qp].trace();

  // if (MooseUtils::absoluteFuzzyEqual(strain_vol_incr, 0.0))
  //   _K[_qp] = 1.0;
  // else
  //   _K[_qp] = p_old / strain_vol_incr * (std::exp(strain_vol_incr / (_compression_idx * (1.0 +
  //   _porosity[_qp]))) - 1.0);
  _K[_qp] = p_old / (_compression_idx * (1.0 + _porosity[_qp]));

  _Cijkl = _K[_qp] * c_hat;
}