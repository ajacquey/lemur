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

#include "LMMechMaterial.h"
#include "LMViscoElasticUpdate.h"
#include "LMViscoPlasticUpdate.h"
#include "Function.h"

registerMooseObject("LemurApp", LMMechMaterial);

InputParameters
LMMechMaterial::validParams()
{
  InputParameters params = LMMechMaterialBase::validParams();
  params.addClassDescription(
      "Class calculating the strain and stress of a material using constant elastic moduli.");
  // Elastic moduli parameters
  params.addRequiredRangeCheckedParam<Real>(
      "bulk_modulus", "bulk_modulus > 0.0", "The bulk modulus of the material.");
  params.addRequiredRangeCheckedParam<Real>(
      "shear_modulus", "shear_modulus > 0.0", "The shear modulus of the material.");
  return params;
}

LMMechMaterial::LMMechMaterial(const InputParameters & parameters)
  : LMMechMaterialBase(parameters),
    // Elastic moduli parameters
    _bulk_modulus(getParam<Real>("bulk_modulus")),
    _shear_modulus(getParam<Real>("shear_modulus"))
{
}

void
LMMechMaterial::computeQpElasticityTensor()
{
  // Bulk modulus
  _K[_qp] = _bulk_modulus;
  // Elasticity tensor
  _Cijkl.fillGeneralIsotropic(_bulk_modulus - 2.0 / 3.0 * _shear_modulus, _shear_modulus, 0.0);
}