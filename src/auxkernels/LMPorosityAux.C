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

#include "LMPorosityAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LemurApp", LMPorosityAux);

InputParameters
LMPorosityAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the porosity.");
  params.addCoupledVar("fluid_pressure", "The fluid pressure variable.");
  return params;
}

LMPorosityAux::LMPorosityAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _pf_dot(coupledDot("fluid_pressure")),
    _biot(getMaterialProperty<Real>("biot_coefficient")),
    _K(getMaterialProperty<Real>("bulk_modulus")),
    _strain_incr(getADMaterialProperty<RankTwoTensor>("strain_increment")),
    _plastic_strain_incr(getADMaterialProperty<RankTwoTensor>("plastic_strain_increment"))
{
}

Real
LMPorosityAux::computeValue()
{
  return _u_old[_qp] +
         (_biot[_qp] - _u[_qp]) * ((1.0 - _biot[_qp]) / _K[_qp] * _pf_dot[_qp] * _dt +
                                   MetaPhysicL::raw_value(_strain_incr[_qp].trace())) +
         (1.0 - _biot[_qp]) * MetaPhysicL::raw_value(_plastic_strain_incr[_qp].trace());
}