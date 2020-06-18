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
    _has_ve(hasADMaterialProperty<Real>("effective_viscosity")),
    _viscous_strain_incr(_has_ve ? &getADMaterialProperty<RankTwoTensor>("viscous_strain_increment")
                                 : nullptr),
    _has_vp(hasADMaterialProperty<Real>("yield_function")),
    _plastic_strain_incr(_has_vp ? &getADMaterialProperty<RankTwoTensor>("plastic_strain_increment")
                                 : nullptr)
{
}

Real
LMPorosityAux::computeValue()
{
  Real ev_ve_incr = _has_ve ? MetaPhysicL::raw_value((*_viscous_strain_incr)[_qp].trace()) : 0.0;
  Real ev_vp_incr = _has_vp ? MetaPhysicL::raw_value((*_plastic_strain_incr)[_qp].trace()) : 0.0;

  return _u_old[_qp] +
         (_biot[_qp] - _u[_qp]) * ((1.0 - _biot[_qp]) / _K[_qp] * _pf_dot[_qp] * _dt +
                                   MetaPhysicL::raw_value(_strain_incr[_qp].trace())) +
         (1.0 - _biot[_qp]) * ev_ve_incr + (1.0 - _biot[_qp]) * ev_vp_incr;
}