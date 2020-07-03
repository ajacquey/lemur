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

#include "LMVelocityBC.h"

registerMooseObject("LemurApp", LMVelocityBC);

InputParameters
LMVelocityBC::validParams()
{
  InputParameters params = DirichletBCBase::validParams();
  params.addClassDescription("Applies a velocity whose value is described by a function.");
  params.addParam<Real>("value", 0.0, "Value of the velocity applied.");
  params.set<bool>("preset") = true;
  return params;
}

LMVelocityBC::LMVelocityBC(const InputParameters & parameters)
  : DirichletBCBase(parameters), _u_old(valueOld()), _value(getParam<Real>("value"))
{
}

Real
LMVelocityBC::computeQpValue()
{
  return _u_old[_qp] + _value * _dt;
}