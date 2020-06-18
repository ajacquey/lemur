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

#include "LMVelocityAux.h"

registerMooseObject("LemurApp", LMVelocityAux);

InputParameters
LMVelocityAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the velocity based on the displacement.");
  params.addCoupledVar("displacement", "The displacement component.");
  return params;
}

LMVelocityAux::LMVelocityAux(const InputParameters & parameters)
  : AuxKernel(parameters), _disp_dot(coupledDot("displacement"))
{
}

Real
LMVelocityAux::computeValue()
{
  return _disp_dot[_qp];
}