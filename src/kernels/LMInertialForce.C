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

#include "LMInertialForce.h"

registerMooseObject("LemurApp", LMInertialForce);

InputParameters
LMInertialForce::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addClassDescription("Inertial term for dynamic deformation.");
  params.set<bool>("use_displaced_mesh") = true;
  params.addRangeCheckedParam<Real>(
      "density", 1.0, "density >= 0.0", "The density of the material.");
  return params;
}

LMInertialForce::LMInertialForce(const InputParameters & parameters)
  : TimeKernel(parameters),
    _u_dot_dot(dotDot()),
    _du_dot_dot_du(dotDotDu()),
    _rho(getParam<Real>("density"))
{
}

Real
LMInertialForce::computeQpResidual()
{
  return _rho * _u_dot_dot[_qp] * _test[_i][_qp];
}

Real
LMInertialForce::computeQpJacobian()
{
  return _rho * _du_dot_dot_du[_qp] * _test[_i][_qp];
}