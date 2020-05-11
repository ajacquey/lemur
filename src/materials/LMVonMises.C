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

#include "LMVonMises.h"

registerMooseObject("LemurApp", LMVonMises);

InputParameters
LMVonMises::validParams()
{
  InputParameters params = LMSingleVarUpdate::validParams();
  params.addClassDescription("Viscoplastic update based on a Von Mises yield function.");
  params.addCoupledVar("coupled_var", "A coupled variable on which the yield depends.");
  params.addRequiredRangeCheckedParam<Real>(
      "yield_strength", "yield_strength >= 0.0", "The yield strength.");
  params.addParam<Real>("hardening_modulus", 0.0, "The hardening modulus of the Von Mises yield.");
  params.addParam<Real>("coupled_hardening_modulus",
                        0.0,
                        "The hardening modulus of the coupled variable for the Von Mises yield.");
  return params;
}

LMVonMises::LMVonMises(const InputParameters & parameters)
  : LMSingleVarUpdate(parameters),
    _coupled_v(isCoupled("coupled_var")),
    _v(_coupled_v ? adCoupledValue("coupled_var") : adZeroValue()),
    _yield_strength(getParam<Real>("yield_strength")),
    _hg(getParam<Real>("hardening_modulus")),
    _ht(getParam<Real>("coupled_hardening_modulus")),
    _has_hardening(_hg != 0.0),
    _intnl(_has_hardening ? &declareADProperty<Real>("deviatoric_plastic_strain") : nullptr),
    _intnl_old(_has_hardening ? &getMaterialPropertyOld<Real>("deviatoric_plastic_strain")
                              : nullptr)
{
}

void
LMVonMises::initQpStatefulProperties()
{
  if (_has_hardening)
    (*_intnl)[_qp] = 0.0;
}

ADReal
LMVonMises::yieldFunction(const ADReal & gamma_vp)
{
  return (_tau_tr - 2.0 * _G * gamma_vp * _dt) - (_yield_strength_tr + _hg * gamma_vp * _dt);
}

ADReal
LMVonMises::yieldFunctionDeriv(const ADReal & /*gamma_vp*/)
{
  return -(2.0 * _G + _hg) * _dt;
}

void
LMVonMises::preReturnMap()
{
  _tau_tr = std::sqrt(0.5) * _stress_tr.deviatoric().L2norm();

  _yield_strength_tr = _yield_strength;
  if (_has_hardening)
  {
    (*_intnl)[_qp] = (*_intnl_old)[_qp];
    _yield_strength_tr = _yield_strength + _hg * (*_intnl)[_qp] + _ht * _v[_qp];
  }
}

void
LMVonMises::postReturnMap(const ADReal & gamma_vp)
{
  if (_has_hardening)
    (*_intnl)[_qp] = (*_intnl_old)[_qp] + gamma_vp * _dt;
}

ADRankTwoTensor
LMVonMises::reformPlasticStrainTensor(const ADReal & gamma_vp)
{
  ADRankTwoTensor flow_dir =
      (_tau_tr != 0.0) ? _stress_tr.deviatoric() / _tau_tr : ADRankTwoTensor();

  return gamma_vp * _dt * flow_dir;
}
