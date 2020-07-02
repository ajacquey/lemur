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

#include "LMPoroMaterial.h"

registerMooseObject("LemurApp", LMPoroMaterial);

InputParameters
LMPoroMaterial::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("Computes properties for fluid flow in a porous material.");
  params.addCoupledVar("porosity", 0.0, "The porosity variable.");
  params.addRequiredRangeCheckedParam<Real>(
      "permeability", "permeability > 0.0", "The permeability of the material.");
  params.addRequiredRangeCheckedParam<Real>(
      "fluid_viscosity", "fluid_viscosity > 0.0", "The fluid viscosity.");
  params.addRangeCheckedParam<Real>(
      "fluid_modulus", "fluid_modulus > 0.0", "The fluid bulk modulus.");
  params.addRangeCheckedParam<Real>(
      "solid_modulus", "solid_modulus > 0.0", "The solid bulk modulus.");
  return params;
}

LMPoroMaterial::LMPoroMaterial(const InputParameters & parameters)
  : ADMaterial(parameters),
    _porosity(coupledValue("porosity")),
    _coupled_mech(hasADMaterialProperty<Real>("bulk_modulus")),
    _perm(getParam<Real>("permeability")),
    _fluid_visco(getParam<Real>("fluid_viscosity")),
    _Kf(isParamValid("fluid_modulus") ? getParam<Real>("fluid_modulus") : 0.0),
    _Ks(isParamValid("solid_modulus") ? getParam<Real>("solid_modulus") : 0.0),
    _K(_coupled_mech ? &getADMaterialProperty<Real>("bulk_modulus") : nullptr),
    _strain_increment(_coupled_mech ? &getADMaterialProperty<RankTwoTensor>("strain_increment")
                                    : nullptr),
    _has_ve(hasADMaterialProperty<Real>("effective_viscosity")),
    _viscous_strain_incr(_has_ve ? &getADMaterialProperty<RankTwoTensor>("viscous_strain_increment")
                                 : nullptr),
    _has_vp(hasADMaterialProperty<Real>("yield_function")),
    _plastic_strain_incr(_has_vp ? &getADMaterialProperty<RankTwoTensor>("plastic_strain_increment")
                                 : nullptr),
    _C_biot(declareADProperty<Real>("biot_compressibility")),
    _fluid_mob(declareProperty<Real>("fluid_mobility")),
    _biot(declareADProperty<Real>("biot_coefficient")),
    _poro_mech(declareADProperty<Real>("poro_mech"))
{
  if (_fe_problem.isTransient() && _coupled_mech && (_Ks == 0.0))
    mooseWarning(
        "LMPoroMaterial: running a transient hydro-mechanical simulation but did not supplied "
        "solid_modulus!");
  if (_fe_problem.isTransient() && (_Kf == 0.0))
    mooseWarning("LMPoroMaterial: running a transient simulation but did not supplied "
                 "fluid_modulus!");
  if (_fe_problem.isTransient() && !isCoupled("porosity"))
    mooseWarning("LMPoroMaterial: running a transient simulation but did not supplied porosity!");
}

void
LMPoroMaterial::computeQpProperties()
{
  // Compressibilities
  Real Cf = (_Kf != 0.0) ? 1.0 / _Kf : 0.0;
  Real Cs = (_Ks != 0.0) ? 1.0 / _Ks : 0.0;
  ADReal Cd = (_coupled_mech && ((*_K)[_qp] != 0.0)) ? 1.0 / (*_K)[_qp] : 0.0;

  // Biot coefficient
  _biot[_qp] = 1.0;
  if (_coupled_mech && (Cd != 0.0))
    _biot[_qp] -= Cs / Cd;

  // Storage
  _C_biot[_qp] = _porosity[_qp] * Cf;
  if (_coupled_mech)
    _C_biot[_qp] += (_biot[_qp] - _porosity[_qp]) * Cs;

  // Fluid mobility
  _fluid_mob[_qp] = _perm / _fluid_visco;

  // Poro-mechanics
  _poro_mech[_qp] = 0.0;
  if (_coupled_mech && _fe_problem.isTransient())
  {
    _poro_mech[_qp] += _biot[_qp] * (*_strain_increment)[_qp].trace() / _dt;
    if (_has_ve)
      _poro_mech[_qp] += (1.0 - _biot[_qp]) * (*_viscous_strain_incr)[_qp].trace() / _dt;
    if (_has_vp)
      _poro_mech[_qp] += (1.0 - _biot[_qp]) * (*_plastic_strain_incr)[_qp].trace() / _dt;
  }
}