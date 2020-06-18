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

#include "LMViscoElasticUpdate.h"

InputParameters
LMViscoElasticUpdate::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("Base class for the viscoelastic correction.");
  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");
  params.addRangeCheckedParam<Real>("abs_tolerance",
                                    1.0e-10,
                                    "abs_tolerance > 0.0",
                                    "The absolute tolerance for the iterative update.");
  params.addRangeCheckedParam<Real>("rel_tolerance",
                                    1.0e-10,
                                    "rel_tolerance > 0.0",
                                    "The relative tolerance for the iterative update.");
  params.addRangeCheckedParam<unsigned int>(
      "max_iterations",
      200,
      "max_iterations >= 1",
      "The maximum number of iterations for the iterative update");
  return params;
}

LMViscoElasticUpdate::LMViscoElasticUpdate(const InputParameters & parameters)
  : ADMaterial(parameters),
    _abs_tol(getParam<Real>("abs_tolerance")),
    _rel_tol(getParam<Real>("rel_tolerance")),
    _max_its(getParam<unsigned int>("max_iterations")),
    _viscosity(declareADProperty<Real>("effective_viscosity")),
    _viscous_strain_incr(declareADProperty<RankTwoTensor>("viscous_strain_increment"))
{
}

void
LMViscoElasticUpdate::setQp(unsigned int qp)
{
  _qp = qp;
}

void
LMViscoElasticUpdate::viscoElasticUpdate(ADRankTwoTensor & stress,
                                         const RankFourTensor & Cijkl,
                                         ADRankTwoTensor & elastic_strain_incr)
{
  // Here we do an iterative update with a single variable (usually scalar viscous strain rate)
  // We are trying to find the zero of the function F which is defined as:
  // F(gamma_v) = yield - eta * gamma_vp^(1 / n)
  // gamma_v: scalar viscous strain rate (scalar)
  // yield: the yield function
  // eta: the viscoplastic viscosity
  // flow rule: gamma_v = (yield / eta)^n

  // Trial stress
  _stress_tr = stress;
  _tau_tr = std::sqrt(0.5) * _stress_tr.deviatoric().L2norm();

  // Elastic moduli
  _G = ElasticityTensorTools::getIsotropicShearModulus(Cijkl);

  // Initialize plastic strain increment
  _viscous_strain_incr[_qp].zero();

  if (MooseUtils::absoluteFuzzyEqual(_stress_tr.deviatoric().L2norm(), 0.0))
  {
    _viscosity[_qp] = effectiveViscosity(0.0);
    return;
  }

  // Pre return map calculations (model specific)
  preReturnMap();

  // Viscoplastic update
  ADReal gamma_v = returnMap();

  // Update quantities
  _viscosity[_qp] = effectiveViscosity(gamma_v);
  _viscous_strain_incr[_qp] = reformViscousStrainTensor(gamma_v);
  elastic_strain_incr -= _viscous_strain_incr[_qp];
  stress -= Cijkl * _viscous_strain_incr[_qp];
  postReturnMap(gamma_v);
}

ADReal
LMViscoElasticUpdate::returnMap()
{
  // Initialize scalar viscous strain rate
  ADReal gamma_v = 0.0;

  // Initial residual
  ADReal res_ini = residual(gamma_v);

  ADReal res = res_ini;
  ADReal jac = jacobian(gamma_v);

  // Newton loop
  for (unsigned int iter = 0; iter < _max_its; ++iter)
  {
    gamma_v -= res / jac;

    res = residual(gamma_v);
    jac = jacobian(gamma_v);

    // Convergence check
    if ((std::abs(res) <= _abs_tol) || (std::abs(res / res_ini) <= _rel_tol))
      return gamma_v;
  }
  throw MooseException(
      "LMViscoElasticUpdate: maximum number of iterations exceeded in 'returnMap'!");
}

ADReal
LMViscoElasticUpdate::residual(const ADReal & gamma_v)
{
  ADReal creep_rate = creepRate(gamma_v);
  ADReal tau = stressInvariant(gamma_v);

  return _tau_tr - tau - 2.0 * _G * creep_rate * _dt;
}

ADReal
LMViscoElasticUpdate::jacobian(const ADReal & gamma_v)
{
  ADReal dcreep_rate = creepRateDeriv(gamma_v);
  ADReal dtau = stressInvariantDeriv(gamma_v);

  return -dtau - 2.0 * _G * dcreep_rate * _dt;
}

ADReal
LMViscoElasticUpdate::stressInvariant(const ADReal & gamma_v)
{
  return _tau_tr - 2.0 * _G * gamma_v * _dt;
}

ADReal
LMViscoElasticUpdate::stressInvariantDeriv(const ADReal & /*gamma_v*/)
{
  return -2.0 * _G * _dt;
}

ADRankTwoTensor
LMViscoElasticUpdate::reformViscousStrainTensor(const ADReal & gamma_v)
{
  ADRankTwoTensor flow_dir =
      (_tau_tr != 0.0) ? _stress_tr.deviatoric() / _tau_tr : ADRankTwoTensor();

  return gamma_v * _dt * flow_dir;
}