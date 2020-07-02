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

#include "LMSingleVarUpdate.h"
#include "ElasticityTensorTools.h"

InputParameters
LMSingleVarUpdate::validParams()
{
  InputParameters params = LMViscoPlasticUpdate::validParams();
  params.addClassDescription("Base class for a single variable viscoplastic update.");
  return params;
}

LMSingleVarUpdate::LMSingleVarUpdate(const InputParameters & parameters)
  : LMViscoPlasticUpdate(parameters)
{
}

void
LMSingleVarUpdate::viscoPlasticUpdate(ADRankTwoTensor & stress,
                                      const ADRankFourTensor & Cijkl,
                                      ADRankTwoTensor & elastic_strain_incr)
{
  // Here we do an iterative update with a single variable (usually scalar viscoplastic strain rate)
  // We are trying to find the zero of the function F which is defined as:
  // F(gamma_vp) = yield - eta * gamma_vp^(1 / n)
  // gamma_vp: scalar viscoplastic strain rate (scalar)
  // yield: the yield function
  // eta: the viscoplastic viscosity
  // n: exponent for Perzyna-like flow rule
  // flow rule: gamma_vp = (yield / eta)^n

  // Trial stress
  _stress_tr = stress;
  // Elastic moduli
  _K = ElasticityTensorTools::getIsotropicBulkModulus(Cijkl);
  _G = ElasticityTensorTools::getIsotropicShearModulus(Cijkl);

  // Initialize plastic strain increment
  _plastic_strain_incr[_qp].zero();

  // Pre return map calculations (model specific)
  preReturnMap();

  // Check yield function
  _yield_function[_qp] = yieldFunction(0.0);
  if (_yield_function[_qp] <= _abs_tol) // Elastic
    return;

  // Viscoplastic update
  ADReal gamma_vp = returnMap();

  // Update quantities
  _yield_function[_qp] = yieldFunction(gamma_vp);
  _plastic_strain_incr[_qp] = reformPlasticStrainTensor(gamma_vp);
  elastic_strain_incr -= _plastic_strain_incr[_qp];
  stress -= Cijkl * _plastic_strain_incr[_qp];
  postReturnMap(gamma_vp);
}

ADReal
LMSingleVarUpdate::returnMap()
{
  // Initialize scalar viscoplastic strain rate
  ADReal gamma_vp = 0.0;

  // Initial residual
  ADReal res_ini = residual(gamma_vp);

  ADReal res = res_ini;
  ADReal jac = jacobian(gamma_vp);

  // Newton loop
  for (unsigned int iter = 0; iter < _max_its; ++iter)
  {
    gamma_vp -= res / jac;

    res = residual(gamma_vp);
    jac = jacobian(gamma_vp);

    // Convergence check
    if ((std::abs(res) <= _abs_tol) || (std::abs(res / res_ini) <= _rel_tol))
      return gamma_vp;
  }
  throw MooseException("LMSingleVarUpdate: maximum number of iterations exceeded in 'returnMap'!");
}

ADReal
LMSingleVarUpdate::residual(const ADReal & gamma_vp)
{
  ADReal res = yieldFunction(gamma_vp);
  if (gamma_vp != 0.0)
    res -= _eta_p * std::pow(gamma_vp, 1.0 / _n);

  return res;
}

ADReal
LMSingleVarUpdate::jacobian(const ADReal & gamma_vp)
{
  ADReal jac = yieldFunctionDeriv(gamma_vp);
  if (gamma_vp != 0.0)
    jac -= _eta_p / _n * std::pow(gamma_vp, 1.0 / _n - 1.0);

  return jac;
}