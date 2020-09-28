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

#include "LMTwoVarUpdate.h"
#include "ElasticityTensorTools.h"

InputParameters
LMTwoVarUpdate::validParams()
{
  InputParameters params = LMViscoPlasticUpdate::validParams();
  params.addClassDescription("Base class for a single variable viscoplastic update.");
  params.addCoupledVar("fluid_pressure", 0, "The fluid pressure variable.");
  params.addRangeCheckedParam<Real>("reference_fluid_pressure",
                                    0.0,
                                    "reference_fluid_pressure>=0.0",
                                    "The reference fluid pressure.");
  params.addRangeCheckedParam<Real>(
      "Arrhenius_coefficient",
      0.0,
      "Arrhenius_coefficient>=0",
      "The Arrhenius coefficient for the fludi pressure activated viscosity.");
  return params;
}

LMTwoVarUpdate::LMTwoVarUpdate(const InputParameters & parameters)
  : LMViscoPlasticUpdate(parameters),
    _pf(adCoupledValue("fluid_pressure")),
    _pf0(getParam<Real>("reference_fluid_pressure")),
    _Ar(getParam<Real>("Arrhenius_coefficient"))
{
}

void
LMTwoVarUpdate::viscoPlasticUpdate(ADRankTwoTensor & stress,
                                   const ADRankFourTensor & Cijkl,
                                   ADRankTwoTensor & elastic_strain_incr)
{
  // Here we do an iterative update with two variables (usually scalar volumetric and deviatoric
  // viscoplastic strain rates)
  // We are trying to find the zero of the two residual functions Rv, Rd which are defined as:
  // Rv(gamma_v, gamma_d) = p - p0 - eta * gamma_v^(1 / n)
  // Rd(gamma_v, gamma_d) = q - q0 - eta * gamma_d^(1 / n)
  // gamma_v: scalar volumetric viscoplastic strain rate
  // gamma_d: scalar deviatoric viscoplastic strain rate
  // p0, q0: the projection of the pressure and deviatoric stress on the yield function
  // eta: the viscoplastic viscosity
  // n: exponent for Perzyna-like flow rule

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
  ADReal chi_v = 0.0, chi_d = 0.0;
  updateDissipativeStress(0.0, 0.0, chi_v, chi_d);
  _yield_function[_qp] = yieldFunction(chi_v, chi_d);
  if (_yield_function[_qp] <= _abs_tol) // Elastic
    return;

  // Viscoplastic update
  ADReal gamma_v = 0.0, gamma_d = 0.0;
  returnMap(gamma_v, gamma_d);

  // Update quantities
  updateDissipativeStress(gamma_v, gamma_d, chi_v, chi_d);
  _yield_function[_qp] = yieldFunction(chi_v, chi_d);
  _plastic_strain_incr[_qp] = reformPlasticStrainTensor(gamma_v, gamma_d);
  elastic_strain_incr -= _plastic_strain_incr[_qp];
  stress -= Cijkl * _plastic_strain_incr[_qp];
  postReturnMap(gamma_v, gamma_d);
}

void
LMTwoVarUpdate::returnMap(ADReal & gamma_v, ADReal & gamma_d)
{
  // Initial residual
  ADReal resv_ini = 0.0, resd_ini = 0.0;
  residual(0.0, 0.0, resv_ini, resd_ini);
  ADReal res_ini = std::sqrt(Utility::pow<2>(resv_ini) + Utility::pow<2>(resd_ini));
  ADReal resv = resv_ini, resd = resd_ini;
  ADReal res = res_ini;

  // Initial jacobian
  ADReal jacvv = 0.0, jacdd = 0.0, jacvd = 0.0, jacdv = 0.0;
  jacobian(0.0, 0.0, jacvv, jacdd, jacvd, jacdv);

  // Useful stuff
  ADReal jac_full = jacvv * jacdd - jacvd * jacdv;
  ADReal resv_full = jacdd * resv - jacvd * resd;
  ADReal resd_full = jacvv * resd - jacdv * resv;

  // Newton loop
  for (unsigned int iter = 0; iter < _max_its; ++iter)
  {
    gamma_v -= resv_full / jac_full;
    gamma_d -= resd_full / jac_full;

    residual(gamma_v, gamma_d, resv, resd);
    jacobian(gamma_v, gamma_d, jacvv, jacdd, jacvd, jacdv);
    jac_full = jacvv * jacdd - jacvd * jacdv;
    resv_full = jacdd * resv - jacvd * resd;
    resd_full = jacvv * resd - jacdv * resv;
    res = std::sqrt(Utility::pow<2>(resv) + Utility::pow<2>(resd));

    // Convergence check
    if ((std::abs(res) <= _abs_tol) || (std::abs(res / res_ini) <= _rel_tol))
      return;
  }
  throw MooseException(
      "LMTwoVarUpdate: maximum number of iterations exceeded in 'returnMap'!\nInitial residual: ",
      res_ini.value(),
      "\nResidual: ",
      res.value(),
      "\n Vol strain rate: ",
      gamma_v.value(),
      "\n Dev strain Rate: ",
      gamma_d.value(),
      "\n");
}

void
LMTwoVarUpdate::residual(const ADReal & gamma_v,
                         const ADReal & gamma_d,
                         ADReal & resv,
                         ADReal & resd)
{
  overStress(gamma_v, gamma_d, resv, resd);
  resv -= std::pow(_eta_p, _n) * gamma_v * std::exp(_Ar * (_pf[_qp] - _pf0));
  resd -= std::pow(_eta_p, _n) * gamma_d * std::exp(_Ar * (_pf[_qp] - _pf0));
}

void
LMTwoVarUpdate::jacobian(const ADReal & gamma_v,
                         const ADReal & gamma_d,
                         ADReal & jacvv,
                         ADReal & jacdd,
                         ADReal & jacvd,
                         ADReal & jacdv)
{
  overStressDerivV(gamma_v, gamma_d, jacvv, jacdv);
  overStressDerivD(gamma_v, gamma_d, jacvd, jacdd);
  jacvv -= std::pow(_eta_p, _n) * std::exp(_Ar * (_pf[_qp] - _pf0));
  jacdd -= std::pow(_eta_p, _n) * std::exp(_Ar * (_pf[_qp] - _pf0));
}