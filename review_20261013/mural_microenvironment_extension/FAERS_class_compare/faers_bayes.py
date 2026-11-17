"""Reproducible BCPNN and two-item Gamma-Poisson shrinker signal statistics.

The GPS here is the two-item signal component used in MGPS, not a higher-order
interaction analysis. Fit on the complete prespecified drug-event background.
BCPNN variants: PhViD 1.0.8 Bate/Gould normal approximation and Noren 2006
Dirichlet-multinomial Monte Carlo, both reported on the log2 IC scale.
GPS posterior formulas are checked against CRAN openEBGM 0.10.0 source.
No posterior hyperparameters are silently clamped.

CLI example:
  py analysis/faers_bayes.py fit background.csv gps_fit.json
  # CSV columns: N,E; optional weight (or weights or weight(s) renamed to W).
  py analysis/faers_bayes.py bcpnn 3 997 97 98903 bcpnn.json
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq, minimize
from scipy.special import digamma, expit, gammaln, logit, polygamma
from scipy.stats import gamma, norm

DEFAULT_STARTS = (
    (0.5, 0.5, 2.0, 2.0, 0.2),
    (1.0, 1.0, 3.0, 3.0, 0.2),
    (0.2, 0.2, 2.0, 0.2, 0.9),
    (1.0, 2.0, 2.0, 0.5, 0.8),
    (0.2, 0.06, 1.4, 1.8, 0.1),
)
TRUNCATIONS = ('openebgm', 'mixture', 'none')


def _theta(theta):
    values = np.asarray(theta, dtype=float)
    if values.shape != (5,) or not np.all(np.isfinite(values)):
        raise ValueError('theta must contain five finite values: alpha1,beta1,alpha2,beta2,P')
    if np.any(values[:4] <= 0) or not 0 < values[4] < 1:
        raise ValueError('Gamma shape/rate must be positive and 0 < P < 1')
    return values


def _arrays(observed, expected, weights=None, positive=False):
    observed = np.asarray(observed, dtype=float)
    expected = np.asarray(expected, dtype=float)
    if observed.ndim > 1 or expected.ndim > 1:
        raise ValueError('Observed and expected inputs must be scalars or one-dimensional arrays')
    n, e = np.broadcast_arrays(np.atleast_1d(observed),np.atleast_1d(expected))
    n, e = n.ravel(), e.ravel()
    w = np.ones_like(n) if weights is None else np.broadcast_to(np.asarray(weights, dtype=float), n.shape).copy()
    if not len(n) or not np.all(np.isfinite(n)) or not np.all(np.isfinite(e)):
        raise ValueError('Counts and expected counts must be nonempty finite arrays')
    if np.any(n < int(positive)) or np.any(n != np.floor(n)) or np.any(e <= 0):
        raise ValueError('Observed counts must be nonnegative integers (positive for truncated fits), and E > 0')
    if not np.all(np.isfinite(w)) or np.any(w <= 0):
        raise ValueError('All weights must be finite and positive')
    return n, e, w


def _log_nb(n, e, alpha, beta):
    # Equivalent to nbinom.logpmf(n, alpha, beta/(beta+E)); stable at tiny E.
    return (gammaln(n + alpha) - gammaln(alpha) - gammaln(n + 1)
            - alpha * np.log1p(e / beta) - n * np.log1p(beta / e))


def _component(n, e, alpha, beta):
    log_ratio = np.log1p(e / beta)
    lp = _log_nb(n, e, alpha, beta)
    lp0 = -alpha * log_ratio
    log_survival = np.log(-np.expm1(lp0))
    d_a = alpha * (digamma(n + alpha) - digamma(alpha) - log_ratio)
    d_b = (alpha * e - n * beta) / (beta + e)
    # p0/(1-p0); exp form avoids overflow of expm1(-lp0) at large -lp0.
    ratio = np.exp(lp0 - log_survival)
    ds_a = ratio * alpha * log_ratio
    ds_b = -ratio * alpha * e / (beta + e)
    return lp, log_survival, d_a, d_b, ds_a, ds_b


def _objective(z, n, e, w, truncation):
    """Average negative log likelihood and analytic transformed-parameter gradient."""
    theta = np.r_[np.exp(z[:4]), expit(z[4])]
    a1, b1, a2, b2, p = theta
    c1 = _component(n, e, a1, b1)
    c2 = _component(n, e, a2, b2)
    lp1, s1, da1, db1, dsa1, dsb1 = c1
    lp2, s2, da2, db2, dsa2, dsb2 = c2
    logp, log1p = np.log(p), np.log1p(-p)
    if truncation == 'openebgm':
        lp1, lp2 = lp1 - s1, lp2 - s2
        da1, db1 = da1 - dsa1, db1 - dsb1
        da2, db2 = da2 - dsa2, db2 - dsb2
    log_mix = np.logaddexp(logp + lp1, log1p + lp2)
    q = np.exp(logp + lp1 - log_mix)
    grad = [q * da1, q * db1, (1-q) * da2, (1-q) * db2, q-p]
    if truncation == 'mixture':
        log_surv_mix = np.logaddexp(logp+s1, log1p+s2)
        qs = np.exp(logp+s1-log_surv_mix)
        log_mix = log_mix - log_surv_mix
        grad[0] -= qs*dsa1
        grad[1] -= qs*dsb1
        grad[2] -= (1-qs)*dsa2
        grad[3] -= (1-qs)*dsb2
        grad[4] -= qs-p
    scale = np.sum(w)
    value = -float(np.dot(w, log_mix))/scale
    gradient = -np.asarray([np.dot(w, g) for g in grad])/scale
    return value, gradient


def gps_nll(theta, observed, expected, weights=None, truncation='mixture'):
    """Weighted negative log likelihood, with explicitly selected truncation.

    openebgm: mix component-wise zero-truncated NB distributions (R source).
    mixture: condition the complete NB mixture on N>0.
    none: unconditional mixture; requires zero cells in the fitting universe.
    """
    if truncation not in TRUNCATIONS:
        raise ValueError(f'truncation must be one of {TRUNCATIONS}')
    th = _theta(theta)
    n, e, w = _arrays(observed, expected, weights, positive=truncation != 'none')
    z = np.r_[np.log(th[:4]), logit(th[4])]
    return _objective(z, n, e, w, truncation)[0]*float(w.sum())


def _canonical(theta):
    th = np.asarray(theta, dtype=float).copy()
    # Component labels are exchangeable; report lower prior mean first.
    if th[0]/th[1] > th[2]/th[3]:
        th = th[[2, 3, 0, 1, 4]]
        th[4] = 1-th[4]
    return th


def fit_gps(observed, expected, weights=None, starts=None, truncation='mixture',
            shape_rate_bounds=(1e-6, 1e6), weight_bounds=(1e-8, 1-1e-8),
            maxiter=1000, ftol=1e-11, gtol=1e-7):
    """Fit two Gamma prior components by multistart NB mixture likelihood.

    Bounds apply only to optimization and are returned in metadata. No data are
    randomly subsampled and no silent data squashing is performed. If weights
    describe squashed rows, the caller must record its squashing algorithm.
    Returns the best successful fit (or failed best candidate with success=False).
    """
    if truncation not in TRUNCATIONS:
        raise ValueError(f'truncation must be one of {TRUNCATIONS}')
    n, e, w = _arrays(observed, expected, weights, positive=truncation != 'none')
    lo, hi = map(float, shape_rate_bounds)
    plo, phi = map(float, weight_bounds)
    if not 0 < lo < hi or not 0 < plo < phi < 1:
        raise ValueError('Invalid shape/rate or mixture-weight optimizer bounds')
    initial = list(DEFAULT_STARTS if starts is None else starts)
    if not initial:
        raise ValueError('At least one start is required')
    bounds = [(np.log(lo), np.log(hi))]*4 + [(logit(plo), logit(phi))]
    attempts = []
    for start in initial:
        th = _theta(start)
        if np.any(th[:4] < lo) or np.any(th[:4] > hi) or not plo <= th[4] <= phi:
            raise ValueError('Start lies outside the specified optimizer bounds')
        z0 = np.r_[np.log(th[:4]), logit(th[4])]
        opt = minimize(_objective, z0, args=(n,e,w,truncation), jac=True,
                       method='L-BFGS-B', bounds=bounds,
                       options={'maxiter':maxiter, 'ftol':ftol, 'gtol':gtol, 'maxls':40})
        fitted = np.r_[np.exp(opt.x[:4]), expit(opt.x[4])]
        boundary = [int(i) for i,(x,(low,high)) in enumerate(zip(opt.x,bounds))
                    if min(x-low,high-x) < 1e-5]
        attempts.append({'start':th.tolist(), 'theta':_canonical(fitted).tolist(),
                         'success':bool(opt.success and np.isfinite(opt.fun)),
                         'status':int(opt.status), 'message':str(opt.message),
                         'negative_log_likelihood':float(opt.fun*w.sum()),
                         'mean_negative_log_likelihood':float(opt.fun),
                         'iterations':int(opt.nit), 'function_evaluations':int(opt.nfev),
                         'gradient_inf_norm':float(np.max(np.abs(opt.jac))),
                         'optimizer_bound_indices':boundary})
    converged = [a for a in attempts if a['success']]
    candidates = converged or attempts
    best = min(candidates, key=lambda a:a['negative_log_likelihood'])
    theta_best = np.asarray(best['theta'])
    close = sum(np.all(np.abs(np.log(np.asarray(a['theta'])[:4]/theta_best[:4])) <= np.log(1.1))
                and abs(a['theta'][4]-theta_best[4]) <= .02 for a in converged)
    warnings = []
    if truncation == 'openebgm':
        warnings.append('Compatibility objective mixes separately truncated components; it differs from conditioning the full prior mixture on N>0.')
    if truncation == 'none' and not np.any(n==0):
        warnings.append('Unconditional fit has no zero-count cells; verify this represents the complete fitting universe.')
    if not converged:
        warnings.append('No start converged: posterior results are not a validated empirical-Bayes fit.')
    if close < 2:
        warnings.append('Fewer than two converged starts agree within the stated parameter tolerance.')
    if any(not a['success'] and a['negative_log_likelihood'] < best['negative_log_likelihood']-1e-6 for a in attempts):
        warnings.append('A nonconverged start achieved a lower objective than the selected converged fit; increase maxiter or inspect that start.')
    if best['optimizer_bound_indices']:
        warnings.append('Best fit reached an optimizer bound; assess boundary sensitivity.')
    if np.unique(n).size == 1:
        warnings.append('All observed counts are identical; mixture parameters may be poorly identified.')
    return {'model':'two-item Gamma-Poisson shrinker; not higher-order MGPS interactions',
            'theta_names':['alpha1','beta1_rate','alpha2','beta2_rate','P'],
            'theta':best['theta'], 'success':bool(converged),
            'negative_log_likelihood':best['negative_log_likelihood'],
            'truncation':truncation, 'fit_rows':int(len(n)), 'sum_weights':float(w.sum()),
            'zero_count_rows':int(np.sum(n==0)), 'zero_count_weight':float(w[n==0].sum()),
            'likelihood_formula':{'mixture':'(P*f1+(1-P)*f2)/(P*S1+(1-P)*S2)',
                                  'openebgm':'P*f1/S1+(1-P)*f2/S2',
                                  'none':'P*f1+(1-P)*f2'}[truncation],
            'weighted':weights is not None, 'shape_rate_bounds':[lo,hi],
            'weight_bounds':[plo,phi], 'optimizer':'scipy L-BFGS-B, analytic gradient',
            'optimizer_settings':{'maxiter':maxiter,'ftol':ftol,'gtol':gtol,'maxls':40},
            'component_order':'ascending alpha/beta prior mean',
            'converged_starts':len(converged), 'parameter_agreement_count':int(close),
            'multistart_parameter_agreement_verified':bool(close>=2),
            'parameter_agreement_rule':'shape/rate within factor 1.1, absolute P difference <=0.02',
            'warnings':warnings, 'attempts':attempts,
            'software':{'numpy':np.__version__,'scipy':scipy.__version__}}


def gps_posterior(observed, expected, theta, quantiles=(.05,.95)):
    """Conditional Gamma-mixture posterior. Returns arrays, including for scalars.

    EB05/EB95 form a central 90% posterior credible interval. Posterior uses the
    supplied theta unchanged; its suitability is the responsibility of fit audit.
    """
    th = _theta(theta)
    n,e,_ = _arrays(observed,expected)
    a1,b1,a2,b2,p = th
    lf1,lf2 = _log_nb(n,e,a1,b1),_log_nb(n,e,a2,b2)
    logden = np.logaddexp(np.log(p)+lf1,np.log1p(-p)+lf2)
    q = np.exp(np.log(p)+lf1-logden)
    logmean = q*(digamma(a1+n)-np.log(b1+e))+(1-q)*(digamma(a2+n)-np.log(b2+e))
    out = {'EBGM':np.exp(logmean), 'posterior_component1_weight':q}
    for probability in quantiles:
        probability = float(probability)
        if not 0 < probability < 1:
            raise ValueError('Quantile probabilities must be between zero and one')
        answers = []
        for count,expct,weight in zip(n,e,q):
            shapes,scales = (a1+count,a2+count),(1/(b1+expct),1/(b2+expct))
            # The mixture quantile lies between the two component quantiles.
            ends = gamma.ppf(probability, shapes, scale=scales)
            left,right = float(min(ends)),float(max(ends))
            if weight == 1:
                answer = float(ends[0])
            elif weight == 0:
                answer = float(ends[1])
            elif left == right:
                answer = left
            else:
                cdf = lambda x: weight*gamma.cdf(x,shapes[0],scale=scales[0])+(1-weight)*gamma.cdf(x,shapes[1],scale=scales[1])-probability
                fl,fr = cdf(left),cdf(right)
                # Quantile endpoints can miss zero by rounding when weights nearly saturate.
                if abs(fl) <= 1e-14:
                    answer = left
                elif abs(fr) <= 1e-14:
                    answer = right
                else:
                    if fl > 0:
                        left = 0.
                    while fr < 0:
                        right *= 2
                        if not np.isfinite(right):
                            raise FloatingPointError('Could not bracket posterior quantile')
                        fr = cdf(right)
                    # Solve on log(lambda), preserving relative accuracy in very small tails.
                    smallest = np.nextafter(0.,1.)
                    lower = max(left,smallest)
                    if cdf(lower) >= 0:
                        answer = 0. if left==0 else lower
                    else:
                        root = brentq(lambda logx:cdf(np.exp(logx)),np.log(lower),np.log(right),
                                      xtol=1e-12,rtol=1e-12)
                        answer = float(np.exp(root))
            answers.append(answer)
        label = f'EB{int(round(probability*100)):02d}' if abs(probability*100-round(probability*100))<1e-10 else f'posterior_quantile_{probability:g}'
        out[label] = np.asarray(answers)
        out[label+'_precision_limited'] = out[label] < np.finfo(float).tiny
    return out


def bcpnn(a, b, c, d, mc_samples=1000000, seed=20261005, batch_size=100000):
    """Both PhViD BCPNN variants for one report-level 2x2 table.

    a=drug+event, b=drug/no event, c=event/no drug, d=neither. Legacy IC025 is
    a normal approximation. Dirichlet IC025 is a Monte Carlo order statistic
    on log2 scale, matching PhViD's rank round(M*.025), converted from ln.
    """
    cells = np.asarray([a,b,c,d],dtype=float)
    if not np.all(np.isfinite(cells)) or np.any(cells<0) or np.any(cells!=np.floor(cells)):
        raise ValueError('a,b,c,d must be finite nonnegative integer report counts')
    if cells.sum()<=0:
        raise ValueError('The 2x2 table must contain at least one report')
    if int(mc_samples)!=mc_samples or mc_samples<1000 or int(batch_size)!=batch_size or batch_size<1:
        raise ValueError('mc_samples must be an integer >=1000 and batch_size a positive integer')
    a,b,c,d = cells
    total,drug,event = cells.sum(),a+b,a+c
    p1,p2,q1,q2 = 1+drug,1+total-drug,1+event,1+total-event
    r1,r2 = 1+a,total-a-1+(total+2)**2/(p1*q1)
    ln2 = np.log(2.)
    mean = (digamma(r1)-digamma(r1+r2)-digamma(p1)+digamma(p1+p2)-digamma(q1)+digamma(q1+q2))/ln2
    variance = (polygamma(1,r1)-polygamma(1,r1+r2)+polygamma(1,p1)-polygamma(1,p1+p2)+polygamma(1,q1)-polygamma(1,q1+q2))/ln2**2
    u,v = (drug+.5)/(total+1),(event+.5)/(total+1)
    prior_total = .5/(u*v)
    prior = prior_total*np.array([u*v,u*(1-v),(1-u)*v,(1-u)*(1-v)])
    posterior = cells+prior
    rng = np.random.default_rng(seed)
    values = np.empty(int(mc_samples))
    batch_quantiles = []
    for start in range(0,int(mc_samples),int(batch_size)):
        stop = min(start+int(batch_size),int(mc_samples))
        draws = rng.gamma(posterior, size=(stop-start,4))
        # IC using unnormalised Gamma draws: avoids loss from normalizing twice.
        ic = (np.log(draws[:,0])+np.log(draws.sum(axis=1))
              -np.log(draws[:,0]+draws[:,1])-np.log(draws[:,0]+draws[:,2]))/ln2
        if not np.all(np.isfinite(ic)):
            raise FloatingPointError('Non-finite Dirichlet IC draw; inspect posterior shape values')
        values[start:stop] = ic
        batch_quantiles.append(float(np.quantile(ic,.025)))
    rank = max(1,int(np.rint(mc_samples*.025)))
    lower = float(np.partition(values,rank-1)[rank-1])
    rankhi = max(1,int(np.rint(mc_samples*.975)))
    upper = float(np.partition(values,rankhi-1)[rankhi-1])
    # Approximate MC standard error from equal-sized independent batch quantiles.
    equal_batches = mc_samples % batch_size == 0 and len(batch_quantiles)>=2
    mcse = float(np.std(batch_quantiles,ddof=1)/np.sqrt(len(batch_quantiles))) if equal_batches else None
    expected = drug*event/total
    warnings = []
    if mcse is not None and abs(lower) <= 3*mcse:
        warnings.append('IC025 is within three estimated Monte Carlo standard errors of zero; increase draws or replicate before threshold classification.')
    return {'a':int(a),'b':int(b),'c':int(c),'d':int(d),'expected':float(expected),
            'IC_observed_log2_OE':float(np.log2(a/expected)) if a>0 and expected>0 else None,
            'IC_mean_legacy_normal':float(mean), 'IC_variance_legacy':float(variance),
            'IC025_legacy_normal':float(mean+norm.ppf(.025)*np.sqrt(variance)),
            'IC975_legacy_normal':float(mean+norm.ppf(.975)*np.sqrt(variance)),
            'IC_mean_dirichlet_MC':float(np.mean(values)),
            'IC025_dirichlet_MC':lower,'IC975_dirichlet_MC':upper,
            'IC025_dirichlet_MCSE_approx':mcse,
            'posterior_probability_IC_le_zero':float(np.mean(values<=0)),
            'dirichlet_prior':prior.tolist(),'dirichlet_posterior':posterior.tolist(),
            'mc_samples':int(mc_samples),'mc_seed':int(seed),'mc_batch_size':int(batch_size),
            'mc_quantile_rank':rank,'mc_quantile_rule':'PhViD order statistic round(M*0.025), one-based',
            'IC_scale':'log2','legacy_method':'Bate beta-binomial / Gould exact moments / normal quantile approximation',
            'MC_method':'Noren 2006 Dirichlet-multinomial; PhViD 1.0.8 formula with log2 correction',
            'warnings':warnings,
            'software':{'numpy':np.__version__,'scipy':scipy.__version__}}


def _read_background(path):
    with Path(path).open(newline='',encoding='utf-8-sig') as stream:
        rows = list(csv.DictReader(stream))
    if not rows or 'N' not in rows[0] or 'E' not in rows[0]:
        raise ValueError('Background CSV requires columns N and E; optional weight, weights, or W')
    weight_key = next((key for key in ('weight','weights','W') if key in rows[0]),None)
    return ([float(row['N']) for row in rows],[float(row['E']) for row in rows],
            [float(row[weight_key]) for row in rows] if weight_key else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action',required=True)
    fit = sub.add_parser('fit'); fit.add_argument('background'); fit.add_argument('output')
    fit.add_argument('--truncation',choices=TRUNCATIONS,default='mixture')
    fit.add_argument('--maxiter',type=int,default=1000)
    bc = sub.add_parser('bcpnn')
    for name in ('a','b','c','d'):bc.add_argument(name,type=int)
    bc.add_argument('output');bc.add_argument('--mc-samples',type=int,default=1000000)
    bc.add_argument('--seed',type=int,default=20261005)
    args = parser.parse_args()
    if args.action=='fit':
        n,e,w = _read_background(args.background)
        result = fit_gps(n,e,w,truncation=args.truncation,maxiter=args.maxiter)
    else:
        result = bcpnn(args.a,args.b,args.c,args.d,mc_samples=args.mc_samples,seed=args.seed)
    Path(args.output).write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    print(args.output)
    if args.action=='fit' and not result['success']:
        raise SystemExit('No GPS fit converged; see output diagnostics')


if __name__=='__main__':
    main()
