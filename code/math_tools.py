
import numpy as np
from stat_tools import proportion_ratio_CI, proportion_sum_CI

def fitness_effect (pN,
                    pM,
                    pS,
                    k_N,
                    n_N,
                    k_M,
                    n_M,
                    pT_S = None,
                    edgotype = 'edgetic',
                    CI = 95,
                    output = True):
    
    # abbreviations for different structural regions
    edgotype_symbol = {'quasi-null':'Q', 'edgetic':'E', 'quasi-wild-type':'W'}
    
    # abbreviation for selected structural region
    T = edgotype_symbol [edgotype]
    
    # prior conditional probabilities
    cond_probs = ['P(%s|%s)' % (T, p) for p in ('N','M','S')]
    
    # posterior conditional probabilities
    posteriors = ['P(%s|%s)' % (p, T) for p in ('N','M','S')]
    
    # Probability for effectively neutral mutations (N) to be of selected edgotype (T)
    pT_N = k_N / n_N
    
    # Probability for mildly deleterious mutations (M) to be of selected edgotype (T)
    pT_M = k_M / n_M
    
    # Assume strongly detrimental mutations are similar to mildly deleterious mutations
    if pT_S is None:
        pT_S = pT_M
    
    # Probability for a new missense mutation to be of selected edgotype
    pT = (pT_N * pN) + (pT_M * pM) + (pT_S * pS)
    
    if pT == 0:
        return {}
    else:
        # Probability for mutations of selected edgotype to be effectively neutral
        pN_T = pT_N * pN / pT
    
        # Probability for mutations of selected edgotype to be mildly deleterious
        pM_T = pT_M * pM / pT
    
        # Probability for mutations of selected edgotype to be strongly detrimental
        pS_T = pT_S * pS / pT
    
        allresults = {prob:p for prob, p in zip(posteriors, (pN_T, pM_T, pS_T))}
        allresults['P(%s)' % T] = pT
    
        if output:
            print()
            print('Fitness effect calculation for %s (%s) mutations:' % (edgotype, T))
            print('P(N) = %.1f %%' % (100 * pN))
            print('P(M) = %.1f %%' % (100 * pM))
            print('P(S) = %.1f %%' % (100 * pS))
            print('%s = %.1f %%' % (cond_probs[0], 100 * pT_N))
            print('%s = %.1f %%' % (cond_probs[1], 100 * pT_M))
            print('%s = %.1f %%' % (cond_probs[2], 100 * pT_S))
            print('P(%s) = %sP(N) + %sP(M) + %sP(S) = %.1f %%' 
                    % (T, cond_probs[0], cond_probs[1], cond_probs[2], 100 * pT))
            print('Probability for %s mutations to be effectively neutral %s = %.1f %%' 
                    % (edgotype, posteriors[0], 100 * pN_T))
            print('Probability for %s mutations to be mildly deleterious %s = %.1f %%' 
                    % (edgotype, posteriors[1], 100 * pM_T))
            print('Probability for %s mutations to be strongly detrimental %s = %.1f %%' 
                    % (edgotype, posteriors[2], 100 * pS_T))
    
        # calculate 95% confidence interval
        if CI:
            for prob, p in zip(posteriors, (pN_T, pM_T, pS_T)):
                if prob is posteriors[2] and p == 0:
                    p_lower, p_upper = 0, 0
                else:
                    if prob is posteriors[0]:
                        pr_lower, pr_upper = proportion_ratio_CI (k_M,
                                                                  n_M,
                                                                  k_N,
                                                                  n_N,
                                                                  a = pM / pN,
                                                                  b = pT_S * pS / pN,
                                                                  conf = CI)
                    elif prob is posteriors[1]:
                        pr_lower, pr_upper = proportion_ratio_CI (k_N,
                                                                  n_N,
                                                                  k_M,
                                                                  n_M,
                                                                  a = pN / pM,
                                                                  b = pT_S * pS / pM,
                                                                  conf = CI)
                    elif prob is posteriors[2]:
                        pr_lower, pr_upper = proportion_sum_CI (k_N,
                                                                n_N,
                                                                k_M,
                                                                n_M,
                                                                a = pN / (pT_S * pS),
                                                                b = pM / (pT_S * pS),
                                                                conf = CI)
                    p_lower = 1 / (1 + pr_upper)
                    p_upper = 1 / (1 + pr_lower)
                if output:
                    print( '%.1f%% confidence interval for %s = (%f, %f)' 
                            % (CI, prob, 100 * p_lower, 100 * p_upper) )
                if (not np.isnan(p_lower)) and (not np.isnan(p_upper)):
                    allresults[prob + '_CI'] = [p - p_lower, p_upper - p]
        return allresults
