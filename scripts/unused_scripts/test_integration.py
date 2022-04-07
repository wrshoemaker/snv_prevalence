import math
from scipy.integrate import quad, quadrature
from scipy.special import gamma
import scipy.stats

def integrand(f, D, theta, f_max):
    return ((1-f)**D) * math.exp(-1*f/f_max) * (f**(theta-1))


theta = 0.0562175907497403
f_max = 0.025
D = 21


min_f_gamma = 0.1
max_f_gamma = 1.9

I = quad(integrand, min_f_gamma, max_f_gamma, args=(D, theta, f_max))


truncated_correction = scipy.stats.gamma.cdf(max_f_gamma, theta, scale=f_max) - scipy.stats.gamma.cdf(min_f_gamma, theta, scale=f_max)

prob_absence = I / ( gamma(theta) * (f_max**theta) * truncated_correction )


print(I[0])

prob_absence_anal = (1 + D*f_max)**(-1*theta)


#print(prob_absence)



#### test integrand that can be solved
def integrand_full_gamma(f, D, theta, f_max):
    return math.exp(-1*f * (D+(f_max**-1) ) ) * (f**(theta-1))


I_full = quad(integrand_full_gamma, 0, 1, args=(D, theta, f_max))

prob_absence_full = I_full / ( gamma(theta) * (f_max**theta))

print(prob_absence[0], prob_absence_full[0])
