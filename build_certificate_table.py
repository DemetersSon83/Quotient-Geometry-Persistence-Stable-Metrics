"""Generate conservatively rounded printed certificates from exact records."""
from fractions import Fraction as F
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent

def rounded(x, digits=4, upper=False):
    x = F(x)
    assert x > 0
    exponent = len(str(x.numerator)) - len(str(x.denominator))
    power = F(10)**exponent
    while power > x:
        exponent -= 1
        power /= 10
    while 10*power <= x:
        exponent += 1
        power *= 10
    unit = power / 10**(digits-1)
    q = x / unit
    integer = -(-q.numerator // q.denominator) if upper else q.numerator // q.denominator
    value = integer*unit
    assert (value >= x) if upper else (value <= x)
    if integer == 10**digits:
        integer //= 10
        exponent += 1
    # The display mantissa is constructed as a string from an exact integer.
    mantissa = str(integer)
    mantissa = mantissa[:1] + '.' + mantissa[1:]
    return value, '$'+mantissa+r'\times10^{'+str(exponent)+'}$'

data = json.loads((HERE/'certificates/local_certificates.json').read_text())
rows, checks = [], []
for r in data['results']:
    eta, eta_tex = rounded(r['eta'])
    delay, delay_tex = rounded(r['tau'])
    theta = F(r['total_error'])
    c0, c1 = F(r['left_block_norm']), F(r['right_block_norm'])
    inverse = 2/(1-theta)*(c0+2*c1/delay)
    k, k_tex = rounded(inverse, 3, upper=True)
    label = r'$n='+str(r['n'])+'$'
    if 'near_centered' in r['label']:
        label = r'Near-centered $n=6$'
    theta_upper = F(-(-theta.numerator*1000//theta.denominator), 1000)
    rows.append(label+' & '+eta_tex+' & '+delay_tex+' & '+f'${float(theta_upper):.3f}$'+' & '+k_tex+r' \\')
    checks.append({'label':r['label'],'eta_down':str(eta),'tau_down':str(delay),
        'theta_up':str(theta_upper),'K_at_printed_tau_up':str(k),
        'exact_K_bound_at_printed_tau':str(inverse),
        'all_rounding_inequalities_verified_exactly':True})
text = r'''\begin{table}[htbp]
\centering\small
\setlength{\tabcolsep}{5pt}
\begin{tabular}{lrrrr}
\toprule
Reference & $\eta_{\rm print}$ & $\tau_{\rm print}$ & $\theta$ upper & $K_{\tau_{\rm print}}$ upper \\
\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule
\end{tabular}
\caption{Validated local recovery for every formation in the ball of radius
$\eta_{\rm print}$ and every $0<|h|\le\tau_{\rm print}$. Radii (radians) and
time endpoints (consensus units) are rounded downward; derivative errors and
inverse bounds are rounded upward. The last column is evaluated at the
printed endpoint. For smaller delays use formula~\eqref{loc:eq:inverse},
not the endpoint constant. Exact records retain the larger unrounded
neighborhoods and delays.}
\label{tab:certificates}
\end{table}
'''
(HERE/'certificate_table.tex').write_text(text)
(HERE/'certificates/printed_bounds_verification.json').write_text(json.dumps(checks,indent=2)+'\n')
print('Five printed certificates rounded in the guaranteed direction.')
