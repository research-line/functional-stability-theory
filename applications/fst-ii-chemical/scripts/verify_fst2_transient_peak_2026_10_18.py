"""Read-only verifier: prints JSON; never writes project files."""
import hashlib, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path
import sympy as sp
import mpmath as mp

BASE = Path(__file__).resolve().parents[1]
PKG = BASE / "_proof-notes/FST2_TRANSIENT_PEAK_2026-10-18"
spec = json.loads((PKG / "CASE_SPEC.json").read_text(encoding="utf-8"))
mp.mp.dps = spec["dps"]
tol = mp.mpf(spec["relative_tolerance"])
exact = []
def check(name, condition):
    if not condition:
        raise AssertionError(name)
    exact.append(name)
def zmat(m):
    return all(sp.simplify(v) == 0 for v in m)
R = sp.Rational
A, B, C = sp.symbols("A B C")
field = sp.Matrix([C/2-A*B, A*B-R(4,5)*B*C, R(4,5)*B*C-C/2])
xstar = sp.Matrix([R(1,6),R(5,8),R(5,24)])
Q = sp.Matrix([[1/sp.sqrt(2),1/sp.sqrt(6)],[-1/sp.sqrt(2),1/sp.sqrt(6)],[0,-2/sp.sqrt(6)]])
Jfull = field.jacobian([A,B,C]).subs(dict(zip([A,B,C],xstar)))
J = sp.simplify(Q.T*Jfull*Q)
Jtarget = sp.Matrix([[-R(13,24),-R(41,24)/sp.sqrt(3)],[1/(4*sp.sqrt(3)),-R(1,12)]])
a,b,D = R(5,16),sp.sqrt(23)/16,R(1588,621)
K = J + a*sp.eye(2)
check("equilibrium_from_polynomial_field", zmat(field.subs(dict(zip([A,B,C],xstar)))))
check("Q_orthonormal", zmat(Q.T*Q-sp.eye(2)))
check("Q_tangent", zmat(sp.ones(1,3)*Q))
check("J_from_polynomial_field", zmat(J-Jtarget))
check("tangent_invariant", zmat(Jfull*Q-Q*J))
check("trace_and_determinant", J.trace()==-2*a and J.det()==a*a+b*b)
check("K_trace_zero", sp.simplify(K.trace())==0)
check("K_squared", zmat(K*K+b*b*sp.eye(2)))
Dcalc = sp.simplify((sp.trace(K.T*K)/b**2-2)/4)
check("D_exact", Dcalc==D)
s,c = sp.symbols("s c", real=True)
M = c*sp.eye(2)+s*K/b
check("det_identity", sp.simplify(M.det()-(c*c+s*s))==0)
check("frobenius_identity", sp.simplify(sp.trace(M.T*M)-(2*(c*c+s*s)+4*D*s*s))==0)
r = R(20999,76224)
check("peak_r_exact", sp.simplify((b*b*D-a*a)/(D*(a*a+b*b))-r)==0)
check("peak_r_in_open_unit_interval", bool(0<r<1))
check("initial_growth_positive", bool(b*b*D>a*a))
check("peak_stationary_squared", sp.simplify(b*b*D*(1-r)-a*a*(1+D*r))==0)
u,v,w = sp.symbols("u v w",real=True)
KK=sp.Matrix([[u,v],[w,-u]])
check("D_nonnegative_identity", sp.expand(sp.trace(KK.T*KK)-2*KK.det()-4*u*u-(v+w)**2)==0)
x=sp.symbols("x",real=True)
der=sp.diff(sp.cos(x)/sp.sqrt(1+D*sp.sin(x)**2),x)
check("derivative_strict_monotonicity_identity", sp.trigsimp(der+(1+D)*sp.sin(x)/(1+D*sp.sin(x)**2)**R(3,2))==0)

def mpv(z): return mp.mpf(str(sp.N(z,mp.mp.dps+5)))
jm = mp.matrix([[mpv(J[i,j]) for j in range(2)] for i in range(2)])
aa,bb,dd = map(mpv,(a,b,D))
tt = mp.asin(mp.sqrt(mpv(r)))/bb
period = mp.pi/bb
def formula(t, av=aa,bv=bb,dv=dd):
    h=mp.sqrt(dv)*abs(mp.sin(bv*t))
    return mp.exp(-av*t)*(mp.sqrt(1+h*h)+h)
def direct(mat,t):
    return mp.svd(mp.expm(mat*t),compute_uv=False)[0]
def rel(x,y): return abs(x-y)/max(mp.mpf(1),abs(x),abs(y))
times=[period*i/120 for i in range(121)]
times += [tt,tt-mp.mpf("1e-8"),tt+mp.mpf("1e-8"),mp.mpf("1e-12"),period-mp.mpf("1e-12"),period+mp.mpf("1e-12")]
times += [n*period+q for n in [1,2,3] for q in [mp.mpf(0),tt,period/2]]
max_error=mp.mpf(0)
for t in times:
    err=rel(direct(jm,t),formula(t))
    max_error=max(max_error,err)
    if err>=tol: raise AssertionError("matrix exponential/SVD crosscheck")
gain=formula(tt)
if not gain>1: raise AssertionError("positive peak")
if not formula(tt-mp.mpf("1e-8"))<gain or not formula(tt+mp.mpf("1e-8"))<gain:
    raise AssertionError("strict local peak")
# Independently obtained maximizing vector and derivative via expm and SVD.
E=mp.expm(jm*tt)
U,S,V=mp.svd(E)
vstar=mp.matrix([V[0,0],V[0,1]])
ustar=mp.matrix([U[0,0],U[1,0]])
witness_error=rel(mp.norm(E*vstar),gain)
derivative=(ustar.T*jm*E*vstar)[0]
if witness_error>=tol or abs(derivative)>=tol:
    raise AssertionError("maximizer and stationary singular-value derivative")
Qm=mp.matrix([[mpv(Q[i,j]) for j in range(2)] for i in range(3)])
direction=Qm*vstar
if abs(sum(direction))>=tol or abs(mp.norm(direction)-1)>=tol:
    raise AssertionError("tangent witness")
normal=mp.matrix([[-aa,-bb],[bb,-aa]])
normal_error=max(rel(direct(normal,t),mp.exp(-aa*t)) for t in [0,tt,period/2,period,2*period])
if normal_error>=tol: raise AssertionError("normal same-spectrum control")
half_error=max(rel(direct(jm,tt+n*period),gain*mp.exp(-aa*n*period)) for n in [1,2,3])
if half_error>=tol: raise AssertionError("halfperiod recurrence")
# Generic stable focus controls: below/on/above threshold; D=0 separately.
controls=[]
for av,bv,k in [(2,1,2),(R(3,4),1,2),(R(1,4),1,2),(1,2,R(1,3)),(1,2,1),(2,1,R(1,2))]:
    av,bv,k=map(mpv,(av,bv,k))
    mat=mp.matrix([[-av,-bv*k],[bv/k,-av]])
    dv=(k-1/k)**2/4
    grows=bv*mp.sqrt(dv)>av
    peak=mp.asin(mp.sqrt((bv*bv*dv-av*av)/(dv*(av*av+bv*bv))))/bv if grows else mp.mpf(0)
    err=max(rel(direct(mat,t),formula(t,av,bv,dv)) for t in [0,mp.mpf("0.37"),mp.pi/(2*bv),mp.pi/bv])
    if err>=tol: raise AssertionError("generic focus control")
    expected=direct(mat,peak)
    if (grows and expected<=1) or (not grows and abs(expected-1)>=tol):
        raise AssertionError("generic threshold control")
    controls.append({"a":str(av),"b":str(bv),"k":str(k),"grows":grows,"peak":str(expected),"error":str(err)})
inputs=[
 PKG/"PROOF.md",PKG/"CASE_SPEC.json",Path(__file__),
 BASE/"_proof-notes/FST2_LYAPUNOV_DOMAIN_2026-10-01/PROOF.md",
 BASE/"papers/FST_II_Chemical_Stability_v3_en.tex",
 BASE/"scripts/autocatalysis_3species.py",
]
result={
 "id":spec["id"],"status":"PASS","claim_status":"PARTIAL_RESULT",
 "completed_utc":datetime.now(timezone.utc).isoformat(),"host":platform.node(),
 "python":sys.version,"sympy":sp.__version__,"mpmath":mp.__version__,"dps":mp.mp.dps,
 "exact_checks":exact,"exact_check_count":len(exact),
 "matrix_exponential_svd_samples":len(times),"max_relative_discrepancy":str(max_error),
 "peak_time":mp.nstr(tt,70),"peak_gain":mp.nstr(gain,70),
 "increase_percent":mp.nstr(100*(gain-1),70),
 "right_singular_vector":[mp.nstr(z,65) for z in vstar],
 "tangent_direction":[mp.nstr(z,65) for z in direction],
 "witness_error":str(witness_error),"peak_derivative":str(derivative),
 "normal_control_error":str(normal_error),"halfperiod_control_error":str(half_error),
 "generic_focus_controls":controls,
 "inputs_sha256":{str(p.relative_to(BASE)).replace(chr(92),"/"):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
 "limits":spec["limitations"],"global_maximum_basis":"analytic proof, not finite grid"
}
print(json.dumps(result,ensure_ascii=False,indent=2))
