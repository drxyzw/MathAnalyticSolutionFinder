import gdb, re
from decimal import Decimal, localcontext

EXP_ZERO = 1 - (1 << 63)   # MPFR special exponents (64-bit long)
EXP_NAN = EXP_ZERO + 1
EXP_INF = EXP_ZERO + 2

MAX_PREC = 1 << 22         # sanity limit: bigger means uninitialized memory
EXP_LIMIT = 1 << 62

def mpfr_to_string(m):
    try:
        return _mpfr_to_string(m)
    except Exception as e:
        return '<unreadable: %s>' % e

def _mpfr_to_string(m):
    prec = int(m['_mpfr_prec']); sign = int(m['_mpfr_sign']); exp = int(m['_mpfr_exp'])
    if not (2 <= prec <= MAX_PREC) or sign not in (1, -1):
        return '<uninitialized>'
    neg = '-' if sign < 0 else ''
    if exp == EXP_ZERO: return neg + '0'
    if exp == EXP_NAN: return 'nan'
    if exp == EXP_INF: return neg + 'inf'
    if abs(exp) >= EXP_LIMIT:
        return '<uninitialized>'
    n = (prec + 63) // 64
    d = m['_mpfr_d']
    ul = gdb.lookup_type('unsigned long')
    mant = 0
    for i in range(n):
        mant |= int(d[i].cast(ul)) << (64 * i)
    digits = int(prec * 0.30103) + 3
    with localcontext() as c:
        c.prec = digits + 20
        c.Emax = 999999999999999999; c.Emin = -999999999999999999
        val = Decimal(mant) * (Decimal(2) ** (exp - 64 * n))
        return neg + f"{val:.{digits}g}"

class MpfrPrinter:
    def __init__(self, val): self.v = val
    def to_string(self): return mpfr_to_string(self.v)

class MprealPrinter:
    def __init__(self, val): self.v = val
    def to_string(self): return mpfr_to_string(self.v['mp'][0])

def lookup(val):
    t = val.type.strip_typedefs()
    name = t.tag or t.name or ''
    if re.fullmatch(r'(mpfr::)?mpreal', name): return MprealPrinter(val)
    if name == '__mpfr_struct': return MpfrPrinter(val)
    return None

gdb.pretty_printers.append(lookup)
