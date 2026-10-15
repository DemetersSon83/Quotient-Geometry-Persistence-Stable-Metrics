"""Outward-rounded dyadic intervals and first-order interval differentiation.

All endpoints are integers times 2**(-BITS). Elementary rounding uses exact
integer floor/ceiling; no floating point arithmetic enters an enclosure.
"""
from fractions import Fraction as F
from math import isqrt

BITS = 96
SCALE = 1 << BITS


def ceildiv(a, b):
    if b <= 0:
        raise ValueError('Positive divisor required')
    return -((-a)//b)


class I:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo=0, hi=None, *, raw=False):
        if raw:
            self.lo, self.hi = int(lo), int(lo if hi is None else hi)
        else:
            a, b = F(lo), F(lo if hi is None else hi)
            self.lo = (a.numerator*SCALE)//a.denominator
            self.hi = ceildiv(b.numerator*SCALE, b.denominator)
        if self.lo > self.hi:
            raise ValueError('Empty interval')

    @staticmethod
    def cast(x):
        return x if isinstance(x,I) else I(x)

    def __add__(self, other):
        if isinstance(other, AD):
            return NotImplemented
        o = I.cast(other)
        return I(self.lo+o.lo,self.hi+o.hi,raw=True)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi,-self.lo,raw=True)

    def __sub__(self, other):
        if isinstance(other, AD):
            return NotImplemented
        return self + (-I.cast(other))

    def __rsub__(self, other):
        return I.cast(other)+(-self)

    def __mul__(self, other):
        if isinstance(other, AD):
            return NotImplemented
        o=I.cast(other)
        v=(self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi)
        return I(min(v)//SCALE,ceildiv(max(v),SCALE),raw=True)

    __rmul__=__mul__

    def reciprocal(self):
        if self.lo<=0<=self.hi:
            raise ZeroDivisionError('Interval crosses zero')
        if self.hi<0:
            return -((-self).reciprocal())
        return I((SCALE*SCALE)//self.hi,ceildiv(SCALE*SCALE,self.lo),raw=True)

    def __truediv__(self, other):
        return self*I.cast(other).reciprocal()

    def __rtruediv__(self, other):
        return I.cast(other)*self.reciprocal()

    def square(self):
        a,b=self.lo*self.lo,self.hi*self.hi
        low=0 if self.lo<=0<=self.hi else min(a,b)
        return I(low//SCALE,ceildiv(max(a,b),SCALE),raw=True)

    def sqrt(self):
        if self.lo<0:
            raise ValueError('Negative square-root interval')
        a=isqrt(self.lo*SCALE)
        b=isqrt(self.hi*SCALE)
        if b*b<self.hi*SCALE:b+=1
        return I(a,b,raw=True)

    def intersect(self, other):
        o=I.cast(other)
        return I(max(self.lo,o.lo),min(self.hi,o.hi),raw=True)

    def lower(self):return F(self.lo,SCALE)
    def upper(self):return F(self.hi,SCALE)
    def absmax(self):return F(max(abs(self.lo),abs(self.hi)),SCALE)
    def contains(self,x):return self.lower()<=F(x)<=self.upper()
    def __repr__(self):return f'I({self.lower()}, {self.upper()})'


class AD:
    __slots__=('v','g')

    def __init__(self,value,gradient):
        self.v=I.cast(value)
        self.g=[I.cast(x) for x in gradient]

    @classmethod
    def constant(cls,x,d):return cls(x,[0]*d)

    @classmethod
    def variable(cls,x,d,k):return cls(x,[int(j==k) for j in range(d)])

    def cast(self,x):
        if isinstance(x,AD):
            if len(x.g)!=len(self.g):raise ValueError('Gradient dimension mismatch')
            return x
        return AD.constant(x,len(self.g))

    def __add__(self,other):
        o=self.cast(other)
        return AD(self.v+o.v,[a+b for a,b in zip(self.g,o.g)])

    __radd__=__add__
    def __neg__(self):return AD(-self.v,[-a for a in self.g])
    def __sub__(self,other):return self+(-self.cast(other))
    def __rsub__(self,other):return self.cast(other)+(-self)

    def __mul__(self,other):
        o=self.cast(other)
        return AD(self.v*o.v,[a*o.v+self.v*b for a,b in zip(self.g,o.g)])

    __rmul__=__mul__

    def square(self):
        return AD(self.v.square(),[2*self.v*a for a in self.g])

    def invsqrt(self):
        z=self.v.sqrt().reciprocal()
        derivative=-F(1,2)*z/self.v
        return AD(z,[derivative*a for a in self.g])


def dot(a,b):return sum((x*y for x,y in zip(a,b)),0)


def norm_upper(matrix):
    """min(Frobenius, sqrt(max row sum * max column sum)) enclosure."""
    p=[[I.cast(x).absmax() for x in row] for row in matrix]
    f=sum((x*x for row in p for x in row),F(0))
    r=max(sum(row,F(0)) for row in p)
    c=max(sum((row[j] for row in p),F(0)) for j in range(len(p[0])))
    return I(min(f,r*c)).sqrt().upper()


def matmul(a,b):
    return [[sum((I.cast(a[i][k])*b[k][j] for k in range(len(b))),I(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def row_sum_norm_upper(a):
    """An upper bound on ||a||_(infinity -> 2)."""
    s=sum((sum(abs(F(x)) for x in row)**2 for row in a),F(0))
    return I(s).sqrt().upper()
