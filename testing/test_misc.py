# miscellaneous tests to figure out how stuff works

from utils import Close

import math,cmath
from math import sqrt, degrees, radians
from cmath import phase

def test_complex():

  r1 = 2 + 3j
  r2 = -1 + 2j

  dr = r2 - r1

  assert Close( dr.real , -3 )
  assert Close( dr.imag , -1 )
  assert Close( abs(dr) , sqrt(10) )

  p1 = 2 + 2j
  p2 = 4 + 4j

  dp = p2 - p1

  # rotations

  assert Close( dp.real, 2 )
  assert Close( dp.imag, 2 )
  assert Close( degrees(phase(dp)), 45. )

  dp *= cmath.rect( 1, radians(90) )
  assert Close( dp.real,-2 )
  assert Close( dp.imag, 2 )
  assert Close( degrees(phase(dp)), 90+45 )

  dp *= cmath.rect( 1, -radians(90) )
  assert Close( dp.real, 2 )
  assert Close( dp.imag, 2 )
  assert Close( degrees(phase(dp)), 45. )

  dp *= cmath.rect( 1, -radians(90) )
  assert Close( dp.real, 2 )
  assert Close( dp.imag,-2 )
  assert Close( degrees(phase(dp)), -45 )

  dp *= cmath.rect( 1, -radians(90) )
  assert Close( dp.real,-2 )
  assert Close( dp.imag,-2 )
  assert Close( degrees(phase(dp)), -90-45 )


  c1 = 2 + 3j
  c2 = 4 + 5j

  c3 = c1*c2
  assert Close(c3.real, 8-15)
  assert Close(c3.imag, 10+12)

  c3 = c1*c2.conjugate()
  assert Close(c3.real, 8+15)
  assert Close(c3.imag, 12-10)
