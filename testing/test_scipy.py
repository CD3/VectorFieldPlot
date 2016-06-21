

import scipy as sc
import numpy as np

from scipy.integrate import odeint


def test_integrator():

  def v(r,s):
    return sc.array( [r,-2*r] )


  r0 = sc.array([0,0])

  s = np.linspace(0,10,100)

  r = odeint( v, r0, s )

  print r
