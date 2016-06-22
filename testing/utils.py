import decimal

def Close( a, b, tol = 0.01 ):
  if isinstance(a,(decimal.Decimal)):
    tol = decimal.Decimal(tol)
  return (a - b)**2 <= 4*tol*tol*(a**2 + b**2)

