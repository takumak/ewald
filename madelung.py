import numpy as np

def madelung_CsCl(l):
  L = 2*l + 1

  A = np.indices([L,L,L])
  I = A[0,:,:,:] - l # x'座標
  J = A[1,:,:,:] - l # y'座標
  K = A[2,:,:,:] - l # z'座標

  e = slice(l%2, L, 2)
  o = slice(1-l%2, L, 2)

  M = np.zeros([L,L,L])

  M[e,e,e] = 1
  M[e,o,o] = 1
  M[o,e,o] = 1
  M[o,o,e] = 1

  M[o,o,o] = -1
  M[o,e,e] = -1
  M[e,o,e] = -1
  M[e,e,o] = -1

  r = np.sqrt(3*(I**2 + J**2 + K**2) - 2*(I*J + I*K + J*K))
  r[l,l,l] = 1

  M[l,l,l] = 0
  M[0::L-1,:,:] /= 2
  M[:,0::L-1,:] /= 2
  M[:,:,0::L-1] /= 2

  print(-(M/np.sqrt(r)).sum()*np.sqrt(3))

madelung_CsCl(1)
madelung_CsCl(2)
madelung_CsCl(3)
madelung_CsCl(4)
madelung_CsCl(5)
madelung_CsCl(6)
madelung_CsCl(7)
madelung_CsCl(8)
madelung_CsCl(9)
madelung_CsCl(10)
madelung_CsCl(100)
