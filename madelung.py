import numpy as np

def madelung_CsCl_block(L, offset):
  offx, offy, offz = offset

  A = np.indices([L,L,L])
  I = A[0,:,:,:] + offx # x'座標
  J = A[1,:,:,:] + offy # y'座標
  K = A[2,:,:,:] + offz # z'座標

  ex = slice(offx%2, L, 2); ey = slice(offy%2, L, 2); ez = slice(offz%2, L, 2)
  ox = slice(1-offx%2, L, 2); oy = slice(1-offy%2, L, 2); oz = slice(1-offz%2, L, 2)

  M = np.zeros([L,L,L])

  M[ex,ey,ez] = 1
  M[ex,oy,oz] = 1
  M[ox,ey,oz] = 1
  M[ox,oy,ez] = 1

  M[ox,oy,oz] = -1
  M[ox,ey,ez] = -1
  M[ex,oy,ez] = -1
  M[ex,ey,oz] = -1




  # 末端の処理

  s0 = 1/4 # 90 deg angle
  s1 = 1/6 # 60 deg angle
  s2 = 1/3 # 120 deg angle

  c0 = 1/8  # pi/2 sr
  c1 = 1/24 # pi/6 sr
  c2 = 5/24 # 5pi/6 sr

  Mc = np.ones([L,L,L])
  p = L-1
  Mc[0::p,:,:] = 1/2
  Mc[:,0::p,:] = 1/2
  Mc[:,:,0::p] = 1/2

  Mc[p,p,:] = s2 # AE
  Mc[0,0,:] = s2 # CG
  Mc[p,:,p] = s0 # AD
  Mc[0,:,0] = s0 # FG
  Mc[:,p,p] = s1 # AB
  Mc[:,0,0] = s1 # GH

  Mc[0,p,:] = s1 # BF
  Mc[p,0,:] = s1 # DH
  Mc[0,:,p] = s0 # BC
  Mc[p,:,0] = s0 # EH
  Mc[:,0,p] = s2 # CD
  Mc[:,p,0] = s2 # EF

  Mc[p,p,p] = c0 # A
  Mc[0,p,p] = c1 # B
  Mc[0,0,p] = c2 # C
  Mc[p,0,p] = c0 # D
  Mc[p,p,0] = c2 # E
  Mc[0,p,0] = c0 # F
  Mc[0,0,0] = c0 # G
  Mc[p,0,0] = c1 # H

  r = np.sqrt(3*(I**2 + J**2 + K**2) - 2*(I*J + I*K + J*K))
  if 0 <= -offx < L and 0 <= -offy < L and 0 <= -offz < L:
    Mc[-offx,-offy,-offz] = 0
    r[-offx,-offy,-offz] = 1

  return -((M*Mc)/np.sqrt(r)).sum()*np.sqrt(3)

def madelung_CsCl(L, N):
  alpha = 0
  for x in range(-N, N):
    for y in range(-N, N):
      for z in range(-N, N):
        alpha += madelung_CsCl_block(L+1, (x*L, y*L, z*L))
  return alpha

for l in 1,2,3,4,5,6,7,8,9,10,100:
  print('% 4d, %f' % (2*l+1, madelung_CsCl(l, 1)))

for c in range(3):
  print('% 4d, %f' % ((c+2)*200+1, madelung_CsCl(100, c+2)))


# Result:
#
# $ python3 madelung.py
#    3, 1.608600
#    5, 1.822236
#    7, 1.801646
#    9, 1.803670
#   11, 1.803144
#   13, 1.803381
#   15, 1.803253
#   17, 1.803329
#   19, 1.803281
#   21, 1.803313
#  201, 1.803300
#  401, 1.803300
#  601, 1.803300
#  801, 1.803300
