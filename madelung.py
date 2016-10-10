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

  from numpy import arccos, pi
  Mc = np.ones([L,L,L])*2*pi
  Mc[0::L-1,:,:] = pi
  Mc[:,0::L-1,:] = pi
  Mc[:,:,0::L-1] = pi

  theta1 = arccos(1/3)
  theta2 = pi - theta1
  Mc[:,0,0] = theta1
  Mc[0,:,0] = theta1
  Mc[0,0,:] = theta1
  Mc[:,L-1,L-1] = theta1
  Mc[L-1,:,L-1] = theta1
  Mc[L-1,L-1,:] = theta1
  Mc[:,0,L-1] = theta2
  Mc[:,L-1,0] = theta2
  Mc[0,:,L-1] = theta2
  Mc[L-1,:,0] = theta2
  Mc[0,L-1,:] = theta2
  Mc[L-1,0,:] = theta2

  Mc[0,0,0] = (Mc[1,0,0] + Mc[0,1,0] + Mc[0,0,1])/4
  Mc[0,0,L-1] = (Mc[1,0,L-1] + Mc[0,1,L-1] + Mc[0,0,1])/4
  Mc[0,L-1,0] = (Mc[1,L-1,0] + Mc[0,1,0] + Mc[0,L-1,1])/4
  Mc[L-1,0,0] = (Mc[1,0,0] + Mc[L-1,1,0] + Mc[L-1,0,1])/4
  Mc[0,L-1,L-1] = (Mc[1,L-1,L-1] + Mc[0,1,L-1] + Mc[0,L-1,1])/4
  Mc[L-1,0,L-1] = (Mc[1,0,L-1] + Mc[L-1,1,L-1] + Mc[L-1,0,1])/4
  Mc[L-1,L-1,0] = (Mc[1,L-1,0] + Mc[L-1,1,0] + Mc[L-1,L-1,1])/4
  Mc[L-1,L-1,L-1] = (Mc[1,L-1,L-1] + Mc[L-1,1,L-1] + Mc[L-1,L-1,1])/4

  r = np.sqrt(3*(I**2 + J**2 + K**2) - 2*(I*J + I*K + J*K))
  if 0 <= -offx < L and 0 <= -offy < L and 0 <= -offz < L:
    Mc[-offx,-offy,-offz] = 0
    r[-offx,-offy,-offz] = 1

  return -((M*Mc/(2*pi))/np.sqrt(r)).sum()*np.sqrt(3)

for l in 1,2,3,4,5,6,7,8,9,10,100,200,400:
  print('% 4d, %f' % (2*l+1, madelung_CsCl_block(2*l+1, (-l, -l, -l))))

# def madelung_CsCl(L, N):
#   alpha = 0
#   for x in range(-N, N):
#     for y in range(-N, N):
#       for z in range(-N, N):
#         alpha += madelung_CsCl_block(L+1, (x*L, y*L, z*L))
#   print(alpha)

# madelung_CsCl(3, 1)
