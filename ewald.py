from numpy import array, indices, ones, zeros, sqrt, pi, exp
from scipy.special import erfc

def ewald(l, eta):
  L = 2*l + 1

  I = indices([L,L,L])[0,:,:,:] - l # x座標
  J = indices([L,L,L])[1,:,:,:] - l # y座標
  K = indices([L,L,L])[2,:,:,:] - l # z座標

  # 格子点の位置ベクトル
  x = zeros([L,L,L,3])
  x[:,:,:,0] = I
  x[:,:,:,1] = J
  x[:,:,:,2] = K

  r2 = I**2 + J**2 + K**2
  r2[l,l,l] = 1 # r^2で割ったときに0除算にならないようにする

  # 逆格子で格子点に向かうベクトルG
  # CsClの場合はxと同じになる
  G = 2*pi*x
  G2 = (2*pi)**2*r2

  V  = exp(0.j - G2/(4*eta))/G2
  V -= exp(-1.j*(G[:,:,:,0]*.5+G[:,:,:,1]*.5+G[:,:,:,2]*.5) - G2/(4*eta))/G2
  V[l,l,l] = 0

  phi1 = (4*pi)*V.sum().real - 2*sqrt(eta/pi)

  Mn = ones([L,L,L])*(-1)
  Mn[0::2,:,:] = 0
  Mn[:,0::2,:] = 0
  Mn[:,:,0::2] = 0
  phi2 = (Mn*erfc(sqrt(eta)*sqrt(r2))/sqrt(r2)).sum()

  return -(phi1 + phi2)*sqrt(3)/2

print(ewald(20, 25))
