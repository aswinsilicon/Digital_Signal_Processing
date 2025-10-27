
def linear_convolution(x_impulse,h):
  N=len(x)
  M=len(h)
  y = np.zeros(N+M-1) #since length of the final result must be of length L =N+M-1
  for i in range(len(y)):
    for j in range(N):
      if 0<= (i-j)<= M:
        y[i]+=x[j]*h[i-j]
  return y
