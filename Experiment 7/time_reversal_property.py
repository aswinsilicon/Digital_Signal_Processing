from numpy.fft import fft, ifft

N = len(y_in)
y_in_rev = y_in.copy()
y_in_rev[1:] = y_in_rev[1:][::-1]  # Reverse all except index 0

Y = fft(y_in)  # DFT of original signal

k = np.arange(N)
Y_N_refl_shift = Y[-k % N] #finding Y[-k (mod N)] for each value of k
Y_rev = fft(y_in_rev) #compute the DFT of time reversed signal y_in_rev
