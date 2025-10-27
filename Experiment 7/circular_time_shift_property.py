from numpy.fft import fft, ifft

n0 = 12 #delaying y_in by 12
y_in_shifted = np.roll(y_in, n0)
Y_shifted_dft = fft(y_in_shifted)

k = np.arange(N)
Y_expected_shift_dft = Y * np.exp(-1j * 2 * np.pi * k * n0 / N)

#PTQ_13
verification_err = np.sum(np.abs(Y_shifted_dft - Y_expected_shift_dft)**2)
print("Circular time shift verification error:", verification_err)
