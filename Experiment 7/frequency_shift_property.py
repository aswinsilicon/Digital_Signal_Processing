from numpy.fft import fft, ifft
k0 = 3 #setting the frequency shift to 3
n = np.arange(N)
y_in_freq_shifted = y_in * np.exp(1j * 2 * np.pi * k0 * n / N)

Y_freq_shifted = fft(y_in_freq_shifted)
Y_expected_freq_shift = np.roll(Y, k0)

#PTQ_14
verification_err = np.sum(np.abs(Y_freq_shifted - Y_expected_freq_shift)**2)
print("Frequency shift verification error:", verification_err)
