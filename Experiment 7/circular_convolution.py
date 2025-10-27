from numpy.fft import fft, ifft
def circular_convolution(x, h):
    N = len(x)
    x = np.array(x)
    h = np.array(h)
    result = np.zeros(N)  # Allocate output array
    
    for n in range(N):
        for m in range(N):
            result[n] += x[m] * h[(n - m) % N]  # Circular index
    return result
    
h = np.random.randn(N)  # Random Gaussian sequence
y_circ = np.real(ifft(fft(y_in) * fft(h)))  # IDFT of product
y_circ_direct = circular_convolution(y_in, h)

verification_err = np.sum(np.abs(y_circ - y_circ_direct)**2)
print("Circular convolution verification error:", verification_err)
