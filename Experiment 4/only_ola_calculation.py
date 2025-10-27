import numpy as np

def overlap_add(x, h, L):
    """
    Overlap-Add convolution (FFT-based) for complex sequences.
    x : 1-D array_like (complex) - input signal
    h : 1-D array_like (complex) - impulse response
    L : int - block length (number of new input samples per block)
    Returns full linear convolution y (dtype=complex).
    """
    x = np.asarray(x, dtype=complex)
    h = np.asarray(h, dtype=complex)

    M = len(h)
    Nx = len(x)
    N = L + M - 1
    P = int(2 ** np.ceil(np.log2(N)))   # FFT length

    H_fft = np.fft.fft(h, P)
    y = np.zeros(Nx + M - 1, dtype=complex)

    # process blocks of length L
    for start in range(0, Nx, L):
        end = min(start + L, Nx)
        x_block = x[start:end]
        x_block_padded = np.pad(x_block, (0, P - len(x_block)), mode='constant')

        X_block = np.fft.fft(x_block_padded)
        Y_block = X_block * H_fft
        y_block = np.fft.ifft(Y_block)   # length P

        # add valid portion of y_block into output (handle edge)
        end_y = min(start + P, len(y))
        valid_len = end_y - start
        y[start:end_y] += y_block[:valid_len]

    return y

# ----------------- example usage -----------------
if __name__ == "__main__":
    # example complex sequences
    x = np.array([1+1j, 2+2j, 3+0j, 4-1j, 2+1j])     # input
    h = np.array([0.5+0.5j, 1-1j, 0.5+0.5j])         # filter
    L = 4                                            # block length

    y = overlap_add(x, h, L)

    print("x:", x)
    print("h:", h)
    print("y (length={}):".format(len(y)))
    for i, v in enumerate(y):
        print(f"y[{i}] = {v.real:.4f} {v.imag:+.4f}j")
