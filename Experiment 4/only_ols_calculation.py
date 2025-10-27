import numpy as np

def overlap_save(x, h, L):
    """
    Overlap-Save convolution (FFT-based) for complex sequences.
    x : input signal (complex)
    h : impulse response (complex)
    L : block length (number of new output samples per block)
    Returns full linear convolution y (complex).
    """
    x = np.asarray(x, dtype=complex)
    h = np.asarray(h, dtype=complex)

    M = len(h)
    Nx = len(x)
    N = L + M - 1
    P = int(2 ** np.ceil(np.log2(N)))

    H_fft = np.fft.fft(h, P)
    x_padded = np.pad(x, (M - 1, 0))
    y = np.zeros(Nx + M - 1, dtype=complex)

    for i in range(0, Nx, L):
        x_block = x_padded[i : i + P]
        if len(x_block) < P:
            x_block = np.pad(x_block, (0, P - len(x_block)))

        X_block = np.fft.fft(x_block)
        Y_block = X_block * H_fft
        y_block = np.fft.ifft(Y_block)

        valid = y_block[M-1 : M-1+L]
        start_y = i
        end_y = min(start_y + L, len(y))
        y[start_y:end_y] = valid[:end_y - start_y]

    return y


# ---------------- Example usage ----------------
if __name__ == "__main__":
    x = np.array([1+1j, 2+2j, 3+0j, 4-1j, 2+1j])     # input
    h = np.array([0.5+0.5j, 1-1j, 0.5+0.5j])         # filter
    L = 4                                            # block length

    y = overlap_save(x, h, L)

    print("x:", x)
    print("h:", h)
    print("y (length={}):".format(len(y)))
    for i, v in enumerate(y):
        print(f"y[{i}] = {v.real:.4f} {v.imag:+.4f}j")
