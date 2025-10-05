import numpy as np
import matplotlib.pyplot as plt
import time

# --------------------------------------------------------------
# Example input: synthetic noisy temperature data
# --------------------------------------------------------------
Fs = 1  # 1 sample per hour
t = np.arange(0, 200, 1/Fs)
true_temperature = 25 + 2 * np.sin(2 * np.pi * 0.01 * t)
noise = np.random.normal(0, 0.8, len(t))
noisy_temperature_data = true_temperature + noise

# --------------------------------------------------------------
# FIR Filter (simple moving-average lowpass filter)
# --------------------------------------------------------------
M = 33
h = np.ones(M) / M   # FIR impulse response

# --------------------------------------------------------------
# Block & FFT Parameters
# --------------------------------------------------------------
L = 1024                           # Block length (input samples per block)
N = L + M - 1                      # Output length per linear convolution block
P = int(2 ** np.ceil(np.log2(N)))  # FFT length = next power of 2 ≥ N

print(f"L = {L}, N = {N}, M = {M}, P = {P}")

# --------------------------------------------------------------
# Overlap-Add FIR Filtering Function
# --------------------------------------------------------------
def overlap_add(x, h, L, P):
    Nx = len(x)        # Length of input signal
    M = len(h)         # Length of filter
    N_blocks = int(np.ceil(Nx / L))  # Number of blocks

    # Step 1: FFT of zero-padded filter
    H_fft = np.fft.fft(h, P)

    # Step 2: Initialize output array
    y_ola = np.zeros(Nx + M - 1)

    # Step 3: Start time measurement
    start_time = time.time()

    # Step 4: Process each block
    for i in range(N_blocks):
        start_idx = i * L
        end_idx_x = min(start_idx + L, Nx)

        # Extract input block and zero-pad to length P
        x_block = x[start_idx:end_idx_x]
        x_block_padded = np.pad(x_block, (0, P - len(x_block)))

        # FFT-based convolution
        X_block = np.fft.fft(x_block_padded)
        Y_block = X_block * H_fft
        y_block = np.real(np.fft.ifft(Y_block))

        # Add block output (overlap-add step)
        end_idx_y = min(start_idx + P, len(y_ola))
        valid_length = end_idx_y - start_idx
        y_ola[start_idx:end_idx_y] += y_block[:valid_length]

    # Step 5: End timing
    end_time = time.time()
    exec_time = end_time - start_time

    print(f"Overlap-Add Execution Time: {exec_time:.6f} seconds")

    return y_ola

# --------------------------------------------------------------
# Apply Overlap-Add Filtering
# --------------------------------------------------------------
filtered_ola_output = overlap_add(noisy_temperature_data, h, L, P)

# --------------------------------------------------------------
# Plot Results
# --------------------------------------------------------------
plt.figure(figsize=(10,5))
plt.plot(noisy_temperature_data, label='Noisy Input', alpha=0.6)
plt.plot(filtered_ola_output, label='Filtered Output (OLA)', linewidth=2)
plt.title('Overlap-Add Convolution Output')
plt.xlabel('Sample Index (Time in Hours)')
plt.ylabel('Temperature (°C)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
