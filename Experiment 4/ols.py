import numpy as np
import matplotlib.pyplot as plt
import time

# --------------------------------------------------------------
# Overlap-Save FIR Filtering Function
# --------------------------------------------------------------

#noisy_temperature_data : input signal
#h : FIR filter impulse response
#L_ols_valid : no. of valid samples per block (without overlap)
#P : FFT length
def overlap_save(noisy_temperature_data, h, L_ols_valid, P):
    Nx = len(noisy_temperature_data)   # Length of input signal
    M = len(h)                         # Length of FIR filter

    #check whether the FFT size is big
    if P < M - 1:
        raise ValueError("FFT length P must be >= M - 1")

    # Pre-compute FFT of FIR filter (for frequency-domain convolution) 
    #with zero padding
    H_fft_ols = np.fft.fft(h, P)
    
    
    # Initialize output array (final filtered signal)
    y_ols = np.zeros(Nx + M - 1)

    #keeps the last M-1 samples of the prev block
    # Initial overlap buffer (M - 1 zeros)
    overlap_buffer = np.zeros(M - 1)

    #BLOCK PROCESSING
    # Start timing
    start_time = time.time()
    num_blocks = int(np.ceil(Nx / L_ols_valid)) 
    # Total number of data blocks

    # Process each block
    for i in range(num_blocks):
        start_idx = i * L_ols_valid
        end_idx = min(start_idx + L_ols_valid, Nx)

        # Get current block and prepend overlap
        current_samples = noisy_temperature_data[start_idx:end_idx]
        x_block = np.concatenate((overlap_buffer, current_samples))

        # Zero-pad if block length < P
        if len(x_block) < P:
            x_block = np.pad(x_block, (0, P - len(x_block)))

        # FFT-based filtering (frequency-domain multiplication)
        X_block = np.fft.fft(x_block)
        Y_block = X_block * H_fft_ols
        y_circular = np.fft.ifft(Y_block).real

        # Extract valid part of output (discard corrupted overlap samples)
        #discards the M-1 samples from the beginning
        valid_output = y_circular[M - 1: M - 1 + L_ols_valid]

        # Store valid output into final output array
        output_start_idx = i * L_ols_valid
        output_end_idx = min(output_start_idx + len(valid_output), len(y_ols))
        y_ols[output_start_idx:output_end_idx] = valid_output[:output_end_idx - output_start_idx]

        # Update overlap buffer for next block
        overlap_buffer = x_block[P - (M - 1):P]

    end_time = time.time()
    print(f"OLS filtering completed in {end_time - start_time:.4f} seconds")

    return y_ols[:Nx + M - 1]


# --------------------------------------------------------------
# Example: Apply Overlap-Save Filtering
# --------------------------------------------------------------

# Generate example noisy signal (temperature data)
Fs = 1  # 1 sample/hour
t = np.arange(0, 200, 1/Fs)
true_temperature = 25 + 2*np.sin(2*np.pi*0.01*t)
noise = np.random.normal(0, 0.8, len(t))
noisy_temperature_data = true_temperature + noise

# Design a simple FIR lowpass filter (e.g., moving average)
M = 33
h = np.ones(M) / M

# Set FFT Length (power of 2) and compute L_ols_valid
P = 512
L_ols_valid = P - (M - 1)

# Apply Overlap-Save Filtering
filtered_ols_output = overlap_save(noisy_temperature_data, h, L_ols_valid, P)

# --------------------------------------------------------------
# Plot Results
# --------------------------------------------------------------
plt.figure(figsize=(10,5))
plt.plot(noisy_temperature_data, label='Noisy Input Signal', alpha=0.6)
plt.plot(filtered_ols_output, label='Filtered Output (OLS)', linewidth=2)
plt.xlabel("Time (Hours)")
plt.ylabel("Temperature (°C)")
plt.title("Temperature Signal FIR Filtering using Overlap-Save")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
