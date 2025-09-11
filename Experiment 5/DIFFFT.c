#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

// Define PI for calculations
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Function to print a complex array
void print_complex_array(const char* title, double complex* x, int N) {
    printf("%s\n", title);
    for (int i = 0; i < N; i++) {
        if (cimag(x[i]) < 0) {
            printf("x[%d] = %7.4f - %7.4fj\n", i, creal(x[i]), -cimag(x[i]));
        } else {
            printf("x[%d] = %7.4f + %7.4fj\n", i, creal(x[i]), cimag(x[i]));
        }
    }
    printf("----------------------------------------\n");
}

// Function to perform bit-reversal permutation
// This is used for the output stage of DIF FFT.
void bit_reversal_reorder(double complex* x, int N) {
    int m = log2(N);
    for (int i = 0; i < N; i++) {
        // Efficiently compute the reversed index
        int r = 0;
        int temp_i = i;
        for (int t = 0; t < m; t++) {
            r <<= 1;
            r |= (temp_i & 1);
            temp_i >>= 1;
        }

        // If i < r, swap the elements
        if (i < r) {
            double complex temp = x[i];
            x[i] = x[r];
            x[r] = temp;
        }
    }
}

// Core DIF FFT function
// This implements the algorithm on pages 10-13 of the PDF.
void dif_fft(double complex* x, int N, int is_inverse) {
    printf("\n### Starting %s Computation ###\n\n", is_inverse ? "IFFT" : "DIF FFT");
    
    print_complex_array("Input to Stage 1 (Original Data):", x, N);

    // --- PHASE I: STAGE-WISE BUTTERFLIES (Top-Down, Pages 10-11) ---
    int m = log2(N);
    // Loop for stages s = 1, 2, ..., m
    for (int s = 1; s <= m; s++) {
        // S1: Stage parameters
        int L = N >> (s - 1); // Block length: L = N, N/2, ..., 2 [cite: 208]
        int M = L / 2;      // Half-block length (butterfly span) [cite: 209]
        
        // Base twiddle factor for this stage [cite: 214]
        double angle = 2.0 * M_PI / L;
        double complex WL = cexp(-I * angle * (is_inverse ? -1 : 1)); // Use conjugate for IFFT [cite: 217]

        // S2: Process all blocks across the array [cite: 218]
        for (int i = 0; i < N; i += L) {
            // B1: Initialize twiddle accumulator [cite: 221]
            double complex w = 1.0 + 0.0 * I;
            
            // B2: Loop for each butterfly pair in the block [cite: 223]
            for (int j = 0; j < M; j++) {
                // P1: Define pair indices [cite: 225, 228]
                int p = i + j;
                int q = i + j + M;
                
                // P2: Load the pair [cite: 227, 229]
                double complex u = x[p];
                double complex v = x[q];
                
                // P3 & P4: Compute DIF butterfly and store in-place [cite: 232, 235]
                x[p] = u + v;
                x[q] = (u - v) * w;
                
                // P5: Update twiddle accumulator [cite: 238]
                w *= WL;
            }
        }
        // Print the array after each stage, as per Page 12, Section 2.5 
        char stage_title[50];
        sprintf(stage_title, "Output of Stage %d (L=%d):", s, L);
        print_complex_array(stage_title, x, N);
    }
    
    // --- PHASE II: OUTPUT BIT-REVERSAL (Page 11) ---
    // The array now holds the spectrum in bit-reversed order [cite: 243]
    bit_reversal_reorder(x, N); // [cite: 245]
    
    // --- FINAL SCALING & OUTPUT ---
    if (is_inverse) {
        // Inverse Transform Scaling (Page 11, Section 2.4) [cite: 247]
        for (int k = 0; k < N; k++) {
            x[k] /= N; // [cite: 248]
        }
        print_complex_array("Final IFFT Output (After Bit-Reversal & Scaling):", x, N);
    } else {
        print_complex_array("Final FFT Output (After Bit-Reversal):", x, N);
    }
}

// Function to check if a number is a power of two
int is_power_of_two(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

int main() {
    int k;
    int N;
    printf("Enter the number of points (N), which must be a power of 2: ");
    scanf("%d", &N);

    // Validate N as per page 10, Section 2.2 
    if (!is_power_of_two(N)) {
        printf("Error: N must be a power of 2. Exiting.\n");
        return 1;
    }

    // Allocate memory for the input sequence
    double complex* x = (double complex*)malloc(N * sizeof(double complex));
    
    printf("Enter the real and imaginary parts of the %d-point sequence x[n]:\n", N);
    for (int i = 0; i < N; i++) {
        double real_part, imag_part;
        printf("x[%d] (real imag): ", i);
        scanf("%lf %lf", &real_part, &imag_part);
        x[i] = real_part + imag_part * I;
    }
     printf("\n Enter 0 to compute FFT and 1 to perform IFFT , Enter 2 to perform both  :\n");
    scanf("%d",&k);
    //dit_fft(x, N, 0); // is_inverse = 0 for FFT
    if(k>2) {
        printf(" Error :Enter a value 0,1 or 2  ");
    }
    else if (k==0){
         dif_fft(x, N, 0);
    }
    else if (k==1){
          dif_fft(x, N, 1); // is_inverse = 1 for IFFT
    }
    else {
          printf("\nPerforming FFT computations\n");
          dif_fft(x, N, 0); 
          printf("\nPerforming IFFT computation\n");
          dif_fft(x, N, 1); // is_inverse = 1 for IFFT
    }
    
    
    
    // Free allocated memory
    free(x);

    return 0;
}
