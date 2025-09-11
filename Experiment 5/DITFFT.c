#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
void bit_reversal_reorder(double complex* x, int N) {
    int m = log2(N);
    for (int i = 0; i < N; i++) {
        // R1-R3: Efficiently compute the reversed index
        int r = 0;
        int temp_i = i;
        for (int t = 0; t < m; t++) {
            r <<= 1;
            r |= (temp_i & 1);
            temp_i >>= 1;
        }

        // B2: If i < r, swap the elements
        if (i < r) {
            double complex temp = x[i];
            x[i] = x[r];
            x[r] = temp;
        }
    }
}

// Core DIT FFT function
void dit_fft(double complex* x, int N, int is_inverse) {
    printf("\n### Starting %s Computation ###\n\n", is_inverse ? "IFFT" : "FFT");
    
    // --- PHASE I: BIT-REVERSAL REORDERING
    bit_reversal_reorder(x, N);
    print_complex_array("After Bit-Reversal (Input to Stage 1):", x, N);

    // --- PHASE II: STAGE-WISE BUTTERFLIES 
    int m = log2(N);
    for (int s = 1; s <= m; s++) {
        // S1: Stage parameters
        int L = 1 << s; // Block length: L = 2^s
        int M = L / 2;  // Half-block length: M = L/2
        
        // Base twiddle factor for this stage
        double angle = 2.0 * M_PI / L * (is_inverse ? -1 : 1);
        double complex WL = cexp(-I * angle);

        // S2: Process blocks across the array
        for (int i = 0; i < N; i += L) {
            // B1: Initialize twiddle accumulator
            double complex w = 1.0 + 0.0 * I;
            
            // B2: Loop for each butterfly pair in the block
            for (int j = 0; j < M; j++) {
                // P1: Define pair indices
                int p = i + j;
                int q = i + j + M;
                
                // P2 & P3: Load pair and apply twiddle
                double complex u = x[p];
                double complex v = w * x[q];
                
                // P4: Compute butterfly outputs (in-place)
                x[p] = u + v;
                x[q] = u - v;
                
                // P5: Update twiddle accumulator
                w *= WL;
            }
        }
        // Print the array after each stage
        char stage_title[50];
        sprintf(stage_title, "Output of Stage %d:", s);
        print_complex_array(stage_title, x, N);
    }
    
        //INVERSE TRANSFORM SCALING 
    if (is_inverse) {
        for (int k = 0; k < N; k++) {
            x[k] /= N;
        }
        print_complex_array("Final IFFT Output (After Scaling):", x, N);
    } else {
        print_complex_array("Final FFT Output:", x, N);
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

    // Validate N as per page 5, Section 1.2
    if (!is_power_of_two(N)) {
        printf("Error: N must be a power of 2. Exiting.\n");
        return 1;
    }

    // Allocate memory for the input sequence and its copy for IFFT
    double complex* x = (double complex*)malloc(N * sizeof(double complex));
    double complex* x_copy = (double complex*)malloc(N * sizeof(double complex));

    printf("Enter the real and imaginary parts of the %d-point sequence x[n]:\n", N);
    for (int i = 0; i < N; i++) {
        double real_part, imag_part;
        printf("x[%d] (real imag): ", i);
        scanf("%lf %lf", &real_part, &imag_part);
        x[i] = real_part + imag_part * I;
        x_copy[i] = x[i]; // Keep a copy for IFFT
    }
    
    printf("\nOriginal Input Sequence:\n");
    print_complex_array("x[n]", x, N);

    // --- Perform FFT ---
    printf("\n Enter 0 to compute FFT and 1 to perform IFFT , Enter 2 to perform both  :\n");
    scanf("%d",&k);
    //dit_fft(x, N, 0); // is_inverse = 0 for FFT
    if(k>2) {
        printf(" Error :Enter a value 0,1 or 2  ");
    }
    else if (k==0){
          dit_fft(x, N, 0); 
    }
    else if (k==1){
          dit_fft(x, N, 1); // is_inverse = 1 for IFFT
    }
    else {
          printf("\nPerforming FFT computations\n");
          dit_fft(x, N, 0); 
          printf("\nPerforming IFFT computation\n");
          dit_fft(x, N, 1); // is_inverse = 1 for IFFT
    }

 
    
  
    
    // Free allocated memory
    free(x);
    free(x_copy);

    return 0;
}

