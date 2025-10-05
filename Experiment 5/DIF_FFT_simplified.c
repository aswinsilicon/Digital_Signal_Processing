#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>

#define PI 3.1415926535897

// Check if N is a power of 2
static bool isPowerOfTwo(unsigned int n) {
    return n && ((n & (n - 1)) == 0);
}

// Bit-reversal reordering (applied AFTER DIF FFT stages)
void bitReverse(double complex x[], int N) {
    int bits = (int)log2(N);
    for (int i = 0; i < N; i++) {
        int rev = 0, temp = i;
        for (int b = 0; b < bits; b++) {
            rev = (rev << 1) | (temp & 1);
            temp >>= 1;
        }
        if (rev > i) {
            double complex tmp = x[i];
            x[i] = x[rev];
            x[rev] = tmp;
        }
    }
}

// Iterative Radix-2 DIF FFT / IFFT
void fft(double complex x[], int N, int inverse) {
    int P = (int)log2(N);

    for (int s = 0; s < P; s++) {
        int m = 1 << (P - s);   // Stage butterfly size
        int half = m / 2;
        double angle = (inverse ? 2.0 : -2.0) * PI / m;
        double complex Wm = cos(angle) + I * sin(angle);

        for (int b = 0; b < N; b += m) {
            double complex w = 1.0;
            for (int j = 0; j < half; j++) {
                double complex u = x[b + j];
                double complex v = x[b + j + half];
                x[b + j]        = u + v;
                x[b + j + half] = (u - v) * w;
                w *= Wm;  // Update twiddle factor
            }
        }

        // Display stage output
        printf("\nStage %d output:\n", s + 1);
        for (int i = 0; i < N; i++) {
            printf("x[%d] = %.2f %+.2fi\n", i, creal(x[i]), cimag(x[i]));
        }
    }

    // Bit-reverse the output to normal order
    bitReverse(x, N);

    // IFFT scaling
    if (inverse) {
        for (int i = 0; i < N; i++)
            x[i] /= N;
    }
}

int main() {
    unsigned int N;
    printf("Enter value for N: ");
    scanf("%u", &N);

    if (!isPowerOfTwo(N)) {
        printf("Error: N must be a power of 2!\n");
        return 1;
    }

    double complex x[N];
    double real, imag;

    // Read input samples
    for (int i = 0; i < N; i++) {
        printf("Enter real part of x[%d]: ", i);
        scanf("%lf", &real);
        printf("Enter imag part of x[%d]: ", i);
        scanf("%lf", &imag);
        x[i] = real + imag * I;
    }

    printf("\nInput:\n");
    for (int i = 0; i < N; i++) {
        printf("x[%d] = %.2f %+.2fi\n", i, creal(x[i]), cimag(x[i]));
    }

    // Perform DIF FFT
    fft(x, N, 0);
    printf("\nFFT output:\n");
    for (int i = 0; i < N; i++) {
        printf("X[%d] = %.2f %+.2fi\n", i, creal(x[i]), cimag(x[i]));
    }

    // Perform DIF IFFT
    fft(x, N, 1);
    printf("\nIFFT output:\n");
    for (int i = 0; i < N; i++) {
        printf("x[%d] = %.2f %+.2fi\n", i, creal(x[i]), cimag(x[i]));
    }

    return 0;
}

