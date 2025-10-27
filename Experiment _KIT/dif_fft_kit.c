#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

// ---------- Bit Reversal ----------
int bit_reverse(int i, int m) {
    int t, r = 0;
    for (t = 0; t < m; t++) {
        r = (r << 1) | (i & 1);
        i >>= 1;
    }
    return r;
}

// ---------- DIF FFT Function ----------
void FFT_DIF(complex double *x, int N, int inverse) {
    int m = log2(N);
    int i, s, L, M, j, k, r;
    complex double tmp, WL, w, u, v;
    double ang;

    // Phase I: butterflies
    for (s = m; s >= 1; s--) {
        L = 1 << s;      // block length
        M = L >> 1;      // half block
        ang = (inverse ? 2.0 : -2.0) * M_PI / L;
        WL = cos(ang) + I * sin(ang);

        for (i = 0; i < N; i += L) {
            w = 1.0 + 0.0 * I;
            for (j = 0; j < M; j++) {
                int p = i + j;
                int q = i + j + M;

                u = x[p];
                v = x[q];

                x[p] = u + v;
                x[q] = (u - v) * w;

                w *= WL;
            }
        }

        // Print stage output
        printf("Stage %d:\n", m - s + 1);
        for (k = 0; k < N; k++) {
            printf("(%lf,%lf)\n", creal(x[k]), cimag(x[k]));
        }
        printf("\n");
    }

    // Phase II: bit reversal (reordering at the end)
    for (i = 0; i < N; i++) {
        r = bit_reverse(i, m);
        if (r > i) {
            tmp = x[i];
            x[i] = x[r];
            x[r] = tmp;
        }
    }

    // Scale for IFFT
    if (inverse) {
        for (k = 0; k < N; k++)
            x[k] /= N;
    }

    // Print final output
    printf("FINAL:\n");
    for (k = 0; k < N; k++) {
        printf("(%lf,%lf)\n", creal(x[k]), cimag(x[k]));
    }
    printf("\n");
}

// ---------- Main ----------
int main() {
    int N, choice, i;
    complex double x[8], y[8];
    double re, im;

    printf("Enter N (power of 2): ");
    scanf("%d", &N);

    if (N <= 0 || (N & (N - 1)) != 0) {
        printf("Error: N must be power of 2.\n");
        return 1;
    }

    printf("Enter real and imaginary parts of input:\n");
    for (i = 0; i < N; i++) {
        scanf("%lf %lf", &re, &im);
        x[i] = re + I * im;
    }

    printf("1) FFT  2) IFFT  3) Both: ");
    scanf("%d", &choice);

    if (choice == 1) {
        printf("\n--- DIF FFT ---\n");
        FFT_DIF(x, N, 0);
    }
    else if (choice == 2) {
        printf("\n--- DIF IFFT ---\n");
        FFT_DIF(x, N, 1);
    }
    else if (choice == 3) {
        for (i = 0; i < N; i++) y[i] = x[i];
        printf("\n--- BOTH (FFT then IFFT) ---\n");
        FFT_DIF(y, N, 0);
        FFT_DIF(y, N, 1);
    }

    return 0;
}
