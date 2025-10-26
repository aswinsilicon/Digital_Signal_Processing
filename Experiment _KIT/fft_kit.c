#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

// ---------- Bit Reversal ----------
int bit_reverse(int i, int m) {

    int t,r;
    r = 0;
    for ( t = 0; t < m; t++) {
        r = (r << 1) | (i & 1);
        i >>= 1;
    }
    return r;
}

// ---------- FFT (DIT) ----------
void FFT(complex double *x, int N, int inverse) {
    int m;
    m = log2(N);
    int i,s,L,M,j,p,q,k,r;
    complex double tmp,WL,w,u,v;
    double ang;


    // Phase I: bit-reversal reordering
    for ( i = 0; i < N; i++) {
        r = bit_reverse(i, m);
        if (i < r) {
            tmp = x[i];
            x[i] = x[r];
            x[r] = tmp;
        }
    }

    // Phase II: butterflies
    for (s = 1; s <= m; s++) {
        L = 1 << s;      // block length
        M = L >> 1;      // half block
        ang = (inverse ? 2.0 : -2.0) * M_PI / L;
        WL = cos(ang) + I * sin(ang);

        for ( i = 0; i < N; i += L) {
             w = 1.0 + 0.0*I;
            for ( j = 0; j < M; j++) {
                 p = i + j;
                 q = i + j + M;
                u = x[p];
                 v = w * x[q];
                x[p] = u + v;
                x[q] = u - v;
                w = w * WL;
            }
        }

        // Print stage output
        printf("Stage %d:\n", s);
        for ( k = 0; k < N; k++){
            printf("(%lf,%lf) ", creal(x[k]), cimag(x[k]));
        printf("\n\n");
        }
    }

    // Scale if IFFT
    if (inverse) {
        for ( k = 0; k < N; k++) x[k] /= N;
    }

    printf("FINAL:\n");
    for ( k = 0; k < N; k++){
            printf("(%lf,%lf) ", creal(x[k]), cimag(x[k]));
        printf("\n\n");
    }
}

//---------- Main ----------
int main() {
    int N, choice,i;
    complex double x[8];
    complex double y[8];
    double re, im;
    printf("Enter N (power of 2): ");
    scanf("%d", &N);

    // Check power of 2
    if (N <= 0 || (N & (N - 1)) != 0) {
        printf("Error: N must be power of 2.\n");
        return 1;
    }

    // Input sequence

    printf("Enter real and imaginary parts of input:\n");
    for ( i = 0; i < N; i++) {

        scanf("%lf %lf", &re, &im);
        x[i] = re + I * im;
    }

    printf("1) FFT  2) IFFT  3) Both: ");
    scanf("%d", &choice);
  if (choice == 1) {

        printf("\n--- FFT ---\n");
        FFT(x, N, 0);
    }
    if (choice == 2) {

        printf("\n--- IFFT ---\n");
        FFT(x, N, 1);
    }

    if (choice == 3) {

        for (i = 0; i < N; i++) y[i] = x[i];
        printf("\n--- BOTH ---\n");
        FFT(y, N, 0);
        FFT(y,N,1);

    }
    return 0;
}
