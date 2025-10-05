#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>

#define PI 3.1415926535897

// Function to check if N is a power of two
static bool isPowerOfTwo(unsigned int n) {
    return n && ((n & (n - 1)) == 0);
}

// Function for bit-reversal reordering
void bitReverse(double complex x[], int N) //takes complex no. x and total no. of samples N
{
    int bits = (int)log2(N);  //fixing the no. of stages/ no. of binary digits involved
    
    for (int i = 0; i < N; i++)
    {
        int rev = 0, temp = i;
        for (int b = 0; b < bits; b++) {
            rev = (rev << 1) | (temp & 1);
            temp >>= 1;
        }
        //next step is similar to swapping using 3rd parameter
        if (rev > i) {
            double complex tmp = x[i];
            x[i] = x[rev];
            x[rev] = tmp;
        }
    }
}

// Iterative Radix-2 DIT FFT / IFFT
void fft(double complex x[], int N, int inverse) {
    int P = (int)log2(N); //no of stages
    bitReverse(x, N); //call the bit reversing function

    for (int s = 0; s < P; s++) //iterate through stages
    {
        int m = 1 << (s + 1); //2^(s+1)
        int half = m / 2;     //half the group size
        
        //twiddle factor setup
        double angle = (inverse ? 2.0 : -2.0) * PI / m;
        //(e^(jPI/m * -2.0 or 2.0))
        //-2.0 for FFT, 2.0 for IFFT
        double complex Wm = cos(angle) + I * sin(angle);
        
        //Loop over each group of size m
        //processes N/m butterfly groups per stage
        for (int b = 0; b < N; b += m) 
        {
            //inner butterfly computation
            double complex w = 1.0;
            //w is the rotating twiddle factor, starting at 1 and multiplying by Wm each iteration.
            for (int j = 0; j < half; j++) 
            {
                //u = first element of the butterfly pair
                double complex u = x[b + j];
                
                //t = twiddle factor × second element of the butterfly pair
                double complex t = w * x[b + j + half];
                x[b + j] = u + t;
                x[b + j + half] = u - t;
                
                
                w *= Wm;
            }
        }

        // Display intermediate stage output
        printf("\nStage %d output:\n", s + 1);
        for (int i = 0; i < N; i++) {
            printf("x[%d] = %.2f %+.2fi\n", i, creal(x[i]), cimag(x[i]));
        }
    }

    // Scale down by N for IFFT
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
        //I is included in complex.h library
    }

    printf("\nInput:\n");
    for (int i = 0; i < N; i++) {
        printf("x[%d] = %.2f %+.2fi\n", i, creal(x[i]), cimag(x[i]));
    }

    // Perform FFT
    fft(x, N, 0);

    printf("\nFFT output:\n");
    for (int i = 0; i < N; i++) {
        printf("X[%d] = %.2f %+.2fi\n", i, creal(x[i]), cimag(x[i]));
    }

    // Perform IFFT
    fft(x, N, 1);

    printf("\nIFFT output:\n");
    for (int i = 0; i < N; i++) {
        printf("x[%d] = %.2f %+.2fi\n", i, creal(x[i]), cimag(x[i]));
    }

    return 0;
}

