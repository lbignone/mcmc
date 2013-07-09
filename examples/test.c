#include <stdlib.h>
#include <stdio.h>

double* func()
{
    double* arr;

    arr = calloc(2, sizeof(double));

    return arr;
}

int main ()
{
    double* arr;

    arr = func();
    arr[0] = 20.0;
    arr[1] = 21.0;
    printf("%f %f\n", arr[0], arr[1]);
}
