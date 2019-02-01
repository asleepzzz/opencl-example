__kernel void mmul(
    const int M,
    const int P,
    const int N,
    __global float* A,
    __global float* B,
    __global float* C)
{
    int k;
    int i = get_global_id(0);
    int j = get_global_id(1);
    if (i<M && j<N) {
        float tmp = 0.0f;
        for (k = 0; k < P; k++){
            tmp += A[i*P + k] * B[k*N + j];
        }
            C[i*N+j] = tmp;
    }
}

