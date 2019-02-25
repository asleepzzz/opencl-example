__kernel void mmul_sub_block(
    const int M,
    const int P,
    const int N,
    __global float* A,
    __global float* B,
    __global float* C)
{
    const int row = get_local_id(0);
    const int col = get_local_id(1);
    const int TS = 16;
    int k;
    const int globalRow = TS*get_group_id(0) + row;
    const int globalCol = TS*get_group_id(1) + col;
    __local float Asub[TS][TS];
    __local float Bsub[TS][TS];

    float acc = 0.0f;
    const int numTiles = P/TS;
    for (int t=0 ; t < numTiles; t++) {
         const int tiledRow = TS*t + row;
         const int tiledCol = TS*t + col;
         Asub[row][col] = A[globalRow*P + tiledCol];
         Bsub[row][col] = B[tiledRow*N + globalCol];
         

         barrier(CLK_LOCAL_MEM_FENCE);
         for (int k=0; k<TS; k++) {
             acc += Asub[row][k] * Bsub[k][col];
         }
         barrier(CLK_LOCAL_MEM_FENCE);
    }
    C[globalRow*N + globalCol] = acc;
}

