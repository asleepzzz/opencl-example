//------------------------------------------------------------------------------
//
//  PROGRAM: Matrix Multiplication driver
//
//  PURPOSE: This is a driver program to test various ways of computing
//           the product:
//
//                C  = A * B
//
//           A and B are set to constant matrices so we
//           can make a quick test of the multiplication.
//
//  USAGE:   The matrices are constant matrices, square and the order is
//           set as a constant, ORDER (see mult.h).
//
//  HISTORY: Written by Tim Mattson, August 2010
//           Modified by Simon McIntosh-Smith, September 2011
//           Modified by Tom Deakin and Simon McIntosh-Smith, October 2012
//           Updated to C++ Wrapper v1.2.6 by Tom Deakin, August 2013
//           Modified to assume square matrices by Simon McIntosh-Smith, Sep 2014
//
//------------------------------------------------------------------------------

#include "matmul.hpp"
#include "matrix_lib.hpp"
#include "util.hpp"
#include <err_code.h>
#include "device_picker.hpp"

int main(int argc, char *argv[])
{
    int M=800;
    int P=3200;
    int N=1600;                  // A[N][N], B[N][N], C[N][N]
    int size;               // Number of elements in each matrix

    int blockSize = 16;

    double start_time;      // Starting time
    double run_time;        // Timing
    util::Timer timer;      // Timing

    size = M * N;

    std::vector<float> h_A(M*P); // Host memory for Matrix A
    std::vector<float> h_B(P*N); // Host memory for Matrix B
    std::vector<float> h_C(size); // Host memory for Matrix C
    std::vector<float> h_seq(size); // memory for seq results
    cl::Buffer d_a, d_b, d_c;   // Matrices in device memory

//--------------------------------------------------------------------------------
// Create a context and queue
//--------------------------------------------------------------------------------

    try
    {

        cl_uint deviceIndex = 0;
        parseArguments(argc, argv, &deviceIndex);

        // Get list of devices
        std::vector<cl::Device> devices;
        unsigned numDevices = getDeviceList(devices);

        // Check device index in range
        if (deviceIndex >= numDevices)
        {
          std::cout << "Invalid device index (try '--list')\n";
          return EXIT_FAILURE;
        }

        cl::Device device = devices[deviceIndex];

        std::string name;
        getDeviceName(device, name);
        std::cout << "\nUsing OpenCL device: " << name << "\n";

        std::vector<cl::Device> chosen_device;
        chosen_device.push_back(device);
        cl::Context context(chosen_device);
        cl::CommandQueue queue(context, device);

//--------------------------------------------------------------------------------
// Run sequential matmul
//--------------------------------------------------------------------------------


        initmat(M,P,N, h_A, h_B, h_seq,true);

        timer.reset();

        printf("\n===== Sequential, matrix mult (dot prod), order %d on host CPU ======\n",N);
        for(int i = 0; i < COUNT; i++)
        {
            zero_mat(M,N, h_seq);

            start_time = static_cast<double>(timer.getTimeMilliseconds()) / 1000.0;

            seq_mat_mul_sdot(M,P,N, h_A, h_B, h_seq);
            run_time  = (static_cast<double>(timer.getTimeMilliseconds()) / 1000.0) - start_time;
            results(M,P,N, h_seq, run_time);
        }

//--------------------------------------------------------------------------------
// Setup the buffers, initialize matrices, and write them into global memory
//--------------------------------------------------------------------------------

        //  Reset A, B and C matrices (just to play it safe)
        initmat(M,P,N, h_A, h_B, h_C,false);

        d_a = cl::Buffer(context, h_A.begin(), h_A.end(), true);

        d_b = cl::Buffer(context, h_B.begin(), h_B.end(), true);

        d_c = cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * size);
//--------------------------------------------------------------------------------
// OpenCL matrix multiplication ... Naive
//--------------------------------------------------------------------------------


        timer.reset();

        // Create the compute program from the source buffer
        cl::Program program(context, util::loadProgram("matrixc_multiply.cl"), true);

        // Create the compute kernel from the program
        cl::make_kernel<int,int,int, cl::Buffer, cl::Buffer, cl::Buffer> naive_mmul(program, "mmul");


        printf("\n===== OpenCL, matrix mult, C(i,j) per work item, order %d ======\n",N);

        // Do the multiplication COUNT times
        for (int i = 0; i < COUNT; i++)
        {
            zero_mat(M,N, h_C);

            start_time = static_cast<double>(timer.getTimeMilliseconds()) / 1000.0;

            // Execute the kernel over the entire range of C matrix elements ... computing
            // a dot product for each element of the product matrix.  The local work
            // group size is set to NULL ... so I'm telling the OpenCL runtime to
            // figure out a local work group size for me.
            cl::NDRange global(M,N);
            cl::NDRange local(blockSize,blockSize);
            //cl::NDRange offset(M/2,N/2);
            naive_mmul(
                    cl::EnqueueArgs(queue,global,local),
                    M,P,N, d_a, d_b, d_c);

            queue.finish();

            run_time  = (static_cast<double>(timer.getTimeMilliseconds()) / 1000.0) - start_time;

            cl::copy(queue, d_c, h_C.begin(), h_C.end());

            results(M,P,N, h_C, run_time);

        } // end for loop

        error(M, N, h_seq,h_C);

//--------------------------------------------------------------------------------
// OpenCL matrix multiplication ... Sub block
//--------------------------------------------------------------------------------

        timer.reset();


        //d_c = cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * size);
        cl::Buffer d_c2 = cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * size);

        // Create the compute program from the source buffer
        cl::Program program2(context, util::loadProgram("matrixc_multiply_subblock.cl"), true);
        //cl::Program program2(context, util::loadProgram("matrixc_multiply.cl"), true);
        // Create the compute kernel from the program
        cl::make_kernel<int,int,int, cl::Buffer, cl::Buffer, cl::Buffer> sub_block_mmul(program2, "mmul_sub_block");
        //cl::make_kernel<int,int,int, cl::Buffer, cl::Buffer, cl::Buffer> sub_block_mmul(program2, "mmul");

        printf("\n===== OpenCL,sub block matrix mult, C(i,j) per work item, order %d ======\n",N);

        // Do the multiplication COUNT times
        for (int i = 0; i < COUNT; i++)
        {
            zero_mat(M,N, h_C);

            start_time = static_cast<double>(timer.getTimeMilliseconds()) / 1000.0;

            // Execute the kernel over the entire range of C matrix elements ... computing
            // a dot product for each element of the product matrix.  The local work
            // group size is set to NULL ... so I'm telling the OpenCL runtime to
            // figure out a local work group size for me.
            cl::NDRange global(M,N);
            cl::NDRange local(blockSize,blockSize);
            //cl::NDRange offset(M/2,N/2);
            sub_block_mmul(
                    cl::EnqueueArgs(queue,global,local),
                    M,P,N, d_a, d_b, d_c2);

            queue.finish();

            run_time  = (static_cast<double>(timer.getTimeMilliseconds()) / 1000.0) - start_time;

            cl::copy(queue, d_c2, h_C.begin(), h_C.end());

            results(M,P,N, h_C, run_time);

        } // end for loop

        error(M, N, h_seq,h_C);


    } catch (cl::Error err)
    {
        std::cout << "Exception\n";
        std::cerr << "ERROR: "
                  << err.what()
                  << "("
                  << err_code(err.err())
                  << ")"
                  << std::endl;
    }

    return EXIT_SUCCESS;
}
