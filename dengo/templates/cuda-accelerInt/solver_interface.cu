/**
 * \file
 * \brief Interface implementation for GPU solvers to be called as a library
 *
 * \author Nicholas Curtis
 * \date 03/09/2015
 *
 * Contains initialization, integration and cleanup functions
 */

#include "solver_interface.cuh"

#ifdef GENERATE_DOCS
namespace genericcu {
#endif

//! Padded # of ODEs to solve
int padded;
//! The solver memory structs
solver_memory* host_solver, *device_solver;
//! The mechanism memory structs
mechanism_memory* host_mech, *device_mech;
//! block and grid sizes
dim3 dimBlock, dimGrid;
//! result flag
int* result_flag;
//! temorary storage
double* y_temp;

/**
 * \brief A convienience method to copy memory between host pointers of different pitches, widths and heights.
 *        Enables easier use of CUDA's cudaMemcpy2D functions.
 *
 * \param[out]              dst             The destination array
 * \param[in]               pitch_dst       The width (in bytes) of the destination array.
                                            This corresponds to the padded number of IVPs to be solved.
 * \param[in]               src             The source pointer
 * \param[in]               pitch_src       The width (in bytes)  of the source array.
                                            This corresponds to the (non-padded) number of IVPs read by read_initial_conditions
 * \param[in]               offset          The offset within the source array (IVP index) to copy from.
                                            This is useful in the case (for large models) where the solver and state vector memory will not fit in device memory
                                            and the integration must be split into multiple kernel calls.
 * \param[in]               width           The size (in bytes) of memory to copy for each entry in the state vector
 * \param[in]               height          The number of entries in the state vector
 */
inline void memcpy2D_in(double* dst, const int pitch_dst, double const * src, const int pitch_src,
                                     const int offset, const size_t width, const int height) {
    for (int i = 0; i < height; ++i)
    {
        memcpy(dst, &src[offset], width);
        dst += pitch_dst;
        src += pitch_src;
    }
}

/**
 * \brief A convienience method to copy memory between host pointers of different pitches, widths and heights.
 *        Enables easier use of CUDA's cudaMemcpy2D functions.
 *
 * \param[out]              dst             The destination array
 * \param[in]               pitch_dst       The width (in bytes)  of the source array.
                                            This corresponds to the (non-padded) number of IVPs read by read_initial_conditions
 * \param[in]               src             The source pointer
 * \param[in]               pitch_src       The width (in bytes) of the destination array.
                                            This corresponds to the padded number of IVPs to be solved.
 * \param[in]               offset          The offset within the destination array (IVP index) to copy to.
                                            This is useful in the case (for large models) where the solver and state vector memory will not fit in device memory
                                            and the integration must be split into multiple kernel calls.
 * \param[in]               width           The size (in bytes) of memory to copy for each entry in the state vector
 * \param[in]               height          The number of entries in the state vector
 */
inline void memcpy2D_out(double* dst, const int pitch_dst, double const * src, const int pitch_src,
                                      const int offset, const size_t width, const int height) {
    for (int i = 0; i < height; ++i)
    {
        memcpy(&dst[offset], src, width);
        dst += pitch_dst;
        src += pitch_src;
    }
}




/**
 * \brief Cleans up the solver
 */
void accelerInt_cleanup() {
    free_gpu_memory(&host_mech, &device_mech);
    cleanup_solver(&host_solver, &device_solver);
    free(y_temp);
    free(host_mech);
    free(host_solver);
    free(result_flag);
    cudaErrorCheck( cudaDeviceReset() );
}




#ifdef GENERATE_DOCS
}
#endif
