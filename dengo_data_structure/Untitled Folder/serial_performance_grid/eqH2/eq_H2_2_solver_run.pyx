cimport numpy as np
import numpy as np
import time
from libc.stdlib cimport malloc, free

cdef extern from "alloca.h":
    void *alloca(int)

##############################################################################

cdef extern from "sundials/sundials_types.h":
    
    ctypedef double realtype
    ctypedef int booleantype
    
    enum: TRUE
    enum: FALSE

cdef extern from "sundials/sundials_nvector.h":

    cdef struct _generic_N_Vector_Ops
    ctypedef _generic_N_Vector_Ops *N_Vector_Ops
    
    cdef struct _generic_N_Vector 
    ctypedef _generic_N_Vector *N_Vector

    ctypedef N_Vector *N_Vector_S

    cdef struct _generic_N_Vector_Ops:
        N_Vector    (*nvclone)(N_Vector)
        N_Vector    (*nvcloneempty)(N_Vector)
        void        (*nvdestroy)(N_Vector)
        void        (*nvspace)(N_Vector, long *, long *)
        realtype*   (*nvgetarraypointer)(N_Vector)
        void        (*nvsetarraypointer)(realtype *, N_Vector)
        void        (*nvlinearsum)(realtype, N_Vector, realtype, N_Vector, N_Vector)
        void        (*nvconst)(realtype, N_Vector)
        void        (*nvprod)(N_Vector, N_Vector, N_Vector)
        void        (*nvdiv)(N_Vector, N_Vector, N_Vector)
        void        (*nvscale)(realtype, N_Vector, N_Vector)
        void        (*nvabs)(N_Vector, N_Vector)
        void        (*nvinv)(N_Vector, N_Vector)
        void        (*nvaddconst)(N_Vector, realtype, N_Vector)
        realtype    (*nvdotprod)(N_Vector, N_Vector)
        realtype    (*nvmaxnorm)(N_Vector)
        realtype    (*nvwrmsnorm)(N_Vector, N_Vector)
        realtype    (*nvwrmsnormmask)(N_Vector, N_Vector, N_Vector)
        realtype    (*nvmin)(N_Vector)
        realtype    (*nvwl2norm)(N_Vector, N_Vector)
        realtype    (*nvl1norm)(N_Vector)
        void        (*nvcompare)(realtype, N_Vector, N_Vector)
        booleantype (*nvinvtest)(N_Vector, N_Vector)
        booleantype (*nvconstrmask)(N_Vector, N_Vector, N_Vector)
        realtype    (*nvminquotient)(N_Vector, N_Vector)
    
    cdef struct _generic_N_Vector:
        void *content
        _generic_N_Vector_Ops *ops   

cdef extern from "sundials/sundials_direct.h":        
#    /*
#     * =================================================================
#     *                C O N S T A N T S
#     * =================================================================
#     */
#    
#    /*
#     *  SUNDIALS_DENSE: dense matrix
#     *  SUNDIALS_BAND:  banded matrix
#     */
#    
    enum: SUNDIALS_DENSE
    enum: SUNDIALS_BAND
#    
#    /*
#     * ==================================================================
#     * Type definitions
#     * ==================================================================
#     */
#    
#    /*
#     * -----------------------------------------------------------------
#     * Type : DlsMat
#     * -----------------------------------------------------------------
#     * The type DlsMat is defined to be a pointer to a structure
#     * with various sizes, a data field, and an array of pointers to
#     * the columns which defines a dense or band matrix for use in 
#     * direct linear solvers. The M and N fields indicates the number 
#     * of rows and columns, respectively. The data field is a one 
#     * dimensional array used for component storage. The cols field 
#     * stores the pointers in data for the beginning of each column.
#     * -----------------------------------------------------------------
#     * For DENSE matrices, the relevant fields in DlsMat are:
#     *    type  = SUNDIALS_DENSE
#     *    M     - number of rows
#     *    N     - number of columns
#     *    ldim  - leading dimension (ldim >= M)
#     *    data  - pointer to a contiguous block of realtype variables
#     *    ldata - length of the data array =ldim*N
#     *    cols  - array of pointers. cols[j] points to the first element 
#     *            of the j-th column of the matrix in the array data.
#     *
#     * The elements of a dense matrix are stored columnwise (i.e columns 
#     * are stored one on top of the other in memory). 
#     * If A is of type DlsMat, then the (i,j)th element of A (with 
#     * 0 <= i < M and 0 <= j < N) is given by (A->data)[j*n+i]. 
#     *
#     * The DENSE_COL and DENSE_ELEM macros below allow a user to access 
#     * efficiently individual matrix elements without writing out explicit 
#     * data structure references and without knowing too much about the 
#     * underlying element storage. The only storage assumption needed is 
#     * that elements are stored columnwise and that a pointer to the 
#     * jth column of elements can be obtained via the DENSE_COL macro.
#     * -----------------------------------------------------------------
#     * For BAND matrices, the relevant fields in DlsMat are:
#     *    type  = SUNDIALS_BAND
#     *    M     - number of rows
#     *    N     - number of columns
#     *    mu    - upper bandwidth, 0 <= mu <= min(M,N)
#     *    ml    - lower bandwidth, 0 <= ml <= min(M,N)
#     *    s_mu  - storage upper bandwidth, mu <= s_mu <= N-1.
#     *            The dgbtrf routine writes the LU factors into the storage 
#     *            for A. The upper triangular factor U, however, may have 
#     *            an upper bandwidth as big as MIN(N-1,mu+ml) because of 
#     *            partial pivoting. The s_mu field holds the upper 
#     *            bandwidth allocated for A.
#     *    ldim  - leading dimension (ldim >= s_mu)
#     *    data  - pointer to a contiguous block of realtype variables
#     *    ldata - length of the data array =ldim*(s_mu+ml+1)
#     *    cols  - array of pointers. cols[j] points to the first element 
#     *            of the j-th column of the matrix in the array data.
#     *
#     * The BAND_COL, BAND_COL_ELEM, and BAND_ELEM macros below allow a 
#     * user to access individual matrix elements without writing out 
#     * explicit data structure references and without knowing too much 
#     * about the underlying element storage. The only storage assumption 
#     * needed is that elements are stored columnwise and that a pointer 
#     * into the jth column of elements can be obtained via the BAND_COL 
#     * macro. The BAND_COL_ELEM macro selects an element from a column
#     * which has already been isolated via BAND_COL. The macro 
#     * BAND_COL_ELEM allows the user to avoid the translation 
#     * from the matrix location (i,j) to the index in the array returned 
#     * by BAND_COL at which the (i,j)th element is stored. 
#     * -----------------------------------------------------------------
#     */
#    
    cdef struct _DlsMat:
        int type
        long int M
        long int N
        long int ldim
        long int mu
        long int ml
        long int s_mu
        realtype *data
        long int ldata
        realtype **cols
        
    ctypedef _DlsMat *DlsMat
    
#    /*
#     * ==================================================================
#     * Data accessor macros
#     * ==================================================================
#     */
#    
#    /*
#     * -----------------------------------------------------------------
#     * DENSE_COL and DENSE_ELEM
#     * -----------------------------------------------------------------
#     *
#     * DENSE_COL(A,j) references the jth column of the M-by-N dense
#     * matrix A, 0 <= j < N. The type of the expression DENSE_COL(A,j) 
#     * is (realtype *). After the assignment in the usage above, col_j 
#     * may be treated as an array indexed from 0 to M-1. The (i,j)-th 
#     * element of A is thus referenced by col_j[i].
#     *
#     * DENSE_ELEM(A,i,j) references the (i,j)th element of the dense 
#     * M-by-N matrix A, 0 <= i < M ; 0 <= j < N.
#     *
#     * -----------------------------------------------------------------
#     */
#    
#    #define DENSE_COL(A,j) ((A->cols)[j])
    void DENSE_COL(DlsMat A, long int j)
#    #define DENSE_ELEM(A,i,j) ((A->cols)[j][i])
    void DENSE_ELEM(DlsMat A, long int i, long int j)    
#    
#    /*
#     * -----------------------------------------------------------------
#     * BAND_COL, BAND_COL_ELEM, and BAND_ELEM
#     * -----------------------------------------------------------------
#     *  
#     * BAND_COL(A,j) references the diagonal element of the jth column 
#     * of the N by N band matrix A, 0 <= j <= N-1. The type of the 
#     * expression BAND_COL(A,j) is realtype *. The pointer returned by 
#     * the call BAND_COL(A,j) can be treated as an array which is 
#     * indexed from -(A->mu) to (A->ml).
#     * 
#     * BAND_COL_ELEM references the (i,j)th entry of the band matrix A 
#     * when used in conjunction with BAND_COL. The index (i,j) should 
#     * satisfy j-(A->mu) <= i <= j+(A->ml).
#     *
#     * BAND_ELEM(A,i,j) references the (i,j)th element of the M-by-N 
#     * band matrix A, where 0 <= i,j <= N-1. The location (i,j) should 
#     * further satisfy j-(A->mu) <= i <= j+(A->ml). 
#     *
#     * -----------------------------------------------------------------
#     */
#     
#    #define BAND_COL(A,j) (((A->cols)[j])+(A->s_mu))
#    #define BAND_COL_ELEM(col_j,i,j) (col_j[(i)-(j)])
#    #define BAND_ELEM(A,i,j) ((A->cols)[j][(i)-(j)+(A->s_mu)])    
        
cdef extern from "sundials/sundials_iterative.h":
    #/*
    # * -----------------------------------------------------------------
    # * enum : types of preconditioning                                
    # * -----------------------------------------------------------------
    # * PREC_NONE  : The iterative linear solver should not use             
    # *              preconditioning.                                       
    # *                                                                
    # * PREC_LEFT  : The iterative linear solver uses preconditioning on    
    # *              the left only.                                         
    # *                                                                
    # * PREC_RIGHT : The iterative linear solver uses preconditioning on    
    # *              the right only.                                        
    # *                                                                
    # * PREC_BOTH  : The iterative linear solver uses preconditioning on    
    # *              both the left and the right.                           
    # * -----------------------------------------------------------------
    # */

    enum: PREC_NONE
    enum: PREC_LEFT
    enum: PREC_RIGHT
    enum: PREC_BOTH 
    
    #/*
    # * -----------------------------------------------------------------
    # * enum : types of Gram-Schmidt routines                          
    # * -----------------------------------------------------------------
    # * MODIFIED_GS  : The iterative solver uses the modified          
    # *                Gram-Schmidt routine ModifiedGS listed in this  
    # *                file.                                           
    # *                                                                
    # * CLASSICAL_GS : The iterative solver uses the classical         
    # *                Gram-Schmidt routine ClassicalGS listed in this 
    # *                file.                                           
    # * -----------------------------------------------------------------
    # */
    #
    enum: MODIFIED_GS
    enum: CLASSICAL_GS
    
    
cdef extern from "sundials/sundials_dense.h":
    #/*
    # * -----------------------------------------------------------------
    # * Functions: DenseGETRF and DenseGETRS
    # * -----------------------------------------------------------------
    # * DenseGETRF performs the LU factorization of the M by N dense
    # * matrix A. This is done using standard Gaussian elimination
    # * with partial (row) pivoting. Note that this applies only
    # * to matrices with M >= N and full column rank.
    # *
    # * A successful LU factorization leaves the matrix A and the
    # * pivot array p with the following information:
    # *
    # * (1) p[k] contains the row number of the pivot element chosen
    # *     at the beginning of elimination step k, k=0, 1, ..., N-1.
    # *
    # * (2) If the unique LU factorization of A is given by PA = LU,
    # *     where P is a permutation matrix, L is a lower trapezoidal
    # *     matrix with all 1's on the diagonal, and U is an upper
    # *     triangular matrix, then the upper triangular part of A
    # *     (including its diagonal) contains U and the strictly lower
    # *     trapezoidal part of A contains the multipliers, I-L.
    # *
    # * For square matrices (M=N), L is unit lower triangular.
    # *
    # * DenseGETRF returns 0 if successful. Otherwise it encountered
    # * a zero diagonal element during the factorization. In this case
    # * it returns the column index (numbered from one) at which
    # * it encountered the zero.
    # *
    # * DenseGETRS solves the N-dimensional system A x = b using
    # * the LU factorization in A and the pivot information in p
    # * computed in DenseGETRF. The solution x is returned in b. This
    # * routine cannot fail if the corresponding call to DenseGETRF
    # * did not fail.
    # * DenseGETRS does NOT check for a square matrix!
    # *
    # * -----------------------------------------------------------------
    # * DenseGETRF and DenseGETRS are simply wrappers around denseGETRF
    # * and denseGETRS, respectively, which perform all the work by
    # * directly accessing the data in the DlsMat A (i.e. the field cols)
    # * -----------------------------------------------------------------
    # */
    long int denseGETRF(realtype **a, long int m, long int n, long int *p)
    void denseGETRS(realtype **a, long int n, long int *p, realtype *b) 
        
    
    #/*
    # * -----------------------------------------------------------------
    # * Functions : DensePOTRF and DensePOTRS
    # * -----------------------------------------------------------------
    # * DensePOTRF computes the Cholesky factorization of a real symmetric
    # * positive definite matrix A.
    # * -----------------------------------------------------------------
    # * DensePOTRS solves a system of linear equations A*X = B with a 
    # * symmetric positive definite matrix A using the Cholesky factorization
    # * A = L*L**T computed by DensePOTRF.
    # *
    # * -----------------------------------------------------------------
    # * DensePOTRF and DensePOTRS are simply wrappers around densePOTRF
    # * and densePOTRS, respectively, which perform all the work by
    # * directly accessing the data in the DlsMat A (i.e. the field cols)
    # * -----------------------------------------------------------------
    # */
    #
    #SUNDIALS_EXPORT long int DensePOTRF(DlsMat A);
    #SUNDIALS_EXPORT void DensePOTRS(DlsMat A, realtype *b);
    
    long int densePOTRF(realtype **a, long int m)
    void densePOTRS(realtype **a, long int m, realtype *b)
    
    #/*
    # * -----------------------------------------------------------------
    # * Functions : DenseGEQRF and DenseORMQR
    # * -----------------------------------------------------------------
    # * DenseGEQRF computes a QR factorization of a real M-by-N matrix A:
    # * A = Q * R (with M>= N).
    # * 
    # * DenseGEQRF requires a temporary work vector wrk of length M.
    # * -----------------------------------------------------------------
    # * DenseORMQR computes the product w = Q * v where Q is a real 
    # * orthogonal matrix defined as the product of k elementary reflectors
    # *
    # *        Q = H(1) H(2) . . . H(k)
    # *
    # * as returned by DenseGEQRF. Q is an M-by-N matrix, v is a vector
    # * of length N and w is a vector of length M (with M>=N).
    # *
    # * DenseORMQR requires a temporary work vector wrk of length M.
    # *
    # * -----------------------------------------------------------------
    # * DenseGEQRF and DenseORMQR are simply wrappers around denseGEQRF
    # * and denseORMQR, respectively, which perform all the work by
    # * directly accessing the data in the DlsMat A (i.e. the field cols)
    # * -----------------------------------------------------------------
    # */
    #
    #SUNDIALS_EXPORT int DenseGEQRF(DlsMat A, realtype *beta, realtype *wrk);
    #SUNDIALS_EXPORT int DenseORMQR(DlsMat A, realtype *beta, realtype *vn, realtype *vm, 
    #			       realtype *wrk);
    #
    int denseGEQRF(realtype **a, long int m, long int n, realtype *beta, realtype *v)
    int denseORMQR(realtype **a, long int m, long int n, realtype *beta,
    			       realtype *v, realtype *w, realtype *wrk)
    #
    #/*
    # * -----------------------------------------------------------------
    # * Function : DenseCopy
    # * -----------------------------------------------------------------
    # * DenseCopy copies the contents of the M-by-N matrix A into the
    # * M-by-N matrix B.
    # * 
    # * DenseCopy is a wrapper around denseCopy which accesses the data
    # * in the DlsMat A and B (i.e. the fields cols)
    # * -----------------------------------------------------------------
    # */
    
    #SUNDIALS_EXPORT void DenseCopy(DlsMat A, DlsMat B);
    void denseCopy(realtype **a, realtype **b, long int m, long int n)
    
    #/*
    # * -----------------------------------------------------------------
    # * Function: DenseScale
    # * -----------------------------------------------------------------
    # * DenseScale scales the elements of the M-by-N matrix A by the
    # * constant c and stores the result back in A.
    # *
    # * DenseScale is a wrapper around denseScale which performs the actual
    # * scaling by accessing the data in the DlsMat A (i.e. the field
    # * cols).
    # * -----------------------------------------------------------------
    # */
    #
    #SUNDIALS_EXPORT void DenseScale(realtype c, DlsMat A);
    void denseScale(realtype c, realtype **a, long int m, long int n)
    
    
    #/*
    # * -----------------------------------------------------------------
    # * Function: denseAddIdentity
    # * -----------------------------------------------------------------
    # * denseAddIdentity adds the identity matrix to the n-by-n matrix
    # * stored in the realtype** arrays.
    # * -----------------------------------------------------------------
    # */
    #
    void denseAddIdentity(realtype **a, long int n)
    #
    #    


##############################################################################


# NSPECIES here is N in the .C.template file
DEF NSPECIES = 10
DEF MAX_NCELLS=1024

cdef extern from "eq_H2_2_solver.h":
    cdef int _MAX_NCELLS  "MAX_NCELLS"
    cdef int _NSPECIES "NSPECIES"

    ctypedef struct eq_H2_2_data:
        double dbin
        double idbin
        double bounds[2]
        int nbins

        double d_zbin
        double id_zbin
        double z_bounds[2]
        int n_zbins

        double current_z
        double zdef
        double dz

        double Ts[MAX_NCELLS]
        double Tdef[MAX_NCELLS]
        double dT[MAX_NCELLS]
        double logTs[MAX_NCELLS]
        double dTs_ge[MAX_NCELLS]
        double r_k01[1024]
        double rs_k01[MAX_NCELLS]
        double drs_k01[MAX_NCELLS]
        double r_k02[1024]
        double rs_k02[MAX_NCELLS]
        double drs_k02[MAX_NCELLS]
        double r_k03[1024]
        double rs_k03[MAX_NCELLS]
        double drs_k03[MAX_NCELLS]
        double r_k04[1024]
        double rs_k04[MAX_NCELLS]
        double drs_k04[MAX_NCELLS]
        double r_k05[1024]
        double rs_k05[MAX_NCELLS]
        double drs_k05[MAX_NCELLS]
        double r_k06[1024]
        double rs_k06[MAX_NCELLS]
        double drs_k06[MAX_NCELLS]
        double r_k07[1024]
        double rs_k07[MAX_NCELLS]
        double drs_k07[MAX_NCELLS]
        double r_k08[1024]
        double rs_k08[MAX_NCELLS]
        double drs_k08[MAX_NCELLS]
        double r_k09[1024]
        double rs_k09[MAX_NCELLS]
        double drs_k09[MAX_NCELLS]
        double r_k10[1024]
        double rs_k10[MAX_NCELLS]
        double drs_k10[MAX_NCELLS]
        double r_k11[1024]
        double rs_k11[MAX_NCELLS]
        double drs_k11[MAX_NCELLS]
        double r_k12[1024]
        double rs_k12[MAX_NCELLS]
        double drs_k12[MAX_NCELLS]
        double r_k13[1024]
        double rs_k13[MAX_NCELLS]
        double drs_k13[MAX_NCELLS]
        double r_k14[1024]
        double rs_k14[MAX_NCELLS]
        double drs_k14[MAX_NCELLS]
        double r_k15[1024]
        double rs_k15[MAX_NCELLS]
        double drs_k15[MAX_NCELLS]
        double r_k16[1024]
        double rs_k16[MAX_NCELLS]
        double drs_k16[MAX_NCELLS]
        double r_k17[1024]
        double rs_k17[MAX_NCELLS]
        double drs_k17[MAX_NCELLS]
        double r_k18[1024]
        double rs_k18[MAX_NCELLS]
        double drs_k18[MAX_NCELLS]
        double r_k19[1024]
        double rs_k19[MAX_NCELLS]
        double drs_k19[MAX_NCELLS]
        double r_k21[1024]
        double rs_k21[MAX_NCELLS]
        double drs_k21[MAX_NCELLS]
        double r_k22[1024]
        double rs_k22[MAX_NCELLS]
        double drs_k22[MAX_NCELLS]
        double c_brem_brem[1024]
        double cs_brem_brem[MAX_NCELLS]
        double dcs_brem_brem[MAX_NCELLS]
        
        double c_ceHeI_ceHeI[1024]
        double cs_ceHeI_ceHeI[MAX_NCELLS]
        double dcs_ceHeI_ceHeI[MAX_NCELLS]
        
        double c_ceHeII_ceHeII[1024]
        double cs_ceHeII_ceHeII[MAX_NCELLS]
        double dcs_ceHeII_ceHeII[MAX_NCELLS]
        
        double c_ceHI_ceHI[1024]
        double cs_ceHI_ceHI[MAX_NCELLS]
        double dcs_ceHI_ceHI[MAX_NCELLS]
        
        double c_ciHeI_ciHeI[1024]
        double cs_ciHeI_ciHeI[MAX_NCELLS]
        double dcs_ciHeI_ciHeI[MAX_NCELLS]
        
        double c_ciHeII_ciHeII[1024]
        double cs_ciHeII_ciHeII[MAX_NCELLS]
        double dcs_ciHeII_ciHeII[MAX_NCELLS]
        
        double c_ciHeIS_ciHeIS[1024]
        double cs_ciHeIS_ciHeIS[MAX_NCELLS]
        double dcs_ciHeIS_ciHeIS[MAX_NCELLS]
        
        double c_ciHI_ciHI[1024]
        double cs_ciHI_ciHI[MAX_NCELLS]
        double dcs_ciHI_ciHI[MAX_NCELLS]
        
        double c_compton_comp_[1024]
        double cs_compton_comp_[MAX_NCELLS]
        double dcs_compton_comp_[MAX_NCELLS]
        
        double c_gammah_gammah[1024]
        double cs_gammah_gammah[MAX_NCELLS]
        double dcs_gammah_gammah[MAX_NCELLS]
        
        double c_gloverabel08_gael[1024]
        double cs_gloverabel08_gael[MAX_NCELLS]
        double dcs_gloverabel08_gael[MAX_NCELLS]
        double c_gloverabel08_gaH2[1024]
        double cs_gloverabel08_gaH2[MAX_NCELLS]
        double dcs_gloverabel08_gaH2[MAX_NCELLS]
        double c_gloverabel08_gaHe[1024]
        double cs_gloverabel08_gaHe[MAX_NCELLS]
        double dcs_gloverabel08_gaHe[MAX_NCELLS]
        double c_gloverabel08_gaHI[1024]
        double cs_gloverabel08_gaHI[MAX_NCELLS]
        double dcs_gloverabel08_gaHI[MAX_NCELLS]
        double c_gloverabel08_gaHp[1024]
        double cs_gloverabel08_gaHp[MAX_NCELLS]
        double dcs_gloverabel08_gaHp[MAX_NCELLS]
        double c_gloverabel08_gphdl[1024]
        double cs_gloverabel08_gphdl[MAX_NCELLS]
        double dcs_gloverabel08_gphdl[MAX_NCELLS]
        double c_gloverabel08_gpldl[1024]
        double cs_gloverabel08_gpldl[MAX_NCELLS]
        double dcs_gloverabel08_gpldl[MAX_NCELLS]
        double c_gloverabel08_h2lte[1024]
        double cs_gloverabel08_h2lte[MAX_NCELLS]
        double dcs_gloverabel08_h2lte[MAX_NCELLS]
        
        double c_h2formation_h2mcool[1024]
        double cs_h2formation_h2mcool[MAX_NCELLS]
        double dcs_h2formation_h2mcool[MAX_NCELLS]
        double c_h2formation_h2mheat[1024]
        double cs_h2formation_h2mheat[MAX_NCELLS]
        double dcs_h2formation_h2mheat[MAX_NCELLS]
        double c_h2formation_ncrd1[1024]
        double cs_h2formation_ncrd1[MAX_NCELLS]
        double dcs_h2formation_ncrd1[MAX_NCELLS]
        double c_h2formation_ncrd2[1024]
        double cs_h2formation_ncrd2[MAX_NCELLS]
        double dcs_h2formation_ncrd2[MAX_NCELLS]
        double c_h2formation_ncrn[1024]
        double cs_h2formation_ncrn[MAX_NCELLS]
        double dcs_h2formation_ncrn[MAX_NCELLS]
        
        double c_reHeII1_reHeII1[1024]
        double cs_reHeII1_reHeII1[MAX_NCELLS]
        double dcs_reHeII1_reHeII1[MAX_NCELLS]
        
        double c_reHeII2_reHeII2[1024]
        double cs_reHeII2_reHeII2[MAX_NCELLS]
        double dcs_reHeII2_reHeII2[MAX_NCELLS]
        
        double c_reHeIII_reHeIII[1024]
        double cs_reHeIII_reHeIII[MAX_NCELLS]
        double dcs_reHeIII_reHeIII[MAX_NCELLS]
        
        double c_reHII_reHII[1024]
        double cs_reHII_reHII[MAX_NCELLS]
        double dcs_reHII_reHII[MAX_NCELLS]
        
        int bin_id[MAX_NCELLS]
        int ncells
    
    # Declare ctype RHS and Jacobian
    ctypedef int(*rhs_f)( realtype, N_Vector , N_Vector , void * )
    ctypedef int(*jac_f)( long int, realtype, N_Vector , N_Vector , DlsMat , void *, N_Vector, N_Vector, N_Vector)

    int eq_H2_2_main(int argc, char **argv)
    eq_H2_2_data *eq_H2_2_setup_data(int *NumberOfFields,
            char ***FieldNames)
    void eq_H2_2_read_rate_tables(eq_H2_2_data*)
    void eq_H2_2_read_cooling_tables(eq_H2_2_data*)
    void eq_H2_2_read_gamma(eq_H2_2_data*)

    double dengo_evolve_eq_H2_2 (double dtf, double &dt, double z,
                                         double *input, double *rtol,
                                         double *atol, int dims,
                                         eq_H2_2_data *data, double *temp)

    # Declare the Jacobian and RHS function
    
    int calculate_jacobian_eq_H2_2(long int N, realtype t,
                   N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    int calculate_rhs_eq_H2_2(realtype t, N_Vector y, 
                    N_Vector ydot, void *user_data);
    
    void eq_H2_2_calculate_temperature(eq_H2_2_data *data, double *input, int nstrip, int nchem);
   
       
    int ensure_electron_consistency(double *input, int nstrip, int nchem);
    

def main_run_eq_H2_2():
    t1 = time.time()
    eq_H2_2_main(0, NULL)
    t2 = time.time()
    print "Total elapsed time: %0.3e" % (t2-t1)

def run_eq_H2_2(ics, double tf, int niter = 10000,
                        int intermediate = 1, z = -1.0, dtarr = -999999):
    assert(_MAX_NCELLS == MAX_NCELLS)
    assert(_NSPECIES == NSPECIES)
    cdef np.ndarray[np.float64_t, ndim=1] H2_1_arr = ics["H2_1"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] H2_1_int
    cdef np.ndarray[np.float64_t, ndim=1] H2_2_arr = ics["H2_2"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] H2_2_int
    cdef np.ndarray[np.float64_t, ndim=1] H_1_arr = ics["H_1"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] H_1_int
    cdef np.ndarray[np.float64_t, ndim=1] H_2_arr = ics["H_2"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] H_2_int
    cdef np.ndarray[np.float64_t, ndim=1] H_m0_arr = ics["H_m0"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] H_m0_int
    cdef np.ndarray[np.float64_t, ndim=1] He_1_arr = ics["He_1"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] He_1_int
    cdef np.ndarray[np.float64_t, ndim=1] He_2_arr = ics["He_2"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] He_2_int
    cdef np.ndarray[np.float64_t, ndim=1] He_3_arr = ics["He_3"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] He_3_int
    cdef np.ndarray[np.float64_t, ndim=1] de_arr = ics["de"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] de_int
    cdef np.ndarray[np.float64_t, ndim=1] ge_arr = ics["ge"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] ge_int
    cdef np.ndarray[np.float64_t, ndim=1] result_int
    cdef np.ndarray[np.float64_t, ndim=2] temp_int
    cdef np.ndarray[np.float64_t, ndim=1] t_int
    cdef np.ndarray[np.float64_t, ndim=1] dt_int

    cdef int i, j, k, iter
    cdef int dims = ge_arr.shape[0]
    cdef int NTOT = NSPECIES * dims
    cdef double *input = <double *> alloca(NTOT * sizeof(double))
    cdef double *prev = <double *> alloca(NTOT * sizeof(double))
    cdef double *atol = <double *> alloca(NTOT * sizeof(double))
    cdef double *rtol = <double *> alloca(NTOT * sizeof(double))
    cdef double *scale = <double *> alloca(NTOT * sizeof(double))
    cdef double *dt   = <double *> alloca( sizeof(double))
    cdef double *ones = <double *> alloca(NTOT * sizeof(double))

    for i in range(NTOT):
        ones[i] = 1.0

    dt[0] = tf/float(niter)

    if intermediate == 1:
        H2_1_int = np.zeros((dims, niter), "float64")
        H2_2_int = np.zeros((dims, niter), "float64")
        H_1_int = np.zeros((dims, niter), "float64")
        H_2_int = np.zeros((dims, niter), "float64")
        H_m0_int = np.zeros((dims, niter), "float64")
        He_1_int = np.zeros((dims, niter), "float64")
        He_2_int = np.zeros((dims, niter), "float64")
        He_3_int = np.zeros((dims, niter), "float64")
        de_int = np.zeros((dims, niter), "float64")
        ge_int = np.zeros((dims, niter), "float64")
        temp_int = np.zeros((dims, niter), "float64")
        result_int = np.zeros(niter, "float64")
        t_int = np.zeros(niter, "float64")
        dt_int = np.zeros(niter, "float64")

    j = 0
    for i in range(dims):
        input[j] = prev[j] = H2_1_arr[i] 
        atol[j] = input[j] * 1e-6
        rtol[j] = 1e-06
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = H2_2_arr[i] 
        atol[j] = input[j] * 1e-6
        rtol[j] = 1e-06
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = H_1_arr[i] 
        atol[j] = input[j] * 1e-6
        rtol[j] = 1e-06
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = H_2_arr[i] 
        atol[j] = input[j] * 1e-6
        rtol[j] = 1e-06
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = H_m0_arr[i] 
        atol[j] = input[j] * 1e-6
        rtol[j] = 1e-06
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = He_1_arr[i] 
        atol[j] = input[j] * 1e-6
        rtol[j] = 1e-06
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = He_2_arr[i] 
        atol[j] = input[j] * 1e-6
        rtol[j] = 1e-06
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = He_3_arr[i] 
        atol[j] = input[j] * 1e-6
        rtol[j] = 1e-06
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = de_arr[i] 
        atol[j] = input[j] * 1e-6
        rtol[j] = 1e-06
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = ge_arr[i] 
        atol[j] = input[j] * 1e-6
        rtol[j] = 1e-06
        scale[j] = input[j]
        j += 1
        

    cdef eq_H2_2_data *data = eq_H2_2_setup_data(NULL, NULL)
    cdef rhs_f f = calculate_rhs_eq_H2_2
    cdef jac_f jf = calculate_jacobian_eq_H2_2

    cdef double ttot = 0.0
    cdef int status
    # Allocate some temporary data
    # Now we manually evolve
    #ttot = dengo_evolve_eq_H2_2(tf, dt, input, rtol, atol, dims, data)
    data.current_z = z
    cdef double *t_now = <double *> malloc( sizeof(double) )

    
    cdef double *dt_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *success_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *ttot_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *temp = <double *> malloc(sizeof(double) * niter)

    cdef int niter_cvodes = niter
    
    cdef double dt_local;  
    dt0 = dt
    cdef double ttt = 0
    # Initialize initial temperature
    for i in range(dims):
        data.Ts[i] = ics['T'][i]
        print("initial  temperature: %.3E" %data.Ts[i])


    for iter in range(niter):
        dt_local = dengo_evolve_eq_H2_2( dt0[0], dt[0], z, input, rtol, atol, dims, data, temp  )
        if dt_local > 0:
            status = 0
            # print( "{} th iterations at time {}".format(iter, dt_local))
        
        j = 0; 

        for i in range(dims):
            H2_1_int[i, iter] = input[j]
            j += 1
            H2_2_int[i, iter] = input[j]
            j += 1
            H_1_int[i, iter] = input[j]
            j += 1
            H_2_int[i, iter] = input[j]
            j += 1
            H_m0_int[i, iter] = input[j]
            j += 1
            He_1_int[i, iter] = input[j]
            j += 1
            He_2_int[i, iter] = input[j]
            j += 1
            He_3_int[i, iter] = input[j]
            j += 1
            de_int[i, iter] = input[j]
            j += 1
            ge_int[i, iter] = input[j]
            j += 1
            temp_int[ i, iter ] = temp[i]

        if status == 0:
            result_int[iter] = 1
            ttot += dt_local
        elif status == 1:
            result_int[iter] = 0
            ttot += dt_local

        t_int[iter] = ttot
        dt_int[iter] = dt_local
        

        if status == 0:
            if iter % 100 == 0:
                print "Successful iteration[% 5i]: (%0.3e) %0.3e / %0.3e" % (iter, dt_local, ttot, tf)

            dt_local = 1.1*dt_local

            dt[0] = dt_local;
            if tf - ttot < dt_local:
                dt_local = tf - ttot
                dt[0] = dt_local;
        elif status == 1:
            dt[0] = dt_local/2.0;
            # copy_array(prev, input, NTOT)
            # Reset the scaling array to match the new values
            # copy_array(input, scale, NTOT)
            if dt[0] < 1e-50 * tf:
                print "dt too small (%0.3e / %0.3e) so breaking" % (dt[0], tf)
                break
            continue
        if ttot >= tf: break
    

    free(dt_arr)
    free(ttot_arr)
    free(success_arr)
    free(data)

    print "End in %s iterations: %0.5e / %0.5e (%0.5e)" % (iter + 1, ttot, tf, tf - ttot)

    rv, rv_t = {}, {}
    H2_1_arr = rv["H2_1"] = np.zeros(dims, "float64")
    H2_2_arr = rv["H2_2"] = np.zeros(dims, "float64")
    H_1_arr = rv["H_1"] = np.zeros(dims, "float64")
    H_2_arr = rv["H_2"] = np.zeros(dims, "float64")
    H_m0_arr = rv["H_m0"] = np.zeros(dims, "float64")
    He_1_arr = rv["He_1"] = np.zeros(dims, "float64")
    He_2_arr = rv["He_2"] = np.zeros(dims, "float64")
    He_3_arr = rv["He_3"] = np.zeros(dims, "float64")
    de_arr = rv["de"] = np.zeros(dims, "float64")
    ge_arr = rv["ge"] = np.zeros(dims, "float64")
    if intermediate:
        rv_t["H2_1"] = H2_1_int[:niter]
        rv_t["H2_2"] = H2_2_int[:niter]
        rv_t["H_1"] = H_1_int[:niter]
        rv_t["H_2"] = H_2_int[:niter]
        rv_t["H_m0"] = H_m0_int[:niter]
        rv_t["He_1"] = He_1_int[:niter]
        rv_t["He_2"] = He_2_int[:niter]
        rv_t["He_3"] = He_3_int[:niter]
        rv_t["de"] = de_int[:niter]
        rv_t["ge"] = ge_int[:niter]
        rv_t["successful"] = result_int.astype("bool")
        rv_t['T'] = temp_int
        rv_t['t'] = t_int
        rv_t['dt'] = dt_int

    j = 0
    for i in range(dims):
        H2_1_arr[i] = input[j] * 2.01588
        j += 1
        H2_2_arr[i] = input[j] * 2.01588
        j += 1
        H_1_arr[i] = input[j] * 1.00794
        j += 1
        H_2_arr[i] = input[j] * 1.00794
        j += 1
        H_m0_arr[i] = input[j] * 1.00794
        j += 1
        He_1_arr[i] = input[j] * 4.002602
        j += 1
        He_2_arr[i] = input[j] * 4.002602
        j += 1
        He_3_arr[i] = input[j] * 4.002602
        j += 1
        de_arr[i] = input[j] * 1.0
        j += 1
        ge_arr[i] = input[j] * 1.0
        j += 1
    return rv, rv_t

cdef copy_array(double *input, double *output, int dims):
    cdef int i
    for i in range(dims):
        output[i] = input[i]

cdef floor_values(double *input, int dims, double floor):
    cdef int i
    for i in range(dims):
        if input[i] < floor:
            input[i] = floor