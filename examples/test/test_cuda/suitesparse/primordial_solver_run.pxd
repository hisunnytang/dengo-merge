# this is the cython header files for the sundials cvode solver

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

