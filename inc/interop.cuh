#if CUDA_VERSION >= 11000

#define cusparseScsrmv(handle, op, rows, cols, nnz, alpha, descr, dval, drow, dcol, x, beta, y)        \
    {                                                                                                  \
        cusparseSpMatDescr_t descrA;                                                                   \
        cusparseDnVecDescr_t descrX, descrY;                                                           \
        cusparseCreateCsr(&descrA, rows, cols, nnz,                                                    \
                          (void *)drow, (void *)dcol, (void *)dval,                                    \
                          CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,                                      \
                          cusparseGetMatIndexBase(descr), CUDA_R_32F); \
        cusparseCreateDnVec(&descrX, cols, x, CUDA_R_32F);                                             \
        cusparseCreateDnVec(&descrY, rows, y, CUDA_R_32F);                                             \
                                                                                                       \
        size_t bufsize;                                                                                \
        void *buf;                                                                                     \
        cusparseSpMV_bufferSize(handle, op,                                                            \
                                (void *)alpha, descrA, descrX, (void *)beta,                           \
                                descrY, CUDA_R_32F, CUSPARSE_CSRMV_ALG1, &bufsize);                    \
        if (bufsize > 0)                                                                               \
            cudaMalloc(&buf, bufsize);                                                                 \
        cusparseSpMV(handle, op,                                                                       \
                     (void *)alpha, descrA, descrX, (void *)beta,                                      \
                     descrY, CUDA_R_32F, CUSPARSE_CSRMV_ALG1, buf);                                    \
        if (bufsize > 0)                                                                               \
            cudaFree(buf);                                                                             \
    }
#endif