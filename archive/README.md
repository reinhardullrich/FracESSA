# Archived source

This directory preserves source that is no longer compiled or used by FracESSA.

`fraction_free_ldlt.hpp` is the former general symmetric fraction-free LDLT factorization. It became unused when the exact
copositivity path stopped performing a separate positive-definiteness factorization. The active candidate solver continues to use
the independent KKT-specialized implementation in `cpp/include/linalg/fraction_free_ldlt_kkt.hpp`.
