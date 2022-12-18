#include "fvMesh.H"

inline scalar Cv
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar cp = Cp(mesh, p, T);
    volatile const scalar cpmcv = EquationOfState::CpMCv(mesh, p, T);
    return cp - cpmcv;
    #else
    return Cp(mesh, p, T) - EquationOfState::CpMCv(mesh, p, T);
    #endif
}

inline scalar Es
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar hs = Hs(mesh, p, T);
    volatile const scalar rho = EquationOfState::rho(mesh, p, T);
    return hs - p/rho;
    #else
    return Hs(mesh, p, T) - p/EquationOfState::rho(mesh, p, T);
    #endif
}

inline scalar Ea
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar ha = Ha(mesh, p, T);
    volatile const scalar rho = EquationOfState::rho(mesh, p, T);
    return ha - p/rho;
    #else
    return Ha(mesh, p, T) - p/EquationOfState::rho(mesh, p, T);
    #endif
}
