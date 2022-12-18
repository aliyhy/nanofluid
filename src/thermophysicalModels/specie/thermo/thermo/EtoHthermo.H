
#include "fvMesh.H"

inline scalar Cp
(
    const fvMesh& mesh, const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar cv = Cv(mesh, p, T);
    volatile const scalar cpmcv = EquationOfState::CpMCv(mesh, p, T);
    return cv + cpmcv;
    #else
    return Cv(mesh, p, T) + EquationOfState::CpMCv(mesh, p, T);
    #endif
}

inline scalar Hs
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar es = Es(mesh, p, T);
    volatile const scalar rho = EquationOfState::rho(mesh, p, T);
    return es + p/rho;
    #else
    return Es(mesh, p, T) + p/EquationOfState::rho(mesh, p, T);
    #endif
}

inline scalar Ha
(
    const fvMesh& mesh,
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar ea = Ea(mesh, p, T);
    volatile const scalar rho = EquationOfState::rho(mesh, p, T);
    return ea + p/rho;
    #else
    return Ea(mesh, p, T) + p/EquationOfState::rho(mesh, p, T);
    #endif
}
