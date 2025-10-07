/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kOmegaLagEB.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

// Compute the turbulent viscosity (BDM: Appendix)
template<class BasicTurbulenceModel>
void kOmegaLagEB<BasicTurbulenceModel>::correctNut()
{   
    volScalarField magS =
        sqrt(2.0)*mag(symm(fvc::grad(this->U_)));

    this->nut_ = k_*
        min
        (
            phit_/omega_,
            0.7/sqrt(3.0)/
            max
            (
                magS,
                dimensionedScalar(pow(dimTime,-1),VSMALL)
            )
        );
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// Compute the turbulent length scale (BDM: Appendix)
template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaLagEB<BasicTurbulenceModel>::Ls() const
{
    return
        CL_*sqrt(
        max(pow3(k_)/sqr(betas_*k_*omega_), 
            dimensionedScalar(sqr(dimLength), Zero))
        + sqr(Ceta_)
        *sqrt(
            max(
                pow3(this->nu())/(betas_*k_*omega_),
                dimensionedScalar(pow(dimLength,4), Zero)
                )
            )
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaLagEB<BasicTurbulenceModel>::kOmegaLagEB
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    gamma1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            0.5
        )
    ),
    gamma2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.6
        )
    ),
    beta_
    (
        dimensionedScalar::getOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.075
        )
    ),
    betas_
    (
        dimensionedScalar::getOrAddToDict
        (
            "betas",
            this->coeffDict_,
            0.09
        )
    ),
    Cws_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cws",
            this->coeffDict_,
            0.05
        )
    ),
    Cp1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cp1",
            this->coeffDict_,
            0.4
        )
    ),
    Cp2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cp2",
            this->coeffDict_,
            0.46
        )
    ),
    C1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.7
        )
    ),
    C5_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C5",
            this->coeffDict_,
            0.2
        )
    ),
    C1s_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C1s",
            this->coeffDict_,
            0.9
        )
    ),
    C1Tilde_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C1Tilde",
            this->coeffDict_,
            1.6
        )
    ),
    C4s_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C4s",
            this->coeffDict_,
            3.41
        )
    ),
    C5s_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C5s",
            this->coeffDict_,
            7.27
        )
    ),
    CL_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CL",
            this->coeffDict_,
            0.164
        )
    ),
    Ceta_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceta",
            this->coeffDict_,
            75.0
        )
    ),
    phiH_
    (
        dimensionedScalar::getOrAddToDict
        (
            "phiH",
            this->coeffDict_,
            0.41
        )
    ),
    sigmaK_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            2.0
        )
    ),
    sigmaOmega_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaOmega",
            this->coeffDict_,
            2.0
        )
    ),
    sigmaPhit_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaPhit",
            this->coeffDict_,
            1.0
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    phit_
    (
        IOobject
        (
            IOobject::groupName("phit", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    ebf_
    (
        IOobject
        (
            IOobject::groupName("ebf", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    phitMin_(dimensionedScalar("phitMin", phit_.dimensions(), VSMALL)),
    ebfMin_(dimensionedScalar("ebfMin", ebf_.dimensions(), Zero)),
    L2Min_(dimensionedScalar("L2Min", sqr(dimLength), SMALL))
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
    bound(phit_, phitMin_);
    bound(ebf_, ebfMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }

    if
    (
        mag(sigmaK_.value()) < VSMALL
     || mag(sigmaOmega_.value()) < VSMALL
     || mag(sigmaPhit_.value()) < VSMALL
    )
    {
        FatalErrorInFunction
            << "Non-zero values are required for the model constants:" << nl
            << "sigmaK = " << sigmaK_ << nl
            << "sigmaOmega = " << sigmaOmega_ << nl
            << "sigmaPhit = " << sigmaPhit_ << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaLagEB<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        betas_.readIfPresent(this->coeffDict());
        Cws_.readIfPresent(this->coeffDict());
        Cp1_.readIfPresent(this->coeffDict());
        Cp2_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C5_.readIfPresent(this->coeffDict());
        C1s_.readIfPresent(this->coeffDict());
        C1Tilde_.readIfPresent(this->coeffDict());
        C4s_.readIfPresent(this->coeffDict());
        C5s_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        Ceta_.readIfPresent(this->coeffDict());
        phiH_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaOmega_.readIfPresent(this->coeffDict());
        sigmaPhit_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kOmegaLagEB<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Construct local convenience references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volSymmTensorField> tS(symm(fvc::grad(U)));
    volScalarField G(this->GName(), nut*(2.0*(dev(tS()) && tS())));
    tS.clear();

    const volScalarField L2(type() + "L2", sqr(Ls()) + L2Min_);
 
    // Compute strain, vorticity and anisotropy tensors
    tmp<volTensorField> tgradU = fvc::grad(U);

    // Mean strain rate tensor
    const volSymmTensorField S
    (
        symm(tgradU())
    );

    volScalarField magS
    (
        sqrt(2.0)*mag(S)
    );

    // Mean vorticity tensor
        const volTensorField W
    (
        0.5*(tgradU() - tgradU().T())
    );
    tgradU.clear();

    const volSymmTensorField DSDiv_(fvc::ddt(S) + fvc::div(this->phi(), S) - S * fvc::div(U));
                                
    const volTensorField SDS 
    (
        (S & DSDiv_.T())/(2.0*magSqr(S))
    );

    // Spalart-Shur curvature correction for vorticity tensor (TLLP:Eq.20)
    const volTensorField WTilde
    (
        W - 2.0*skew(SDS)
    );

    // Anisotropy tensor (TLLP:Eq.18)
    const dimensionedScalar beta2_ = scalar(2.0)*(scalar(1.0) - C5_)/(C1_ + C1s_ + scalar(1.0));
    volTensorField A
    (
        -scalar(2.0)*nut/k_*(S + scalar(2.0)*beta2_*((S & WTilde) - (WTilde & S))/
        (mag(S + WTilde)))
    );

    volScalarField taus(1.0/magS);
    forAll(taus, celli)
    {
        if(magS[celli] < 1.0e-8)
        {
            taus[celli] = 0.0;
        }
    }

    volScalarField Cp3Lim(ebf_*omega_/magS);
    volScalarField Cp3(0.0*Cp3Lim);
    forAll(Cp3, celli)
    {
        if(Cp3Lim[celli] < 10.0)
        {
            Cp3[celli] = 1.7;
        }
        else
        {
            Cp3[celli] = 0.64;
        }
    }

    // Wall-normal vectors defined through the elliptic blending factor
    volVectorField n 
    (
        fvc::grad(ebf_)/max(
            mag(fvc::grad(ebf_)), dimensionedScalar(dimless/dimLength, SMALL)
            )
    );
    
    // Update epsilon and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // Turbulent specific dissipation rate equation  (BDM:Eq.19)
    // omega
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
      ==
      - fvm::SuSp
        (   
            alpha()*rho()*((1-pow3(ebf_()))*gamma1_ + pow3(ebf_())*gamma2_)*
            (
                - G()/k_()
                + 2.0/3.0*divU 
            ),
            omega_
        )
      - fvm::Sp(alpha()*rho()*beta_*omega_, omega_)
      + fvOptions(alpha, rho, omega_)
    );
    
    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation (BDM:Eq.19)
    // k
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
      ==
        alpha()*rho()*G()
      - fvm::SuSp(2.0/3.0*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*betas_*omega_(), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


// Elliptic blending function equation (BDM:Eq.20)
    tmp<fvScalarMatrix> ebfEqn
   (
    -fvm::laplacian(ebf_)
    ==
    -fvm::SuSp(1.0/L2(), ebf_) + 1.0/L2()
   );

    ebfEqn.ref().relax();
    fvOptions.constrain(ebfEqn.ref());
    solve(ebfEqn);
    fvOptions.correct(ebf_);
    bound(ebf_, this->ebfMin_);

// Stress-strain lag equation (BDM:Eq.21)
    tmp<fvScalarMatrix> phitEqn
    (
        fvm::ddt(alpha, rho, phit_)
      + fvm::div(alphaRhoPhi, phit_)
      - fvm::laplacian(alpha*rho*DphitEff(), phit_)
      ==
      - fvm::SuSp
        (
            alpha()*rho()*
            (
                (1.0 - pow3(ebf_()))*Cws_*betas_*omega_()
                + pow3(ebf_())*(C1Tilde_ + C1s_*(G() - (2.0/3.0)*k_()*divU)
                /(betas_*k_()*omega_()))*betas_*omega_()
                + Cp1_*(G()/k_() - (2.0/3.0)*divU)
                - pow3(ebf_())*Cp2_*sqrt(2.0)*mag(S())
            )
          , phit_
        )
        + alpha()*rho()*
        (
            pow3(ebf_())*Cp3/phiH_*betas_*omega_()
            + pow3(ebf_())*(betas_*omega_()/phiH_)*sqr(taus())
            *((C4s_*(A() & S()) - C5s_*(A() & WTilde())) && S())
        )
      + fvOptions(alpha, rho, phit_)
    );

    phitEqn.ref().relax();
    fvOptions.constrain(phitEqn.ref());
    solve(phitEqn);
    fvOptions.correct(phit_);
    bound(phit_, phitMin_);
    
    
    // Bounding phit
    
    forAll(phit_, celli)
    {
        if(phit_[celli] > 2.0)
        {
            phit_[celli] = 2.0;
        }
    }
    
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
