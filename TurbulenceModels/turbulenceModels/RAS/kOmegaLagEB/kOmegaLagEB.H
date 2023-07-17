/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

Class
    Foam::RASModels::kOmegaLagEB

Group
    grpRASTurbulence

Description
    The k-epsilon Lag Elliptic-Blending turbulence closure model for incompressible and
    compressible flows.

    The model is a three-transport-equation linear-eddy-viscosity turbulence
    closure model alongside an elliptic blending equation.

    \heading Input fields
    \plaintable
        k         | Turbulent kinetic energy [m2/s2]
        omega        
        phit      | Lag between stress-strain [-]
        ebf         | Elliptic blending factor [-]
    \endplaintable

    Reference:
    \verbatim
        Standard model (Tag:TLLP):
            Laurence
            Tunstall, R., Lardeau, S., Laurence, D., Prosser, R. (2016).
            An Elliptic Blending Lag Model for Flows in Thermal-Hydraulics Systems
            ETMM11 Conference
            https://www.researchgate.net/publication/308413285_AN_ELLIPTIC_BLENDING_LAG_MODEL_FOR_FLOWS_IN_THERMAL-HYDRAULICS_SYSTEMS
    \endverbatim

    The default model coefficients are:
    \verbatim
        kOmegaLagEBCoeffs
        {
            Cmu         0.22;       // Turbulent viscosity constant
            betas       0.09;
            C1          1.7;
            C3          0.8;
            C4          0.625;
            C5          0.2;
            C1s         0.9;
            C2s         0.65;       
            CK          2.3;        // Model constant for E
            Ct          4.0;        // Model constant for T
            Ceps1       1.44;       // Model constant for epsilon
            Ceps2       1.9;        // Model constant for epsilon
            CL          0.164;      // Model constant for L
            Ceta        75.0;       // Model constant for L
            sigmaK      1.0;        // Turbulent Prandtl number for k
            sigmaEps    1.2;        // Turbulent Prandtl number for epsilon
            sigmaPhit   1.0;        // Turbulent Prandtl number for phit
        }
    \endverbatim

Note
    The name of the variable representing the lag between stress and strain is 'phi'.
    However, the name 'phi' preexisted in OpenFOAM; therefore, this name was
    replaced by 'phit' herein.
    The same applies for the elliptic blending factor (called 'alpha'). 
    Therefore, it has been replaced by 'ebf' herein.

SourceFiles
    kOmegaLagEB.C

SeeAlso
    kOmega.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaLagEB_H
#define kOmegaLagEB_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class kOmegaLagEB Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaLagEB
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        kOmegaLagEB(const kOmegaLagEB&) = delete;

        //- No copy assignment
        void operator=(const kOmegaLagEB&) = delete;


protected:

    // Protected Data
        // Model coefficients

            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;
            dimensionedScalar beta_;
            dimensionedScalar betas_;
            dimensionedScalar Cws_;
            dimensionedScalar Cp1_;
            dimensionedScalar Cp2_;
            dimensionedScalar C1_;
            dimensionedScalar C5_;
            dimensionedScalar C1s_;
            dimensionedScalar C1Tilde_;
            dimensionedScalar C4s_;
            dimensionedScalar C5s_;
            dimensionedScalar CL_;
            dimensionedScalar Ceta_;
            dimensionedScalar phiH_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaOmega_;
            dimensionedScalar sigmaPhit_;


        // Fields

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Turbulent kinetic energy dissipation rate [m2/s3]
            volScalarField omega_;

            //- Lag between stress-strain [-]
            volScalarField phit_;

            //- Elliptic blending factor [-]
            volScalarField ebf_;


        // Bounding values

            dimensionedScalar phitMin_;
            dimensionedScalar ebfMin_;
            dimensionedScalar L2Min_;

    // Protected Member Functions

        //- Update nut with the latest available k, phit, and T
        virtual void correctNut();

        //- Return the turbulent length scale, L
        tmp<volScalarField> Ls() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaLagEB");


    // Constructors

        //- Construct from components
        kOmegaLagEB
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaLagEB() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k (TLLP:Eq.5)
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    this->nut_/sigmaK_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for epsilon (TLLP:Eq.6)
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    this->nut_/sigmaOmega_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for phit (TLLP:Eq.7)
        tmp<volScalarField> DphitEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DphitEff",
                    this->nut_/sigmaPhit_ + this->nu()
                )
            );
        }

        //- Return the turbulent kinetic energy field
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulent kinetic energy dissipation rate field
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the normalised wall-normal fluctuating velocity scale field
        virtual tmp<volScalarField> phit() const
        {
            return phit_;
        }

        //- Return the elliptic relaxation factor field
        virtual tmp<volScalarField> ebf() const
        {
            return ebf_;
        }

        //- Solve the transport equations and correct the turbulent viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaLagEB.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
