/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "Kinzel.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixturesNuc
{
    defineTypeNameAndDebug(Kinzel, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeTwoPhaseMixtureNuc,
        Kinzel,
        components
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::Kinzel
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixtureNuc(typeName, U, phi),

   // n_("n", dimless/dimVolume, phaseChangeTwoPhaseMixtureCoeffs_),
    dNuc_("dNuc", dimLength, phaseChangeTwoPhaseMixtureCoeffs_),
    Cc_("Cc", dimless, phaseChangeTwoPhaseMixtureCoeffs_),
    Cv_("Cv", dimless, phaseChangeTwoPhaseMixtureCoeffs_),
    Pref_("Pref",pSat().dimensions(),phaseChangeTwoPhaseMixtureCoeffs_),
    p0_("0", pSat().dimensions(), 0.0)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::C_infty() const
{
    const volScalarField& n_ = alpha1_.db().lookupObject<volScalarField>("n_");
    volScalarField alphaNuc(this->alphaNuc(n_));
    //const volScalarField& rRb_set = alpha1_.db().lookupObject<volScalarField>("rRb_set");
    return pow
    (
        (constant::mathematical::pi*n_)/36,
        1.0/3.0
    );
}


Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::alphaNuc
(
    const volScalarField& n_
) const
{
    volScalarField Vnuc = n_*constant::mathematical::pi*pow3(dNuc_)/6;
    return Vnuc/(1 + Vnuc);
}

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::alphaGammaTermVap() const
{
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    const volScalarField& n_ = alpha1_.db().lookupObject<volScalarField>("n_");
    volScalarField alphaNuc(this->alphaNuc(n_));
    return pow
    (
     alphaNuc/(1+alphaNuc -limitedAlpha1),
     7.0/5.0
    );
}


Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::pCoeffCond
(
    const volScalarField& p
) const
{
    const volScalarField& n_ = alpha1_.db().lookupObject<volScalarField>("n_");
    volScalarField alphaNuc(this->alphaNuc(n_));
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    volScalarField rho
    (
        limitedAlpha1*rho1() + (scalar(1) - limitedAlpha1)*rho2()
    );

    return sqrt(mag((2.0/3.0)*(p-pSat())/rho1()+((2.0/3.0)*Pref_/rho1()*1.0/(1.0-1.4)*(alphaGammaTermVap()))));
}

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::pCoeffVap
(
   const volScalarField& p
) const
{
    const volScalarField& n_ = alpha1_.db().lookupObject<volScalarField>("n_");
    volScalarField alphaNuc(this->alphaNuc(n_));
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    return sqrt(mag(((2.0/3.0)*((p-pSat())))/rho1()+(((2.0/3.0)*Pref_/rho1())*(1.0/(1.0-1.4))*(alphaGammaTermVap()-(alphaNuc/(1+alphaNuc-limitedAlpha1))))));
}

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::alphaLiqCoeff() const

{
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    const volScalarField& n_ = alpha1_.db().lookupObject<volScalarField>("n_");
    volScalarField alphaNuc(this->alphaNuc(n_)); 
    return pow
    (
       1.0+alphaNuc-limitedAlpha1,
       1.0/6.0
    );
}

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::alphaLiqCoeffVap() const
{
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    const volScalarField& n_ = alpha1_.db().lookupObject<volScalarField>("n_");
    volScalarField alphaNuc(this->alphaNuc(n_));
    return pow
    (
      1.0+alphaNuc-limitedAlpha1,
      2.0/3.0
    );
}




Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::mDotAlphal() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeffCond(this->pCoeffCond(p));
    volScalarField pCoeffVap(this->pCoeffVap(p));
    const volScalarField& n_ = alpha1_.db().lookupObject<volScalarField>("n_");
    volScalarField alphaNuc(this->alphaNuc(n_));

    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    volScalarField rho
    (
        limitedAlpha1*rho1() + (scalar(1) - limitedAlpha1)*rho2()
    );
   
    return Pair<tmp<volScalarField>>
    (
        Cc_*pCoeffCond*pos0((((2.0/3.0)*((p-pSat())))/rho1()+(((2.0/3.0)*Pref_/rho1())*(1.0/(1.0-1.4))*(alphaGammaTermVap()))))*alphaLiqCoeff()*C_infty()*rho1()*rho2()/rho,
       	-Cv_*pCoeffVap*neg((((2.0/3.0)*((p-pSat())))/rho1()+(((2.0/3.0)*Pref_/rho1())*(1.0/(1.0-1.4))*(alphaGammaTermVap()-(alphaNuc/(1+alphaNuc-limitedAlpha1))))))*alphaLiqCoeffVap()*C_infty()*rho1()*rho2()/rho
    );
    
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::pCoeffChange() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeffCond(this->pCoeffCond(p));
    volScalarField pCoeffVap(this->pCoeffVap(p));
    return Pair<tmp<volScalarField>>
    (
     1.0*pCoeffCond,
     1.0*pCoeffVap
    );
}
      

    
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::mDotP() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeffCond(this->pCoeffCond(p));
    volScalarField pCoeffVap(this->pCoeffVap(p));
    const volScalarField& n_ = alpha1_.db().lookupObject<volScalarField>("n_");
    volScalarField alphaNuc(this->alphaNuc(n_));

    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    volScalarField apCoeffCond(limitedAlpha1*pCoeffCond);
    volScalarField apCoeffVap(limitedAlpha1*pCoeffVap);
    volScalarField rho
    (
	limitedAlpha1*rho1() + (scalar(1) - limitedAlpha1)*rho2()
    );
    return Pair<tmp<volScalarField>>
    (
        Cc_*(1.0 - limitedAlpha1)*pos0(p - pSat())*apCoeffCond*rho*C_infty()/(mag(p-pSat())+0.01*pSat()),

        (-Cv_)*(1.0 + alphaNuc - limitedAlpha1)*neg(p - pSat())*apCoeffVap*rho*C_infty()/(mag(p-pSat())+0.01*pSat())
    );
}


void Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::correct()
{}


bool Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::read()
{
    if (phaseChangeTwoPhaseMixtureNuc::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = optionalSubDict(type() + "Coeffs");

       // phaseChangeTwoPhaseMixtureCoeffs_.lookup("n") >> n_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("dNuc") >> dNuc_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cc") >> Cc_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv") >> Cv_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
