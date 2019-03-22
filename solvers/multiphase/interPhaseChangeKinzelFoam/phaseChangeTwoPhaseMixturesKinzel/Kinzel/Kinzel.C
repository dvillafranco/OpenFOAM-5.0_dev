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

    return pow
    (
	((2/3)*(max(p-pSat(), p0_))/rho1()+((2/3)*100000/rho1()*1/(1-1.4)*(alphaNuc/(1+alphaNuc-limitedAlpha1)))),
        1.0/2.0
    );
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
    return pow
    (
        ((2/3)*(min(p-pSat(), p0_))/rho1()+((2/3)*100000/rho1()*1/(1-1.4)*(alphaNuc/(1+alphaNuc-limitedAlpha1)))),
        1.0/2.0
    );
}

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::alphaLiqCoeff() const

{
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
     
    return pow
    (
       1-limitedAlpha1,
       1.0/6.0
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
   
    return Pair<tmp<volScalarField>>
    (
        Cc_*pCoeffCond*alphaLiqCoeff()*C_infty(),
       	Cv_*pCoeffVap*alphaLiqCoeff()*C_infty()
    );
    
}

//else
//{
//Foam::Pair<Foam::tmp<Foam:volScalarField>>
//Foam::phaseChangeTwoPhaseMixtures::SchnerrSauer::mDotAlphal() const
//{
//   const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
//   volScalarField pCoeff(this->pCoeff(p));
//
//   volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
//
//  return Pair<tmp<volScalarField>>
//   (
//     Cc_limitedAlpha1*pCoeff*p0_*max(p-pSat(),p0_),
//    
//     Cv_*(1.0+alphaNuc()-limitedAlpha1)*pCoeff*min(p-pSat(), p0_)
//   );
//}
//}
//}    
//	return Pair<tmp<volScalarField>>
//	(
//		p0_,
//
//		Cv_*(1.0 + alphaNuc() - limitedAlpha1)*pCoeff*min(p-pSat(),p0_)
//	);
    



Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixturesNuc::Kinzel::mDotP() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeffCond(this->pCoeffCond(p));
    const volScalarField& n_ = alpha1_.db().lookupObject<volScalarField>("n_");
    volScalarField alphaNuc(this->alphaNuc(n_));

    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    volScalarField apCoeff(limitedAlpha1*pCoeffCond);

    return Pair<tmp<volScalarField>>
    (
        Cc_*(1.0 - limitedAlpha1)*pos0(p - pSat())*apCoeff,

        (-Cv_)*(1.0 + alphaNuc - limitedAlpha1)*neg(p - pSat())*apCoeff
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
