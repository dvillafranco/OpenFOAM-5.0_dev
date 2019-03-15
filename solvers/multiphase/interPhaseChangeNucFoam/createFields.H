Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field pSat\n" << endl;
volScalarField pSat
(
    IOobject
    (
	"pSat",
        runTime.timeName(),
 	mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field n\n" << endl;
volScalarField n_
(
    IOobject
    (
	"n_",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field rRb_set\n" << endl;
volScalarField rRb_set
(
    IOobject
    (
	"rRb_set",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Creating phaseChangeTwoPhaseMixtureNuc\n" << endl;
autoPtr<phaseChangeTwoPhaseMixtureNuc> mixture =
    phaseChangeTwoPhaseMixtureNuc::New(U, phi);

volScalarField& alpha1(mixture->alpha1());
volScalarField& alpha2(mixture->alpha2());

const dimensionedScalar& rho1 = mixture->rho1();
const dimensionedScalar& rho2 = mixture->rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();
#include "phaseChangeTwoPhaseMixtureNuc.H"

volScalarField vDotcAlpha1Write
        (
                IOobject
                        (
                                "vDotcAlpha1Write",
                                mesh.time().timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                        ),
			mesh,
			dimensionedScalar("vDotcAlpha1Write",dimless/dimTime,0.0)
			
        );

volScalarField vDotvAlpha1Write
    	(
	 	IOobject
			(
				"vDotvAlpha1Write",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("vDotvAlpha1Write",dimless/dimTime,0.0)
	);

volScalarField vDotvmcAlpha1Write
	(
		IOobject
			(
				"vDotvmcAlpha1Write",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("vDotvmcAlpha1Write",dimless/dimTime,0.0)
	);

volScalarField mDotcAlphalWrite
	(
		IOobject
			(
				"mDotcAlphalWrite",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("mDotcAlphalWrite",dimless/dimTime,0.0)
	);

volScalarField mDotvAlphalWrite
	(
		IOobject
			(
				"mDotvAlphalWrite",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("mDotvAlphalWrite",dimless/dimTime,0.0)
	);

// Construct interface from alpha1 distribution
interfaceProperties interface(alpha1, U, mixture());

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture())
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());
