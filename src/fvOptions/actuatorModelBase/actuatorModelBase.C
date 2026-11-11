/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of turbinesFoam, which is based on OpenFOAM.

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

#include "actuatorModelBase.H"
#include "actuatorLineSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorModelBase, 0);
    addToRunTimeSelectionTable
    (
        option,
        actuatorModelBase,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::actuatorModelBase::calculateALData
(
    label fieldI
)
{
    // check if rho exists -> compressible simulation
    if (rhoPtr_ == nullptr && mesh_.foundObject<volScalarField>("rho"))
    {
        rhoPtr_ = &mesh_.lookupObject<volScalarField>("rho");
        // check if mu field exists -> compressible simulation
        if (muPtr_ == nullptr && mesh_.foundObject<volScalarField>("thermo:mu"))
        {
            muPtr_ = &mesh_.lookupObject<volScalarField>("thermo:mu");
        }
    }

    bool update = false;

    if (searchEnginePtr.valid() == false)
    {
        searchEnginePtr.reset(new meshSearch(mesh_));
    }

    // set up lists to hold the cell indices and processor numbers for all
    // points in all elements
    if (initialized_ == false)
    {
        allocateAL();
        label nTotal = 0;
        label nEpsilon = 0;
        forAll(actuatorLines_, i)
        {
            forAll(actuatorLines_[i]->elements(), j)
            {
                nTotal += actuatorLines_[i]->elements()[j].locationCount();
                nEpsilon += actuatorLines_[i]->elements()[j].epsilonCount();
            }
        }

        cellI_.resize(nTotal);
        procI_.resize(nTotal);
        locations_.resize(nTotal);
        velocities_.resize(nTotal);
        epsilon_.resize(nEpsilon);
        if (rhoPtr_ != nullptr)
        {
            rho_.resize(nEpsilon);
            if (muPtr_ != nullptr)
            {
                mu_.resize(nEpsilon);
            }
        }

        // We need to do an initial sweep to find the center points first for
        // proper interpolation of the epsilon points.
        List<label> centerCellI(nEpsilon, -1);
        List<label> centerProcI(nEpsilon, -1);
        List<point> centerLocations(nEpsilon, vector::zero);
        
        setupPositions(false);
        
        label index = 0;
        label myProcNo = Pstream::myProcNo();
        forAll(actuatorLines_, i)
        {
            forAll(actuatorLines_[i]->elements(), j)
            {
                actuatorLines_[i]->elements()[j].collectLocationData
                (
                    centerCellI,
                    centerProcI,
                    centerLocations,
                    index,
                    false
                );
            }
        }
        if (index > nEpsilon)
        {
            FatalErrorInFunction
                << "index = " << index
                << " nEpsilon = " << nEpsilon
                << abort(FatalError);
        }

        // We actually only need to collect data on the cell index to see
        // failed cells, processes only need to know if they own it or not
        reduce(centerProcI, maxOp<List<label>>());

        // check the cells we couldn't find in first sweep
        forAll(centerProcI, i)
        {
            if (centerProcI[i] < 0)
            {
                if (meshBoundBox_.containsInside(centerLocations[i]))
                {
                    centerCellI[i] = searchEnginePtr->findCell
                    (
                        centerLocations[i]
                    );
                    if (centerCellI[i] >= 0)
                    {
                        centerProcI[i] = myProcNo;
                    }
                }
            }
        }
        index = 0;
        forAll(actuatorLines_, i)
        {
            forAll(actuatorLines_[i]->elements(), j)
            {
                actuatorLines_[i]->elements()[j].distributeCenterCellI
                (
                    centerCellI,
                    centerProcI,
                    index
                );
            }
        }
        initializeAL();
        initialized_ = true;
    }

    // virtual function that should perform the first step of
    // finding mesh locations of the relevant points
    if (lastTime_ != mesh_.time().value())
    {

        setupPositions(true);
        update = true;
        lastTime_ = mesh_.time().value();
        firstField_ = fieldI;
    }

    // for outer iterations, we only want to update the inflow velocities
    // if this is the first field to be called
    if (firstField_ == fieldI)
    {
        update = true;
    }

    if (update)
    {
        // in case mesh is updated
        if (mesh_.changing())
        {
            searchEnginePtr.reset(new meshSearch(mesh_));
        }
        // Gather data from the elements
        label index = 0;
        label myProcNo = Pstream::myProcNo();
        forAll(actuatorLines_, i)
        {
            forAll(actuatorLines_[i]->elements(), j)
            {
                actuatorLines_[i]->elements()[j].collectLocationData
                (
                    cellI_,
                    procI_,
                    locations_,
                    index
                );
            }
        }
        // We need to ensure that only one processor own one cell
        // (multiple processes can find the same cell)
        reduce(procI_, maxOp<List<label>>());

        // check the cells we couldn't find in first sweep
        forAll(procI_, i)
        {
            if (procI_[i] < 0)
            {
                if (meshBoundBox_.containsInside(locations_[i]))
                {
                    cellI_[i] = searchEnginePtr->findCell(locations_[i]);
                    if (cellI_[i] >= 0)
                    {
                        procI_[i] = myProcNo;
                    }
                }
            }
        }

        // This part is only a safety check to ensure that the point is in mesh
        // can be removed if this is not needed
        reduce(procI_, maxOp<List<label>>());
        forAll(procI_, i)
        {
            if (procI_[i] < 0)
            {
                // Raise fatal error since inflow velocity cannot be detected
                FatalErrorIn("void actuatorModelBase::calculateALData()")
                    << "Inflow velocity point for position: "
                    << locations_[i] << " not found in mesh"
                    << abort(FatalError);
            }
        }


        // Velocity interpolation
        velocities_ = vector::zero;
        interpolationCellPoint<vector> UInterp
        (
            mesh_.lookupObject<volVectorField>("U")
        );

        forAll(cellI_, i)
        {
            // Only interpolate if owner
            if (procI_[i] == myProcNo)
            {
                velocities_[i] = UInterp.interpolate
                (
                    locations_[i],
                    cellI_[i]
                );
            }
        }
        reduce(velocities_, sumOp<List<vector>>());
        index = 0;
        forAll(actuatorLines_, i)
        {
            forAll(actuatorLines_[i]->elements(), j)
            {
                actuatorLines_[i]->elements()[j].distributeVelocityData
                (
                    cellI_,
                    procI_,
                    velocities_,
                    index
                );
            }
        }

        calculateForces();

        distributeEpsilon();

        // Compressible simulation, multiply force by density
        if (rhoPtr_ != nullptr)
        {
            index = 0;
            forAll(cellI_, i)
            {
                // Only check if owner, otherwise set to large value so minOp
                // will ignore it
                if (procI_[i] == myProcNo)
                {
                    // only small variations in rho, pick local cell value
                    rho_[i] = (*rhoPtr_)[cellI_[i]];
                }
                else
                {
                    rho_[i] = VGREAT;
                }
            }
            reduce(rho_, minOp<List<scalar>>());
            index = 0;
            forAll(actuatorLines_, i)
            {
                forAll(actuatorLines_[i]->elements(), j)
                {
                    actuatorLines_[i]->elements()[j].distributeRhoData
                    (
                        rho_,
                        index
                    );
                }
            }

            if (muPtr_ != nullptr)
            {
                index = 0;
                forAll(cellI_, i)
                {
                    // Only check if owner, otherwise set to large value
                    // so minOp will ignore it
                    if (procI_[i] == myProcNo)
                    {
                        // only small variations in rho, pick local cell value
                        mu_[i] = (*muPtr_)[cellI_[i]];
                    }
                    else
                    {
                        mu_[i] = VGREAT;
                    }
                }
                reduce(mu_, minOp<List<scalar>>());
                index = 0;
                forAll(actuatorLines_, i)
                {
                    forAll(actuatorLines_[i]->elements(), j)
                    {
                        actuatorLines_[i]->elements()[j].distributeMuData
                        (
                            mu_,
                            index
                        );
                    }
                }
            }
        }
    }
}

void Foam::fv::actuatorModelBase::distributeEpsilon()
{
    label index = 0;
    forAll(actuatorLines_, i)
    {
        forAll(actuatorLines_[i]->elements(), j)
        {
            actuatorLines_[i]->elements()[j].collectEpsilonData
            (
                epsilon_,
                index
            );
        }
    }
    reduce(epsilon_, minOp<List<scalar>>());

    index = 0;
    forAll(actuatorLines_, i)
    {
        forAll(actuatorLines_[i]->elements(), j)
        {
            actuatorLines_[i]->elements()[j].distributeEpsilonData
            (
                epsilon_,
                index
            );
        }
    }
}

void Foam::fv::actuatorModelBase::createForceField()
{
    // Already created, do nothing
    if (forceFieldPtr_.valid())
    {
        return;
    }

    forceFieldPtr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                name() + forceFieldName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                writeForceField_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("Force", dimForce/dimVolume, Zero)
        )
    );
    //forceFieldPtr_().write();
}

// Dummy functions unless one overloads it
void Foam::fv::actuatorModelBase::allocateAL()
{
}

void Foam::fv::actuatorModelBase::initializeAL()
{
}

void Foam::fv::actuatorModelBase::setupPositions(bool includeRing)
{
}

void Foam::fv::actuatorModelBase::calculateForces()
{
}

const List<Foam::fv::actuatorLineSource*>&
Foam::fv::actuatorModelBase::actuatorLines() const
{
    return actuatorLines_;
}

void Foam::fv::actuatorModelBase::writeFieldsAtEnd()
{
    if
    (
        writeFieldsAtEnd_
     && not finalFieldsWritten_
     && (
             mesh_.time().value()
           + 0.5*mesh_.time().deltaT().value()
           >= mesh_.time().endTime().value()
        )
    )
    {
        //Info<< "Writing fields at final time" << endl;
        Info<< "turbinesFoam: forcing field write at time "
            << mesh_.time().value()
            << endl;

        Time& runTime = const_cast<Time&>(mesh_.time());
        runTime.writeNow();

        finalFieldsWritten_ = true;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::actuatorModelBase::actuatorModelBase(
    const word &name,
    const word &modelType,
    const dictionary &dict,
    const fvMesh &mesh)
    : cellSetOption(name, modelType, dict, mesh),
      meshBoundBox_(mesh_.points(), false),
      velocities_(0),
      locations_(0),
      cellI_(0),
      procI_(0),
      epsilon_(0),
      rho_(0),
      mu_(0),
      actuatorLines_(0),
      rhoPtr_(nullptr),
      muPtr_(nullptr),
      initialized_(false),
      writeFieldsAtEnd_(false),
      finalFieldsWritten_(false)
{
    meshBoundBox_.inflate(1e-6);
    read(dict);
    if (writeForceField_)
    {
        forceFieldPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    name + forceFieldName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector("Force", dimForce/dimVolume, Zero)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::actuatorModelBase::~actuatorModelBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::actuatorModelBase::addForce
(
    volVectorField &forceField,
    scalar scale,
    bool compressible
)
{
}

void Foam::fv::actuatorModelBase::addTurbulence(
    fvMatrix<scalar> &eqn,
    const word fieldName)
{
}

// Used to determine if force should be saved to a local force field or directy
// to the owners force field. (Do we really need this functionality though?)
void Foam::fv::actuatorModelBase::addForceFromChild
(
    volVectorField& forceField,
    scalar scale,
    bool compressible
)
{
    // Select if we should use local force field or input argument
    const bool useLocal = forceFieldPtr_.valid();
    volVectorField& target = useLocal ? forceFieldPtr_() : forceField;

    // zero our own force field if we use it
    if (useLocal)
    {
        target.primitiveFieldRef() = vector::zero;
        if (target.dimensions() != forceField.dimensions()/dimVolume)
        {
            target.dimensions().reset(forceField.dimensions());
        }
    }

    addForce(target, scale, compressible);

    // If we use local field, add it to the input field
    if (useLocal)
    {
        forceField += target;
        // we multiply with rho after, as this is handled by addSup
        // for the main field
        if (compressible)
        {
            // First time, we need to save rho
            if (rhoPtr_ == nullptr && mesh_.foundObject<volScalarField>("rho"))
            {
                rhoPtr_ = &mesh_.lookupObject<volScalarField>("rho");
            }
            // rho should exist for compressible simulations, but check anyway
            if (rhoPtr_ != nullptr)
            {
                target *= *rhoPtr_;
            }
        }
    }
}

void Foam::fv::actuatorModelBase::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    if (initialized_ == false)
    {
        createForceField();
    }
    volVectorField& forceField = forceFieldPtr_();
    forceField.primitiveFieldRef() = vector::zero;
    
    // Should not be needed?
    if (forceField.dimensions() != eqn.dimensions()/dimVolume)
    {
        forceField.dimensions().reset(eqn.dimensions()/dimVolume);
    }

    calculateALData(fieldI);

    addForce(forceField, 1.0, false);

    forceField.correctBoundaryConditions();

    eqn += forceField;

    writeFieldsAtEnd();
}


void Foam::fv::actuatorModelBase::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    if (initialized_ == false)
    {
        createForceField();
    }
    volVectorField& forceField = forceFieldPtr_();
    forceField.primitiveFieldRef() = vector::zero;

    if (forceField.dimensions() != eqn.dimensions()/dimVolume/dimDensity)
    {
        forceField.dimensions().reset(eqn.dimensions()/dimVolume/dimDensity);
    }

    calculateALData(fieldI);

    addForce(forceField, 1.0, true);

    // multiply with local density
    forceField *= rho;

    forceField.correctBoundaryConditions();
    
    eqn += forceField;

    writeFieldsAtEnd();
}


void Foam::fv::actuatorModelBase::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    calculateALData(fieldI);

    word fieldName = fieldNames_[fieldI];
    Info<< endl << "Adding " << fieldName << " from " << name_ << endl << endl;
    addTurbulence(eqn, fieldName);

    writeFieldsAtEnd();
}

bool Foam::fv::actuatorModelBase::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        writeForceField_ = coeffs_.lookupOrDefault
        (
            "writeForceField",
            true
        );
        writeFieldsAtEnd_ = coeffs_.lookupOrDefault
        (
            "writeFieldsAtEnd",
            false
        );
        printPerf_ = coeffs_.lookupOrDefault
        (
            "printPerf",
            true
        );

        forceFieldName_ = coeffs_.lookupOrDefault<word>
        (
            "forceFieldName",
            ":force"
        );
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
