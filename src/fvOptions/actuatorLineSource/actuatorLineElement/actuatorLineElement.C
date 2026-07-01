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

#include "actuatorLineElement.H"
#include "addToRunTimeSelectionTable.H"
#include "geometricOneField.H"
#include "fvMatrices.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorLineElement, 0);
    defineRunTimeSelectionTable(actuatorLineElement, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::actuatorLineElement::read()
{
    // Parse dictionary
    dict_.lookup("position") >> position_;
    dict_.lookup("chordLength") >> chordLength_;
    dict_.lookup("chordDirection") >> chordDirection_;
    dict_.lookup("chordRefDirection") >> chordRefDirection_;
    dict_.lookup("chordMount") >> chordMount_;
    dict_.lookup("spanLength") >> spanLength_;
    dict_.lookup("spanDirection") >> spanDirection_;
    dict_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
    freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
    dict_.lookup("rootDistance") >> rootDistance_;
    dict_.lookup("cone") >> cone_;
    dict_.lookup("velocitySampleRadius") >> velocitySampleRadius_;
    dict_.lookup("nVelocitySamples") >> nVelocitySamples_;

    // Create dynamic stall model if found
    if (dict_.found("dynamicStall"))
    {
        dictionary dsDict = dict_.subDict("dynamicStall");
        word dsName;
        dsDict.lookup("dynamicStallModel") >> dsName;
        dynamicStall_ = dynamicStallModel::New
        (
            dsDict,
            dsName,
            mesh_.time(),
            profileData_
        );
        dsDict.lookup("active") >> dynamicStallActive_;
    }

    // Read flow curvature correction subdictionary
    if (dict_.found("flowCurvature"))
    {
        dictionary fcDict = dict_.subDict("flowCurvature");
        flowCurvatureActive_ = fcDict.lookupOrDefault("active", false);
        word defaultName = "none";
        flowCurvatureModelName_ = fcDict.lookupOrDefault
        (
            "flowCurvatureModel",
            defaultName
        );
    }

    // Read nu from object registry
    const dictionary& transportProperties = mesh_.lookupObject<IOdictionary>
    (
        "transportProperties"
    );
    dimensionedScalar nu;
    transportProperties.lookup("nu") >> nu;
    nu_ = nu.value();

    // Read writePerf switch
    dict_.lookup("writePerf") >> writePerf_;
    dict_.lookup("writePerfEnd") >> writePerfEnd_;

    if (debug)
    {
        Info<< "actuatorLineElement properties:" << endl;
        Info<< "Position: " << position_ << endl;
        Info<< "chordLength: " << chordLength_ << endl;
        Info<< "chordDirection: " << chordDirection_ << endl;
        Info<< "spanLength: " << spanLength_ << endl;
        Info<< "spanDirection: " << spanDirection_ << endl;
        Info<< "cone: " << cone_ << endl;
        Info<< "writePerf: " << writePerf_ << endl;
        Info<< "writePerfEnd: " << writePerfEnd_ << endl;
    }
}


void Foam::fv::actuatorLineElement::rotateVector
(
    vector& vectorToRotate,
    vector rotationPoint,
    vector axis,
    scalar radians
)
{
    // Declare and define the rotation matrix (from SOWFA)
    tensor RM;
    scalar angle = radians;
    RM.xx() = Foam::sqr(axis.x())
            + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle);
    RM.xy() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle);
    RM.xz() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle);
    RM.yy() = Foam::sqr(axis.y())
            + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z())
            + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    // Rotation matrices make a rotation about the origin, so need to subtract
    // rotation point off the point to be rotated.
    vectorToRotate -= rotationPoint;

    // Perform the rotation.
    vectorToRotate = RM & vectorToRotate;

    // Return the rotated point to its new location relative to the rotation
    // point
    vectorToRotate += rotationPoint;
}


Foam::label Foam::fv::actuatorLineElement::findCell
(
    const point& location
)
{
    label localCell = -1;

    // Fastest path: Same point as last time
    bool reuse =
        previousLocationValid_
        && magSqr(location - previousLocation_) < SMALL;

    reduce(reuse, andOp<bool>());
    if (reuse)
    {
        return
        (
            Pstream::myProcNo() == previousProcI_
          ? previousCellI_
          : -1
        );
    }

    // Special actuator disc path, use cached index
    if (centerCellI_.size() > 0)
    {
        if (activeRingIndex_ >= 0) // ring sampling active
        {
            if (ringProcI_[azimuthIndex_][activeRingIndex_] != -2) // first use
            {
                previousProcI_ = ringProcI_[azimuthIndex_][activeRingIndex_];
                previousCellI_ = ringCellI_[azimuthIndex_][activeRingIndex_];
            }
            else
            {
                ringProcI_[azimuthIndex_][activeRingIndex_] = -1;
                ringCellI_[azimuthIndex_][activeRingIndex_] = -1;
            }
        }
        else
        {
            if (centerCellI_[azimuthIndex_] != -2)
            {
                previousProcI_ = centerProcI_[azimuthIndex_];
                previousCellI_ = centerCellI_[azimuthIndex_];
            }
            else
            {
                centerProcI_[azimuthIndex_] = -1;
                centerCellI_[azimuthIndex_] = -1;
            }
        }
    }
    // Fast path:
    // Previous owner checks cached cell + neighbors
    bool cacheMiss = false;
    if
    (
        previousProcI_ == Pstream::myProcNo()
            && previousCellI_ >= 0
            && previousCellI_ < mesh_.nCells()
    )
    {
        // Check cached cell
        if (mesh_.pointInCell(location, previousCellI_))
        {
            localCell = previousCellI_;
        }
        else
        {
            // If not previous cell, check neighboring cells
            const labelList& nbrs =
                mesh_.cellCells()[previousCellI_];

            forAll(nbrs, nbrI)
            {
                label testCell = nbrs[nbrI];

                if
                (
                    testCell >= 0
                    && testCell < mesh_.nCells()
                    && mesh_.pointInCell(location, testCell)
                )
                {
                    localCell = testCell;
                    break;
                }
            }
            cacheMiss = true;
        }
    }

    // Fallback: expensive global/local search
    if (localCell == -1)
    {
        if (meshBoundBox_.containsInside(location))
        {
            localCell = mesh_.findCell(location);
            cacheMiss = true;
        }
    }

    // Determine owning processor
    label ownerProc = (localCell >= 0)
      ? Pstream::myProcNo() : -1;

    reduce(ownerProc, maxOp<label>());

    // Update cache
    if (localCell >= 0)
    {
        previousCellI_ = localCell;
        previousProcI_ = Pstream::myProcNo();
        previousLocation_ = location;
        previousLocationValid_ = true;
        if (centerCellI_.size() > 0)
        {
            if (cacheMiss)
            {
                if (activeRingIndex_ >= 0) // ring sampling active
                {
                   ringProcI_[azimuthIndex_][activeRingIndex_] = previousProcI_;
                   ringCellI_[azimuthIndex_][activeRingIndex_] = previousCellI_;
                }
                else
                {
                    centerProcI_[azimuthIndex_] = previousProcI_;
                    centerCellI_[azimuthIndex_] = previousCellI_;
                }
            }
        }
    }
    else
    {
        previousCellI_ = -1;
        previousProcI_ = -1;
        previousLocationValid_ = false;
    }

    return localCell;
}


void Foam::fv::actuatorLineElement::lookupCoefficients()
{
    liftCoefficient_ = profileData_.liftCoefficient(angleOfAttack_);
    dragCoefficient_ = profileData_.dragCoefficient(angleOfAttack_);
    momentCoefficient_ = profileData_.momentCoefficient(angleOfAttack_);
}


Foam::scalar Foam::fv::actuatorLineElement::calcProjectionEpsilon()
{
    // Lookup Gaussian coeffs from profileData dict if present
    dictionary GaussianCoeffs = profileData_.dict().subOrEmptyDict
    (
        "GaussianCoeffs"
    );
    scalar chordFactor = GaussianCoeffs.lookupOrDefault("chordFactor", 0.25);
    scalar dragFactor = GaussianCoeffs.lookupOrDefault("dragFactor", 1.0);
    scalar meshFactor = GaussianCoeffs.lookupOrDefault("meshFactor", 2.0);

    // Provide ideal epsilon target for lift based on chord length
    scalar epsilonLift = chordFactor*chordLength_;

    // Epsilon based on drag/momentum thickness
    scalar epsilonDrag = dragFactor*dragCoefficient_*chordLength_/2.0;

    // Threshold is based on lift or drag, whichever is larger
    scalar epsilonThreshold = Foam::max(epsilonLift, epsilonDrag);

    scalar epsilon = VGREAT;
    scalar epsilonMesh = VGREAT;
    const scalarField& V = mesh_.V();
    label posCellI = findCell(position_);

    if (posCellI >= 0)
    {
        // Projection width based on local cell size (from Troldborg (2008))
        epsilonMesh = 2.0*Foam::cbrt(V[posCellI]);
        epsilonMesh *= meshFactor; // Cell could have non-unity aspect ratio

        if (epsilonMesh > epsilonThreshold)
        {
            epsilon = epsilonMesh;
        }
        else
        {
            epsilon = epsilonThreshold;
        }
    }

    // Reduce epsilon over all processors
    reduce(epsilon, minOp<scalar>());

    // If epsilon is not reduced, position is not in the mesh
    if (not (epsilon < VGREAT))
    {
        // Raise fatal error since mesh size cannot be detected
        FatalErrorIn("void actuatorLineElement::applyForceField()")
            << "Position of " << name_ << " Position: " << position_
            << " not found in mesh" << abort(FatalError);
    }

    if (debug)
    {
        reduce(epsilonMesh, minOp<scalar>());
        word epsilonMethod;
        if (epsilon == epsilonLift)
        {
            epsilonMethod = "lift-based";
        }
        else if (epsilon == epsilonDrag)
        {
            epsilonMethod = "drag-based";
        }
        else if (epsilon == epsilonMesh)
        {
            epsilonMethod = "mesh-based";
        }
        Info<< "    epsilon (" << epsilonMethod << "): " << epsilon << endl;
    }

    return epsilon;
}


void Foam::fv::actuatorLineElement::correctFlowCurvature
(
    scalar& angleOfAttackRad
)
{
    if (debug)
    {
        Info<< "    Correcting for flow curvature with "
            << flowCurvatureModelName_ << " model" << endl;
    }

    if (flowCurvatureModelName_ == "Goude")
    {
        angleOfAttackRad +=
            omega_*chordLength_/(2*mag(relativeVelocity_))*cos(cone_);
    }
    else if (flowCurvatureModelName_ == "MandalBurton")
    {
        // Calculate relative velocity at leading and trailing edge
        vector relativeVelocityLE = inflowVelocity_ - velocityLE_;
        vector relativeVelocityTE = inflowVelocity_ - velocityTE_;

        // Calculate angle of attack at leading and trailing edge
        scalar alphaLE = asin((planformNormal_ & relativeVelocityLE)
                       / (mag(planformNormal_)*mag(relativeVelocityLE)));
        scalar alphaTE = asin((planformNormal_ & relativeVelocityTE)
                       / (mag(planformNormal_)*mag(relativeVelocityTE)));

        scalar beta = alphaTE - alphaLE;

        angleOfAttackRad += atan2((1.0 - cos(beta/2.0)), sin(beta/2.0));
    }
    else if (flowCurvatureModelName_ == "constantOffset")
    {
        dictionary fcDict = dict_.subDict("flowCurvature");
        dictionary coeffs = fcDict.subDict(flowCurvatureModelName_ + "Coeffs");
        scalar offsetDeg = 0.0;
        coeffs.lookup("offsetDeg") >> offsetDeg;
        angleOfAttackRad += degToRad(offsetDeg);
    }
}


void Foam::fv::actuatorLineElement::multiplyForceRho
(
    const volScalarField& rho
)
{
    // Lookup local density
    label cellI = findCell(position_);
    scalar localRho = VGREAT;
    if (cellI >= 0)
    {
        localRho = rho[cellI];
    }

    reduce(localRho, minOp<scalar>());
    forceVector_ *= localRho;
}


void Foam::fv::actuatorLineElement::applyForceField
(
    volVectorField& forceField,
    scalar scale
)
{
    // Calculate projection width
    scalar epsilon = calcProjectionEpsilon();
    scalar projectionRadius = (epsilon*Foam::sqrt(Foam::log(1.0/0.001)));

    // Apply force to the cells within the element's sphere of influence
    scalar sphereRadius = chordLength_ + projectionRadius;
    scalar sphereRadiusSqr = sphereRadius*sphereRadius;
    scalar invepsilonSqr = 1.0/(epsilon*epsilon);
    scalar internalFactor = scale/(Foam::pow(epsilon, 3)
                          * Foam::pow(Foam::constant::mathematical::pi, 1.5));
    // forceField is opposite forceVector
    const vector scaledForce = -forceVector_*internalFactor;

    const vectorField& C = mesh_.C();

    vectorField& force =
        forceField.primitiveFieldRef();

    const scalar px = position_.x();
    const scalar py = position_.y();
    const scalar pz = position_.z();

    if (activePositionsPtr_ != nullptr)
    {
        if (influenceCells_.empty())
        {
            forAll(*activePositionsPtr_, cellI)
            {
                const vector& c = (*activePositionsPtr_)[cellI];

                scalar dx = c.x() - px;
                scalar dy = c.y() - py;
                scalar dz = c.z() - pz;

                scalar dis = dx*dx + dy*dy + dz*dz;
                if (dis <= sphereRadiusSqr)
                {
                    (*activeForceFieldPtr_)[cellI] +=
                        scaledForce*Foam::exp(-dis*invepsilonSqr);
                }
            }
        }
        else
        {
            const List<label>& cells = influenceCells_[azimuthIndex_];

            forAll(cells, i)
            {
                label cellI = cells[i];
                const vector& c = (*activePositionsPtr_)[cellI];

                scalar dx = c.x() - px;
                scalar dy = c.y() - py;
                scalar dz = c.z() - pz;

                scalar dis = dx*dx + dy*dy + dz*dz;
                if (dis <= sphereRadiusSqr)
                {
                    (*activeForceFieldPtr_)[cellI] +=
                        scaledForce*Foam::exp(-dis*invepsilonSqr);
                }
            }
        }
    }
    else
    {
        forAll(mesh_.cells(), cellI)
        {
            const vector& c = C[cellI];

            scalar dx = c.x() - px;
            scalar dy = c.y() - py;
            scalar dz = c.z() - pz;

            scalar dis = dx*dx + dy*dy + dz*dz;
            if (dis <= sphereRadiusSqr)
            {
                force[cellI] += scaledForce*Foam::exp(-dis*invepsilonSqr);
            }
        }
    }

    if (debug)
    {
        Info<< "    sphereRadius: " << sphereRadius << endl;
    }
}


void Foam::fv::actuatorLineElement::allocateInfluenceCells
(
    label count,
    bool cacheInteractions
)
{
    tmpInfluenceCells_.setSize(8000);
    if (cacheInteractions)
    {
        influenceCells_.setSize(count);
    }
    centerCellI_.setSize(count, -2);
    centerProcI_.setSize(count, -2);
    if (velocitySampleRadius_ > 0.0)
    {
        ringCellI_.setSize(count);
        ringProcI_.setSize(count);
        forAll(ringCellI_, azimuthI)
        {
            ringCellI_[azimuthI].setSize
            (
                nVelocitySamples_,
                -2
            );

            ringProcI_[azimuthI].setSize
            (
                nVelocitySamples_,
                -2
            );
        }
    }
}


void Foam::fv::actuatorLineElement::constructInfluenceCellList
(
    label azimuthIndex,
    labelList& globalToLocal,
    label& nActive
)
{
    // Calculate projection width
    dragCoefficient_ = 4; // Should realistically not be larger than this
    scalar epsilon = calcProjectionEpsilon();
    scalar projectionRadius = (epsilon*Foam::sqrt(Foam::log(1.0/0.001)));

    // Apply force to the cells within the element's sphere of influence
    scalar sphereRadius = chordLength_ + projectionRadius;
    scalar sphereRadiusSqr = sphereRadius*sphereRadius;

    const vectorField& C = mesh_.C();

    label nCells = 0;

    forAll(C, cellI)
    {
        scalar dis = magSqr(C[cellI] - position_);
        if (dis <= sphereRadiusSqr)
        {
            // Grow buffer if needed
            if (nCells >= tmpInfluenceCells_.size())
            {
                tmpInfluenceCells_.setSize
                (
                    2*tmpInfluenceCells_.size()
                );
            }
            if (globalToLocal[cellI] == -1)
            {
                globalToLocal[cellI] = nActive++;
            }
            if (influenceCells_.size() > 0)
            {
                tmpInfluenceCells_[nCells++] = globalToLocal[cellI];
            }
        }
    }
    if (influenceCells_.size() > 0)
    {
        influenceCells_[azimuthIndex].setSize(nCells);

        for (label i = 0; i < nCells; i++)
        {
            influenceCells_[azimuthIndex][i] =
                tmpInfluenceCells_[i];
        }
    }
}

void Foam::fv::actuatorLineElement::setAzimuthIndex
(
    label azimuthIndex
)
{
    azimuthIndex_ = azimuthIndex;
    if (azimuthIndex == 0)
    {
        stringBuffer_.str("");
        stringBuffer_.clear();
    }
}

void Foam::fv::actuatorLineElement::calculateInflowVelocity
(
    const interpolationCellPoint<vector>& UInterp
)
{
    // Find local flow velocity by interpolating to element location
    vector inflowVelocityPoint = position_;
    
    vector localVelocitySum = vector::zero;
    label localNSamples = 0;
    
    // If the flow only is sampled in the center
    if (velocitySampleRadius_ <= 0.0)
    {
        label inflowCellI = findCell(inflowVelocityPoint);
        if (inflowCellI >= 0)
        {
            localVelocitySum = UInterp.interpolate
            (
                inflowVelocityPoint,
                inflowCellI
            );
            localNSamples = 1;
        }
    }
    // If the flow is sampled by using a circle around position_
    else
    {
        // Circle radius should be normalized with epsilon
        scalar sampleRadius = calcProjectionEpsilon()*velocitySampleRadius_;

        // Unit vector in chordwise direction
        vector chordNormal = chordDirection_ / mag(chordDirection_);
        
        // Calculate mean value over all circle points
        for (label point = 0; point < nVelocitySamples_; point++)
        {
            // distribute the points evenly in terms of angular distance
            scalar pointAngle = Foam::constant::mathematical::pi * 2.0 * point/
                                nVelocitySamples_;
            scalar chordDist = sampleRadius * Foam::cos(pointAngle);
            scalar normalDist = sampleRadius * Foam::sin(pointAngle);
            vector samplePoint = inflowVelocityPoint +
                                 chordDist * chordNormal +
                                 normalDist * planformNormal_;

            // Sample the velocity
            activeRingIndex_ = point; // actuator disc, use ring cache
            label sampleCellI = findCell(samplePoint);
            if (sampleCellI >= 0)
            {
                localVelocitySum += UInterp.interpolate
                (
                    samplePoint,
                    sampleCellI
                );
                localNSamples++;
            }
        }
        // no longer doing ring sampling, do not use this cache
        activeRingIndex_ = -1;
    }
    
    // Reduce inflow velocity over all processors
    reduce(localVelocitySum, sumOp<vector>());
    reduce(localNSamples, sumOp<label>());

    // Set inflow Velocity as the mean value
    inflowVelocity_ = 1.0 / localNSamples * localVelocitySum;

    // If inflow velocity is not detected, position is not in the mesh
    if (localNSamples == 0)
    {
        // Raise fatal error since inflow velocity cannot be detected
        FatalErrorIn("void actuatorLineElement::calculateForce()")
            << "Inflow velocity point for " << name_ << " Position: "
            << position_ << " not found in mesh"
            << abort(FatalError);
    }
}


void Foam::fv::actuatorLineElement::createOutputFile()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = mesh_.time().path()/"../postProcessing/actuatorLineElements"
            / mesh_.time().timeName();
    }
    else
    {
        dir = mesh_.time().path()/"postProcessing/actuatorLineElements"
            / mesh_.time().timeName();
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_.open(dir/name_ + ".csv", std::ios::out);

    if (outputFile_.is_open())
    {
        outputFile_<< "time,root_dist,x,y,z,rel_vel_mag,Re,alpha_deg,"
                   << "alpha_geom_deg,cl,cd,fx,fy,fz,end_effect_factor,"
                   << "c_ref_t,c_ref_n,f_ref_t,f_ref_n" << std::endl;
    }
}


void Foam::fv::actuatorLineElement::writePerf()
{
    scalar time = mesh_.time().value();

    // write time,root_dist,x,y,z,rel_vel_mag,Re,alpha_deg,alpha_geom_deg,cl,cd,
    // fx,fy,fz,end_effect_factor,c_ref_t,c_ref_n,f_ref_t,f_ref_n
    stringBuffer_<< time << "," << rootDistance_ << "," << position_.x() << ","
                << position_.y() << "," << position_.z() << ","
                << mag(relativeVelocity_) << "," << Re_ << "," << angleOfAttack_
                << "," << angleOfAttackGeom_ << "," << liftCoefficient_ << ","
                << dragCoefficient_ << "," << forceVector_.x() << ","
                << forceVector_.y() << "," << forceVector_.z() << ","
                << endEffectFactor_ << "," << tangentialRefCoefficient() << ","
                << normalRefCoefficient() << "," << tangentialRefForce() << ","
                << normalRefForce() << std::endl;

    // only write to file with writePerf_, writePerfEnd_ writes in destructor
    if (writePerf_ && outputFile_.is_open())
    {
        outputFile_ << stringBuffer_.str();
        stringBuffer_.str("");
        stringBuffer_.clear();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineElement::actuatorLineElement
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    name_(name),
    mesh_(mesh),
    meshBoundBox_(mesh_.points(), false),
    planformNormal_(vector::zero),
    velocity_(vector::zero),
    forceVector_(vector::zero),
    previousCellI_(-1),
    previousProcI_(-1),
    previousLocation_(point::zero),
    previousLocationValid_(false),
    centerCellI_(0),
    centerProcI_(0),
    ringCellI_(0),
    ringProcI_(0),
    activeRingIndex_(-1),
    relativeVelocity_(vector::zero),
    relativeVelocityGeom_(vector::zero),
    angleOfAttack_(0.0),
    angleOfAttackGeom_(0.0),
    liftCoefficient_(0.0),
    dragCoefficient_(0.0),
    momentCoefficient_(0.0),
    profileName_(dict.lookup("profileName")),
    profileData_(profileName_, dict.subDict("profileData"), debug),
    dynamicStallActive_(false),
    omega_(0.0),
    chordMount_(0.25),
    flowCurvatureActive_(false),
    flowCurvatureModelName_("none"),
    velocityLE_(vector::zero),
    velocityTE_(vector::zero),
    writePerf_(false),
    writePerfEnd_(false),
    outputFile_(nullptr),
    rootDistance_(0.0),
    endEffectFactor_(1.0),
    addedMassActive_(dict.lookupOrDefault("addedMass", false)),
    addedMass_(mesh.time(), dict.lookupOrDefault("chordLength", 1.0), debug),
    influenceCells_(0),
    tmpInfluenceCells_(0),
    activePositionsPtr_(nullptr),
    activeForceFieldPtr_(nullptr),
    azimuthIndex_(0)
{
    meshBoundBox_.inflate(1e-6);
    read();
    if (writePerf_ || writePerfEnd_)
    {
        createOutputFile();
    }
    mesh_.cellCells();
    int precision = 6;
    if (mesh_.time().controlDict().found("writePrecision"))
    {
        mesh_.time().controlDict().lookup("writePrecision") >> precision;
    }
    stringBuffer_.precision(precision);
}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineElement::~actuatorLineElement()
{
    if (writePerfEnd_ && writePerf_ == false)
    {
        if (outputFile_.is_open())
        {
            outputFile_ << stringBuffer_.str();
        }
    }
    if (outputFile_.is_open())
    {
       outputFile_.close();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::fv::actuatorLineElement::name() const
{
    return name_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::chordLength() const
{
    return chordLength_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::spanLength()
{
    return spanLength_;
}


const Foam::vector& Foam::fv::actuatorLineElement::position()
{
    return position_;
}


const Foam::vector& Foam::fv::actuatorLineElement::velocity()
{
    return velocity_;
}


const Foam::vector& Foam::fv::actuatorLineElement::relativeVelocity()
{
    return relativeVelocity_;
}


const Foam::vector& Foam::fv::actuatorLineElement::relativeVelocityGeom()
{
    return relativeVelocityGeom_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::angleOfAttack()
{
    return angleOfAttack_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::angleOfAttackGeom()
{
    return angleOfAttackGeom_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::liftCoefficient()
{
    return liftCoefficient_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::dragCoefficient()
{
    return dragCoefficient_;
}


const Foam::scalar& Foam::fv::actuatorLineElement::momentCoefficient()
{
    return momentCoefficient_;
}


Foam::scalar Foam::fv::actuatorLineElement::tangentialRefCoefficient()
{
    return profileData_.convertToCRT
    (
        liftCoefficient_,
        dragCoefficient_,
        inflowRefAngle()
    );
}


Foam::scalar Foam::fv::actuatorLineElement::tangentialRefForce()
{
    return 0.5 * chordLength_ * tangentialRefCoefficient()
        * magSqr(relativeVelocity_);
}


Foam::scalar Foam::fv::actuatorLineElement::normalRefCoefficient()
{
    return profileData_.convertToCRN
    (
        liftCoefficient_,
        dragCoefficient_,
        inflowRefAngle()
    );
}


Foam::scalar Foam::fv::actuatorLineElement::normalRefForce()
{
    return 0.5 * chordLength_ * normalRefCoefficient()
        * magSqr(relativeVelocity_);
}


Foam::scalar Foam::fv::actuatorLineElement::inflowRefAngle()
{
    // Calculate inflow velocity angle in degrees (AFTAL Phi)
    scalar arg =
        (-relativeVelocity_ & chordRefDirection_)
        / (mag(relativeVelocity_) * mag(chordRefDirection_));
    scalar inflowVelAngleRad =
        acos(sign(arg)*min(Foam::scalar(1.0), mag(arg)));
    return radToDeg(inflowVelAngleRad);
}


const Foam::scalar& Foam::fv::actuatorLineElement::rootDistance()
{
    return rootDistance_;
}


void Foam::fv::actuatorLineElement::calculateForce
(
    const interpolationCellPoint<vector>& UInterp
)
{
    scalar pi = Foam::constant::mathematical::pi;

    // Calculate vector normal to chord--span plane
    planformNormal_ = -chordDirection_ ^ spanDirection_;
    planformNormal_ /= mag(planformNormal_);

    if (debug)
    {
        Info<< "Calculating force contribution from actuatorLineElement "
            << name_ << endl;
        Info<< "    position: " << position_ << endl;
        Info<< "    chordDirection: " << chordDirection_ << endl;
        Info<< "    spanDirection: " << spanDirection_ << endl;
        Info<< "    elementVelocity: " << velocity_ << endl;
        Info<< "    planformNormal: " << planformNormal_ << endl;
    }

    // Find local flow velocity by interpolating to element location
    calculateInflowVelocity(UInterp);
    

    // Subtract spanwise component of inflow velocity
    vector spanwiseVelocity = spanDirection_
                            * (inflowVelocity_ & spanDirection_)
                            / magSqr(spanDirection_);
    inflowVelocity_ -= spanwiseVelocity;

    // Calculate relative velocity and Reynolds number
    relativeVelocity_ = inflowVelocity_ - velocity_;
    Re_ = mag(relativeVelocity_)*chordLength_/nu_;

    // Calculate angle of attack (radians)
    scalar arg =
        (planformNormal_ & relativeVelocity_)
        / (mag(planformNormal_) * mag(relativeVelocity_));
    scalar angleOfAttackRad = asin(sign(arg)*min(Foam::scalar(1.0), mag(arg)));
    scalar angleOfAttackUncorrected = radToDeg(angleOfAttackRad);
    relativeVelocityGeom_ = freeStreamVelocity_ - velocity_;
    scalar argGeom =
        (planformNormal_ & relativeVelocityGeom_)
        / (mag(planformNormal_) * mag(relativeVelocityGeom_));
    angleOfAttackGeom_ =
        asin(sign(argGeom)*min(Foam::scalar(1.0), mag(argGeom)));
    angleOfAttackGeom_ *= 180.0/pi;

    // Apply flow curvature correction to angle of attack
    if (flowCurvatureActive_)
    {
        correctFlowCurvature(angleOfAttackRad);
    }

    // Calculate angle of attack in degrees
    angleOfAttack_ = radToDeg(angleOfAttackRad);

    // Update Reynolds number of profile data
    
    //calculateInflowVelocity(UInterp); // Is this result ever used?

    profileData_.updateRe(Re_);

    // Lookup lift and drag coefficients
    lookupCoefficients();

    if (debug)
    {
        Info<< "    inflowVelocity: " << inflowVelocity_ << endl;
        Info<< "    relativeVelocity: " << relativeVelocity_ << endl;
        Info<< "    Reynolds number: " << Re_ << endl;
        Info<< "    Geometric angle of attack (degrees): "
            << angleOfAttackGeom_ << endl;
        Info<< "    Angle of attack (uncorrected, degrees): "
            << angleOfAttackUncorrected << endl;
        Info<< "    Angle of attack (corrected, degrees): "
            << angleOfAttack_ << endl;
    }

    // Correct coefficients with dynamic stall model
    if (dynamicStallActive_)
    {
        dynamicStall_->correct
        (
            mag(relativeVelocity_),
            angleOfAttack_,
            liftCoefficient_,
            dragCoefficient_,
            momentCoefficient_
        );
    }

    // Correct for added mass effects
    if (addedMassActive_)
    {
        addedMass_.correct
        (
            liftCoefficient_,
            dragCoefficient_,
            momentCoefficient_,
            degToRad(angleOfAttack_),
            mag(chordDirection_ & relativeVelocity_),
            mag(planformNormal_ & relativeVelocity_)
        );
    }

    // Apply end effect correction factor to lift coefficient
    liftCoefficient_ *= endEffectFactor_;

    // Calculate force per unit density
    scalar area = chordLength_ * spanLength_;
    scalar magSqrU = magSqr(relativeVelocity_);
    scalar lift = 0.5*area*liftCoefficient_*magSqrU;
    scalar drag = 0.5*area*dragCoefficient_*magSqrU;
    vector liftDirection = relativeVelocity_ ^ spanDirection_;
    liftDirection /= mag(liftDirection);
    vector dragDirection = relativeVelocity_/mag(relativeVelocity_);
    forceVector_ = lift*liftDirection + drag*dragDirection;

    if (debug)
    {
        Info<< "    liftDirection: " << liftDirection << endl;
        Info<< "    dragDirection: " << dragDirection << endl;
        Info<< "    force (per unit density): " << forceVector_ << endl;
    }
}


void Foam::fv::actuatorLineElement::rotate
(
    vector rotationPoint,
    vector axis,
    scalar radians,
    bool rotateVelocity=true
)
{
    // Declare and define the rotation matrix (from SOWFA)
    tensor RM;
    scalar angle = radians;
    RM.xx() = Foam::sqr(axis.x())
            + (1.0 - Foam::sqr(axis.x())) * Foam::cos(angle);
    RM.xy() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) - axis.z() * Foam::sin(angle);
    RM.xz() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.y() * Foam::sin(angle);
    RM.yx() = axis.x() * axis.y()
            * (1.0 - Foam::cos(angle)) + axis.z() * Foam::sin(angle);
    RM.yy() = Foam::sqr(axis.y())
            + (1.0 - Foam::sqr(axis.y())) * Foam::cos(angle);
    RM.yz() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.x() * Foam::sin(angle);
    RM.zx() = axis.x() * axis.z()
            * (1.0 - Foam::cos(angle)) - axis.y() * Foam::sin(angle);
    RM.zy() = axis.y() * axis.z()
            * (1.0 - Foam::cos(angle)) + axis.x() * Foam::sin(angle);
    RM.zz() = Foam::sqr(axis.z())
            + (1.0 - Foam::sqr(axis.z())) * Foam::cos(angle);

    if (debug)
    {
        Info<< "Rotating actuatorLineElement: " << name_ << endl;
        Info<< "Rotation point: " << rotationPoint << endl;
        Info<< "Rotation axis: " << axis << endl;
        Info<< "Rotation angle (radians): " << radians << endl;
        Info<< "Rotation matrix:" << endl << RM << endl;
        Info<< "Initial position: " << position_ << endl;
        Info<< "Initial chordDirection: " << chordDirection_ << endl;
        Info<< "Initial spanDirection: " << spanDirection_ << endl;
        Info<< "Initial velocity: " << velocity_ << endl;
    }

    // Rotation matrices make a rotation about the origin, so need to subtract
    // rotation point off the point to be rotated.
    vector point = position_;
    point -= rotationPoint;

    // Perform the rotation.
    point = RM & point;

    // Return the rotated point to its new location relative to the rotation
    // point
    point += rotationPoint;

    // Set the position of the element
    position_ = point;

    // Rotate the span and chord vectors of the element
    chordDirection_ = RM & chordDirection_;
    spanDirection_ = RM & spanDirection_;

    // Rotate the element's velocity vector if specified
    if (rotateVelocity)
    {
        velocity_ = RM & velocity_;
        chordRefDirection_ = RM & chordRefDirection_;
    }

    if (debug)
    {
        Info<< "Final position: " << position_ << endl;
        Info<< "Final chordDirection: " << chordDirection_ << endl;
        Info<< "Final chordRefDirection: " << chordRefDirection_ << endl;
        Info<< "Final spanDirection: " << spanDirection_ << endl;
        Info<< "Final velocity: " << velocity_ << endl << endl;
    }
}


void Foam::fv::actuatorLineElement::pitch
(
    scalar radians,
    scalar chordFraction
)
{
    vector rotationPoint = position_;
    rotationPoint += chordDirection_*(chordMount_ - chordFraction);
    rotate(rotationPoint, spanDirection_, radians, false);
}


void Foam::fv::actuatorLineElement::translate(vector translationVector)
{
    position_ += translationVector;
}


void Foam::fv::actuatorLineElement::setVelocity(vector velocity)
{
    if (debug)
    {
        Info<< "Changing velocity of " << name_ << " from "
            << velocity_ << " to " << velocity << endl << endl;
    }
    velocity_ = velocity;
}


void Foam::fv::actuatorLineElement::setSpeed(scalar speed)
{
    if (mag(velocity_) > 0)
    {
        velocity_ /= mag(velocity_);
        velocity_ *= speed;
    }
}


void Foam::fv::actuatorLineElement::setSpeed
(
    vector point,
    vector axis,
    scalar omega
)
{
    if (debug)
    {
        Info<< "Setting speed of " << name_ << " from rotation" << endl;
        Info<< "    Initial velocity: " << velocity_ << endl;
    }

    // First find radius from axis to element position -- formula from
    // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    vector point2 = point + axis;
    scalar radius = mag((position_ - point) ^ (position_ - point2))
                  / mag(point2 - point);
    scalar speed = omega*radius;
    setSpeed(speed);

    scalar angleLE = 0.0;
    scalar angleTE = 0.0;
    if (radius > 0.0)
    {
        // Set velocity at leading edge
        scalar radiusLE = sqrt(magSqr(0.25*chordLength_) + magSqr(radius));
        angleLE = atan2(0.25*chordLength_, radius);
        velocityLE_ = velocity_*radiusLE/radius;
        rotateVector(velocityLE_, vector::zero, spanDirection_, angleLE);

        // Set velocity at trailing edge
        scalar radiusTE = sqrt(magSqr(0.75*chordLength_) + magSqr(radius));
        angleTE = atan2(-0.75*chordLength_, radius);
        velocityTE_ = velocity_*radiusTE/radius;
        rotateVector(velocityTE_, vector::zero, spanDirection_, angleTE);
    }

    // Also set omega for flow curvature correction
    setOmega(omega);

    if (debug)
    {
        Info<< "    Radius: " << radius << endl;
        Info<< "    Final velocity: " << velocity_ << endl;
        Info<< "    Leading edge velocity: " << velocityLE_ << endl;
        Info<< "    Trailing edge velocity: " << velocityTE_ << endl;
        Info<< "    Leading edge velocity angle (radians): "
            << angleLE << endl;
        Info<< "    Trailing edge velocity angle (radians): "
            << angleTE << endl;
    }
}


void Foam::fv::actuatorLineElement::scaleVelocity(scalar scale)
{
    velocity_ *= scale;
}


const Foam::vector& Foam::fv::actuatorLineElement::force()
{
    return forceVector_;
}


Foam::vector Foam::fv::actuatorLineElement::moment(vector point)
{
    // Calculate radius vector
    vector radius = position_ - point;
    vector moment = radius ^ forceVector_;
    vector pitchingMoment = 0.5*chordLength_*chordLength_*spanLength_
                          * momentCoefficient_*magSqr(relativeVelocity_)
                          * spanDirection_;
    return moment + pitchingMoment;
}


void Foam::fv::actuatorLineElement::addForce
(
    const interpolationCellPoint<vector>& UInterp,
    volVectorField& forceField,
    scalar scale
)
{
    /*volVectorField forceFieldI
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            forceField.dimensions(),
            vector::zero
        )
    );*/

    //const volVectorField& Uin(eqn.psi());

    calculateForce(UInterp);
    applyForceField(forceField, scale);

    // Add force to total actuator line force
    //forceField += forceFieldI;

    // Write performance to file
    if (Pstream::master() && (writePerf_ || writePerfEnd_))
    {
        writePerf();
    }
}


void Foam::fv::actuatorLineElement::addForce
(
    const volScalarField& rho,
    const interpolationCellPoint<vector>& UInterp,
    volVectorField& forceField,
    scalar scale
)
{
    /*volVectorField forceFieldI
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            forceField.dimensions()/rho.dimensions(),
            vector::zero
        )
    );*/

    //const volVectorField& Uin(eqn.psi());
    calculateForce(UInterp);
    applyForceField(forceField, scale);

    // Multiply force vector by local density
    multiplyForceRho(rho);

    // Multiply this element's force field by density field
    // (moved to actuatorLineSource as it is the same rho for all)
    //forceFieldI *= rho;

    // Add force to total actuator line force
    // forceField += forceFieldI;

    // Write performance to file
    if (Pstream::master() && (writePerf_ || writePerfEnd_))
    {
        writePerf();
    }
}

void Foam::fv::actuatorLineElement::addTurbulence
(
    fvMatrix<scalar>& eqn,
    word fieldName
)
{
    // Calculate projection radius
    scalar epsilon = calcProjectionEpsilon();
    scalar projectionRadius = (epsilon*Foam::sqrt(Foam::log(1.0/0.001)));

    // Add turbulence to the cells within the element's sphere of influence
    scalar sphereRadius = chordLength_ + projectionRadius;
    scalar sphereRadiusSqr = sphereRadius*sphereRadius;
    scalar invepsilonSqr = 1.0/(epsilon*epsilon);
    scalar internalFactor = 1.0/(Foam::pow(epsilon, 3)
                          * Foam::pow(Foam::constant::mathematical::pi, 1.5));

    const volScalarField& kField =
        mesh_.lookupObject<volScalarField>("k");
    const volScalarField& epsilonField =
        mesh_.lookupObject<volScalarField>("epsilon");

    // Calculate TKE injection rate
    scalar area = chordLength_ * spanLength_;
    scalar magSqrU = magSqr(relativeVelocity_);
    scalar lift = 0.5*area*liftCoefficient_*magSqrU;
    scalar drag = 0.5*area*dragCoefficient_*magSqrU;
    scalar Ceps1 = 1.44;
    scalar magU = mag(relativeVelocity_);
    scalar liftBasedTurbInjection = -1.0;
    scalar dragBasedTurbInjection = -1.0;
    scalar turbRed =
      liftBasedTurbInjection * lift * magU +
      dragBasedTurbInjection * drag * magU;
    
    const vectorField& C = mesh_.C();
    scalarField& src = eqn.source();
    const scalarField& V = mesh_.V();

    if (influenceCells_.empty())
    {
        forAll(mesh_.cells(), cellI)
        {
            scalar dis = magSqr(C[cellI] - position_);
            if (dis <= sphereRadiusSqr)
            {
                scalar factor = Foam::exp(-dis*invepsilonSqr)*internalFactor;

                if (fieldName == "k")
                {
		    src[cellI] -= factor * turbRed * V[cellI];
                }
                else if (fieldName == "epsilon")
                {
		    scalar kCell   = kField[cellI];
                    scalar epsCell = epsilonField[cellI];
		    src[cellI] -= factor * Ceps1 * turbRed
		                  * (epsCell/kCell) * V[cellI];
                }
            }
        }
    }
    else
    {
        const labelList& cells = influenceCells_[azimuthIndex_];
        forAll(cells, i)
        {
            label cellI = cells[i];
            scalar dis = magSqr(C[cellI] - position_);
            if (dis <= sphereRadiusSqr)
            {
                scalar factor = Foam::exp(-dis*invepsilonSqr)*internalFactor;
                if (fieldName == "k")
                {
		    src[cellI] -= factor * turbRed * V[cellI];
                }
                else if (fieldName == "epsilon")
                {
		    scalar kCell   = kField[cellI];
                    scalar epsCell = epsilonField[cellI];
		    src[cellI] -= factor * Ceps1 * turbRed
		                  * (epsCell/kCell) * V[cellI];
                }
            }
        }
    }
}


void Foam::fv::actuatorLineElement::setDynamicStallActive(bool active)
{
    dynamicStallActive_ = active;
}


void Foam::fv::actuatorLineElement::setOmega(scalar omega)
{
    omega_ = omega;
}


void Foam::fv::actuatorLineElement::setEndEffectFactor(scalar factor)
{
    endEffectFactor_ = factor;
}


void Foam::fv::actuatorLineElement::setVelocitySampleRadius(scalar radius)
{
    velocitySampleRadius_ = radius;
}


void Foam::fv::actuatorLineElement::setNVelocitySamples(label nSamples)
{
    nVelocitySamples_ = nSamples;
}


void Foam::fv::actuatorLineElement::setCustomTime
(
    scalar time,
    scalar deltaT,
    bool useCustomTime
)
{
    if (dynamicStallActive_)
    {
        dynamicStall_->setCustomTime(time, deltaT, useCustomTime);
    }
    if (addedMassActive_)
    {
        addedMass_.setCustomTime(time, deltaT, useCustomTime);
    }
}


void Foam::fv::actuatorLineElement::setCompactFields
(
    vectorField& activePositions,
    vectorField& activeForceField
)
{
    activePositionsPtr_ = &activePositions;
    activeForceFieldPtr_ = &activeForceField;
}
// ************************************************************************* //
