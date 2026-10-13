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
    dict_.lookup("position") >> position_[0];
    dict_.lookup("chordLength") >> chordLength_;
    dict_.lookup("chordDirection") >> chordDirection_[0];
    dict_.lookup("chordRefDirection") >> chordRefDirection_[0];
    dict_.lookup("chordMount") >> chordMount_;
    dict_.lookup("spanLength") >> spanLength_;
    dict_.lookup("spanDirection") >> spanDirection_[0];
    dict_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
    freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
    dict_.lookup("rootDistance") >> rootDistance_;
    dict_.lookup("cone") >> cone_;
    dict_.lookup("velocitySampleRadius") >> velocitySampleRadius_;
    dict_.lookup("nVelocitySamples") >> nVelocitySamples_;
    dict_.lookup("endEffectsCorrectAoA") >> endEffectsCorrectAngleOfAttack_;
    dict_.lookup("endEffectsCorrectWindAngle") >> endEffectsCorrectAngleOfRelWind_;


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
    if (mesh_.foundObject<IOdictionary>("transportProperties"))
    {
        const dictionary& transportProperties = mesh_.lookupObject<IOdictionary>
        (
            "transportProperties"
        );
        dimensionedScalar nu;
        transportProperties.lookup("nu") >> nu;
        nu_ = nu.value();
    }
    else if (mesh_.foundObject<volScalarField>("thermo:mu"))
    {
        // get the dynamic viscosity and density fields
        const volScalarField& mu =
            mesh_.lookupObject<volScalarField>("thermo:mu");
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>("rho");

        // for simplicity, assume that nu is approximately constant
        // and use value from first cell
        nu_ = mu[0] / rho[0];
    }
    else
    {
        FatalErrorIn("actuatorLineElement::read()")
                << "Could not find transportProperties,"
                << " nor thermophysicalProperties in simulation"
                << abort(FatalError);
    }

    // Read writePerf switch
    dict_.lookup("writePerf") >> writePerf_;
    dict_.lookup("writePerfEnd") >> writePerfEnd_;

    if (debug)
    {
        Info<< "actuatorLineElement properties:" << endl;
        Info<< "Position: " << position_[azimuthIndex_] << endl;
        Info<< "chordLength: " << chordLength_ << endl;
        Info<< "chordDirection: " << chordDirection_[azimuthIndex_] << endl;
        Info<< "spanLength: " << spanLength_ << endl;
        Info<< "spanDirection: " << spanDirection_[azimuthIndex_] << endl;
        Info<< "cone: " << cone_ << endl;
        Info<< "writePerf: " << writePerf_ << endl;
        Info<< "writePerfEnd: " << writePerfEnd_ << endl;
    }

    // Lookup Gaussian coeffs from profileData dict if present
    dictionary GaussianCoeffs = profileData_.dict().subOrEmptyDict
    (
        "GaussianCoeffs"
    );
    chordFactor_ = GaussianCoeffs.lookupOrDefault("chordFactor", 0.25);
    dragFactor_ = GaussianCoeffs.lookupOrDefault("dragFactor", 1.0);
    meshFactor_ = GaussianCoeffs.lookupOrDefault("meshFactor", 2.0);

    // Make it pssible to override these factors for an individual line
    if (dict_.found("chordFactor"))
    {
        dict_.lookup("chordFactor") >> chordFactor_;
    }
    if (dict_.found("dragFactor"))
    {
        dict_.lookup("dragFactor") >> dragFactor_;
    }
    if (dict_.found("meshFactor"))
    {
        dict_.lookup("meshFactor") >> meshFactor_;
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

void Foam::fv::actuatorLineElement::lookupCoefficients()
{
    liftCoefficient_[azimuthIndex_] =
        profileData_.liftCoefficient(angleOfAttack_);
    dragCoefficient_[azimuthIndex_] =
        profileData_.dragCoefficient(angleOfAttack_);
    momentCoefficient_[azimuthIndex_] =
        profileData_.momentCoefficient(angleOfAttack_);
}


void Foam::fv::actuatorLineElement::calcProjectionEpsilon()
{
    // Provide ideal epsilon target for lift based on chord length
    scalar epsilonLift = chordFactor_*chordLength_;

    // Epsilon based on drag/momentum thickness
    scalar epsilonDrag =
        dragFactor_*dragCoefficient_[azimuthIndex_]*chordLength_/2.0;

    // Threshold is based on lift or drag, whichever is larger
    scalar epsilonThreshold = Foam::max(epsilonLift, epsilonDrag);

    epsilon_[azimuthIndex_] = VGREAT;
    scalar epsilonMesh = VGREAT;
    const scalarField& V = mesh_.V();

    if (centerProcI_[azimuthIndex_] == Pstream::myProcNo()
            && centerCellI_[azimuthIndex_] >= 0)
    {
        // Projection width based on local cell size (from Troldborg (2008))
        epsilonMesh = 2.0*Foam::cbrt(V[centerCellI_[azimuthIndex_]]);
        epsilonMesh *= meshFactor_; // Cell could have non-unity aspect ratio

        if (epsilonMesh > epsilonThreshold)
        {
            epsilon_[azimuthIndex_] = epsilonMesh;
        }
        else
        {
            epsilon_[azimuthIndex_] = epsilonThreshold;
        }
    }

    if (debug)
    {
        //reduce(epsilonMesh, minOp<scalar>());
        word epsilonMethod;
        if (epsilon_[azimuthIndex_] == epsilonLift)
        {
            epsilonMethod = "lift-based";
        }
        else if (epsilon_[azimuthIndex_] == epsilonDrag)
        {
            epsilonMethod = "drag-based";
        }
        else if (epsilon_[azimuthIndex_] == epsilonMesh)
        {
            epsilonMethod = "mesh-based";
        }
        Info<< "    epsilon (" << epsilonMethod
            << "): " << epsilon_[azimuthIndex_] << endl;
    }
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
            omega_*chordLength_/
            (2*mag(relativeVelocity_[azimuthIndex_]))*cos(cone_);
    }
    else if (flowCurvatureModelName_ == "MandalBurton")
    {
        // Calculate relative velocity at leading and trailing edge
        vector relativeVelocityLE =
            inflowVelocity_[azimuthIndex_] - velocityLE_;
        vector relativeVelocityTE =
            inflowVelocity_[azimuthIndex_] - velocityTE_;

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


void Foam::fv::actuatorLineElement::multiplyForceRho()
{
    forceVector_[azimuthIndex_] *= localRho_[azimuthIndex_];
}


void Foam::fv::actuatorLineElement::applyForceField
(
    volVectorField& forceField,
    scalar scale
)
{
    // Calculate projection width
    //scalar epsilon[azimuthIndex_] = calcProjectionEpsilon();

    scalar projectionRadius =
        (epsilon_[azimuthIndex_]*Foam::sqrt(Foam::log(1.0/0.001)));

    // Apply force to the cells within the element's sphere of influence
    scalar sphereRadius = chordLength_ + projectionRadius;
    scalar sphereRadiusSqr = sphereRadius*sphereRadius;
    scalar invepsilonSqr =
        1.0/(epsilon_[azimuthIndex_]*epsilon_[azimuthIndex_]);
    scalar internalFactor = scale/(Foam::pow(epsilon_[azimuthIndex_], 3)
                          * Foam::pow(Foam::constant::mathematical::pi, 1.5));

    // forceField is opposite forceVector
    const vector scaledForce = -forceVector_[azimuthIndex_]*internalFactor;

    const vectorField& C = mesh_.C();

    vectorField& force =
        forceField.primitiveFieldRef();

    const scalar px = position_[azimuthIndex_].x();
    const scalar py = position_[azimuthIndex_].y();
    const scalar pz = position_[azimuthIndex_].z();

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
        // Check if the sphere ever will be within the mesh
        // (May not be the case for parallel runs)
        scalar distSqr = 0.0;

        for (direction dir=0; dir<3; dir++)
        {
            if (position_[azimuthIndex_][dir] < meshBoundBox_.min()[dir])
            {
                scalar d = meshBoundBox_.min()[dir] -
                           position_[azimuthIndex_][dir];
                distSqr += d*d;
            }
            else if (position_[azimuthIndex_][dir] > meshBoundBox_.max()[dir])
            {
                scalar d = position_[azimuthIndex_][dir] -
                           meshBoundBox_.max()[dir];
                distSqr += d*d;
            }
        }
        // Only run the loop if sphere may overlap with mesh
        if (distSqr <= sphereRadiusSqr)
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
    centerCellI_.setSize(count);
    centerCellI_ = -1;

    centerProcI_.setSize(count);
    centerProcI_ = -1;

    chordDirection_.setSize(count);
    chordDirection_ = chordDirection_[0];

    spanDirection_.setSize(count);
    spanDirection_ = spanDirection_[0];

    chordRefDirection_.setSize(count);
    chordRefDirection_ = chordRefDirection_[0];

    localRho_.setSize(count);
    localRho_ = localRho_[0];

    localMu_.setSize(count);
    localMu_ = localMu_[0];

    position_.setSize(count);
    position_ = position_[0];

    previousLocation_.setSize(count);
    previousLocation_ = previousLocation_[0];

    velocity_.setSize(count);
    velocity_ = velocity_[0];

    forceVector_.setSize(count);
    forceVector_ = forceVector_[0];

    inflowVelocity_.setSize(count);
    inflowVelocity_ = inflowVelocity_[0];

    epsilon_.setSize(count);
    epsilon_ = chordLength_;

    relativeVelocity_.setSize(count);
    relativeVelocity_ = relativeVelocity_[0];

    liftCoefficient_.setSize(count);
    liftCoefficient_ = liftCoefficient_[0];

    dragCoefficient_.setSize(count);
    dragCoefficient_ = dragCoefficient_[0];

    momentCoefficient_.setSize(count);
    momentCoefficient_ = momentCoefficient_[0];

    if (velocitySampleRadius_ > 0.0)
    {
        ringCellI_.setSize(count);
        ringProcI_.setSize(count);
        velocitiesRing_.setSize(count);
        previousRingLocation_.setSize(count);
        forAll(ringCellI_, azimuthI)
        {
            ringCellI_[azimuthI].setSize
            (
                nVelocitySamples_
            );
            ringCellI_[azimuthI] = -1;

            ringProcI_[azimuthI].setSize
            (
                nVelocitySamples_
            );
            ringProcI_[azimuthI] = -1;

            velocitiesRing_[azimuthI].setSize
            (
                nVelocitySamples_
            );
            velocitiesRing_[azimuthI] = vector::zero;

            previousRingLocation_[azimuthI].setSize
            (
                nVelocitySamples_
            );
            previousRingLocation_[azimuthI] = vector::zero;
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
    scalar projectionRadius
        = (epsilon_[azimuthIndex_]*Foam::sqrt(Foam::log(1.0/0.001)));

    // Apply force to the cells within the element's sphere of influence
    scalar sphereRadius = chordLength_ + projectionRadius;
    scalar sphereRadiusSqr = sphereRadius*sphereRadius;

    const vectorField& C = mesh_.C();

    label nCells = 0;

    forAll(C, cellI)
    {
        scalar dis = magSqr(C[cellI] - position_[azimuthIndex_]);
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
    label azimuthIndex,
    bool clearBuffer
)
{
    azimuthIndex_ = azimuthIndex;
    if (clearBuffer && azimuthIndex == 0)
    {
        stringBuffer_.str("");
        stringBuffer_.clear();
    }
}

label Foam::fv::actuatorLineElement::findNearbyCell
(
    const point &location,
    const label previousCell
)
{
    label localCell = -1;

    // Check if previousCell still is valid
    if (mesh_.pointInCell(location, previousCell))
    {
        localCell = previousCell;
    }
    else
    {
        // If not previous cell, check neighboring cells
        const labelList& nbrs =
            mesh_.cellCells()[previousCell];

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
    }
    return localCell;
}

void Foam::fv::actuatorLineElement::findCells(bool includeRing)
{
    // Find local flow velocity by interpolating to element location
    vector inflowVelocityPoint = position_[azimuthIndex_];

    if (centerProcI_[azimuthIndex_] == Pstream::myProcNo()
                && centerCellI_[azimuthIndex_] >= 0
                && centerCellI_[azimuthIndex_] < mesh_.nCells())
    {
        if
        (
            magSqr
            (
                position_[azimuthIndex_] - previousLocation_[azimuthIndex_]
            )
            > SMALL
        )
        {
            // if we cannot reuse the previous value directly
            centerCellI_[azimuthIndex_] =
                findNearbyCell
                (
                    position_[azimuthIndex_],
                    centerCellI_[azimuthIndex_]
                );
            if (centerCellI_[azimuthIndex_] < 0) // if lookup failed
            {
                centerProcI_[azimuthIndex_] = -1;
            }
        }
    }
    else
    {
        centerCellI_[azimuthIndex_] = -1;
        centerProcI_[azimuthIndex_] = -1;
    }
    previousLocation_[azimuthIndex_] = position_[azimuthIndex_];

    // If the flow is sampled by using a circle around position_
    if (includeRing && velocitySampleRadius_ > 0.0)
    {
        // Circle radius should be normalized with epsilon
        // Use old value of epsilon for speed

        scalar sampleRadius = epsilon_[azimuthIndex_]*velocitySampleRadius_;

        // Unit vector in chordwise direction
        vector chordNormal = chordDirection_[azimuthIndex_] /
                                mag(chordDirection_[azimuthIndex_]);

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


            if (ringProcI_[azimuthIndex_][point] == Pstream::myProcNo()
                && ringCellI_[azimuthIndex_][point] >= 0
                && ringCellI_[azimuthIndex_][point] < mesh_.nCells())
            {
                if
                (
                    magSqr
                    (
                        samplePoint -
                            previousRingLocation_[azimuthIndex_][point]
                    )
                    > SMALL
                )
                {
                    ringCellI_[azimuthIndex_][point] =
                        findNearbyCell
                        (
                            samplePoint,
                            ringCellI_[azimuthIndex_][point]
                        );

                    // if lookup failed
                    if (ringCellI_[azimuthIndex_][point] < 0)
                    {
                        ringProcI_[azimuthIndex_][point] = -1;
                    }
                }
            }
            else
            {
                ringCellI_[azimuthIndex_][point] = -1;
                ringProcI_[azimuthIndex_][point] = -1;
            }
            previousRingLocation_[azimuthIndex_][point] = samplePoint;
        }
    }
}

void Foam::fv::actuatorLineElement::calculateRingVelocity()
{
    vector localVelocitySum = vector::zero;
    
    // If the flow is sampled by using a circle around position_, then
    // overwrite the inflow velocity with the mean value over all circle points
    if (velocitySampleRadius_ > 0.0)
    {
        // Calculate mean value over all circle points
        for (label point = 0; point < nVelocitySamples_; point++)
        {
            localVelocitySum += velocitiesRing_[azimuthIndex_][point];
        }
        // Set inflow Velocity as the mean value
        inflowVelocity_ = 1.0 / nVelocitySamples_ * localVelocitySum;
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
        outputFile_<< "time,root_dist,x,y,z,rel_vel_mag,Ux,Re,alpha_deg,"
                   << "alpha_geom_deg,phi,cl,cd,fx,fy,fz,end_effect_factor,"
                   << "c_ref_t,c_ref_n,f_ref_t,f_ref_n" << std::endl;
    }
}


void Foam::fv::actuatorLineElement::writePerf()
{
    scalar time = mesh_.time().value();

    // write time,root_dist,x,y,z,rel_vel_mag,Ux,Re,alpha_deg,alpha_geom_deg,
    // phi,cl,cd,fx,fy,fz,end_effect_factor,c_ref_t,c_ref_n,f_ref_t,f_ref_n
    stringBuffer_<< time << ","
            << rootDistance_ << ","
            << position_[azimuthIndex_].x() << ","
            << position_[azimuthIndex_].y() << ","
            << position_[azimuthIndex_].z() << ","
            << mag(relativeVelocity_[azimuthIndex_]) << ","
            << relativeVelocity_[azimuthIndex_].x() << ","
            << Re_ << ","
            << angleOfAttack_ << ","
            << angleOfAttackGeom_ << ","
            << inflowRefAngle() << ","
            << liftCoefficient_[azimuthIndex_] << ","
            << dragCoefficient_[azimuthIndex_] << ","
            << forceVector_[azimuthIndex_].x()*localRho_[azimuthIndex_] << ","
            << forceVector_[azimuthIndex_].y()*localRho_[azimuthIndex_] << ","
            << forceVector_[azimuthIndex_].z()*localRho_[azimuthIndex_] << ","
            << endEffectFactor_ << "," << tangentialRefCoefficient() << ","
            << normalRefCoefficient() << ","
            << tangentialRefForce()*localRho_[azimuthIndex_] << ","
            << normalRefForce()*localRho_[azimuthIndex_] << std::endl;

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
    chordDirection_(1, vector::zero),
    spanDirection_(1, vector::zero),
    planformNormal_(vector::zero),
    chordRefDirection_(1, vector::zero),
    position_(1, vector::zero),
    velocity_(1, vector::zero),
    freeStreamVelocity_(vector::zero),
    forceVector_(1, vector::zero),
    localRho_(1, 1.0),
    localMu_(1, -1.0),
    inflowVelocity_(1, vector::zero),
    epsilon_(1, 0.0),
    previousLocation_(1, point::zero),
    previousRingLocation_(0),
    centerCellI_(1, -1),
    centerProcI_(1, -1),
    ringCellI_(0),
    ringProcI_(0),
    relativeVelocity_(1, vector::zero),
    relativeVelocityGeom_(vector::zero),
    angleOfAttack_(0.0),
    angleOfAttackGeom_(0.0),
    liftCoefficient_(1, 0.0),
    dragCoefficient_(1, 0.0),
    momentCoefficient_(1, 0.0),
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
    angleOfAttackDiff_(0.0),
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
    return position_[azimuthIndex_];
}


const Foam::vector& Foam::fv::actuatorLineElement::velocity()
{
    return velocity_[azimuthIndex_];
}


const Foam::vector& Foam::fv::actuatorLineElement::relativeVelocity()
{
    return relativeVelocity_[azimuthIndex_];
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
    return liftCoefficient_[azimuthIndex_];
}


const Foam::scalar& Foam::fv::actuatorLineElement::dragCoefficient()
{
    return dragCoefficient_[azimuthIndex_];
}


const Foam::scalar& Foam::fv::actuatorLineElement::momentCoefficient()
{
    return momentCoefficient_[azimuthIndex_];
}


Foam::scalar Foam::fv::actuatorLineElement::tangentialRefCoefficient()
{
    return profileData_.convertToCRT
    (
        liftCoefficient_[azimuthIndex_],
        dragCoefficient_[azimuthIndex_],
        inflowRefAngle()
    );
}


Foam::scalar Foam::fv::actuatorLineElement::tangentialRefForce()
{
    return 0.5 * chordLength_ * tangentialRefCoefficient()
        * magSqr(relativeVelocity_[azimuthIndex_]);
}


Foam::scalar Foam::fv::actuatorLineElement::normalRefCoefficient()
{
    return profileData_.convertToCRN
    (
        liftCoefficient_[azimuthIndex_],
        dragCoefficient_[azimuthIndex_],
        inflowRefAngle()
    );
}


Foam::scalar Foam::fv::actuatorLineElement::normalRefForce()
{
    return 0.5 * chordLength_ * normalRefCoefficient()
        * magSqr(relativeVelocity_[azimuthIndex_]);
}


Foam::scalar Foam::fv::actuatorLineElement::inflowRefAngle()
{
    // Calculate inflow velocity angle in degrees (AFTAL Phi)
    scalar arg =
        (-relativeVelocity_[azimuthIndex_] & chordRefDirection_[azimuthIndex_])
        / (mag(relativeVelocity_[azimuthIndex_])
        * mag(chordRefDirection_[azimuthIndex_]));
    scalar inflowVelAngleRad =
        acos(sign(arg)*min(Foam::scalar(1.0), mag(arg)));
    return radToDeg(inflowVelAngleRad - angleOfAttackDiff_);
}


const Foam::scalar& Foam::fv::actuatorLineElement::rootDistance()
{
    return rootDistance_;
}


void Foam::fv::actuatorLineElement::calculateForce()
{
    scalar pi = Foam::constant::mathematical::pi;
    
    // compressible case, localMu is -1 for the incompressible case
    if (localMu_[azimuthIndex_] > 0)
    {
        nu_ = localMu_[azimuthIndex_]/localRho_[azimuthIndex_];
    }

    // Calculate vector normal to chord--span plane
    planformNormal_ =
        -chordDirection_[azimuthIndex_] ^ spanDirection_[azimuthIndex_];
    planformNormal_ /= mag(planformNormal_);

    if (debug)
    {
        Info<< "Calculating force contribution from actuatorLineElement "
            << name_ << endl;
        Info<< "    position: " << position_[azimuthIndex_] << endl;
        Info<< "    chordDirection: " << chordDirection_[azimuthIndex_] << endl;
        Info<< "    spanDirection: " << spanDirection_[azimuthIndex_] << endl;
        Info<< "    elementVelocity: " << velocity_[azimuthIndex_] << endl;
        Info<< "    planformNormal: " << planformNormal_ << endl;
    }

    // Find local flow velocity by interpolating to element location
    calculateRingVelocity();
    

    // Subtract spanwise component of inflow velocity
    vector spanwiseVelocity =
        spanDirection_[azimuthIndex_]
            * (inflowVelocity_[azimuthIndex_] & spanDirection_[azimuthIndex_])
            / magSqr(spanDirection_[azimuthIndex_]);
    inflowVelocity_[azimuthIndex_] -= spanwiseVelocity;

    // Calculate relative velocity and Reynolds number
    relativeVelocity_[azimuthIndex_] =
        inflowVelocity_[azimuthIndex_] - velocity_[azimuthIndex_];
    Re_ = mag(relativeVelocity_[azimuthIndex_])*chordLength_/nu_;

    // Calculate angle of attack (radians)
    scalar arg =
        (planformNormal_ & relativeVelocity_[azimuthIndex_])
        / (mag(planformNormal_) * mag(relativeVelocity_[azimuthIndex_]));
    scalar angleOfAttackRad = asin(sign(arg)*min(Foam::scalar(1.0), mag(arg)));
    scalar angleOfAttackUncorrected = radToDeg(angleOfAttackRad);
    relativeVelocityGeom_ = freeStreamVelocity_ - velocity_[azimuthIndex_];
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

    // Update Reynolds number of profile data
    profileData_.updateRe(Re_);

    // Calculate angle of attack in degrees
    angleOfAttack_ = radToDeg(angleOfAttackRad);

    // If we correct the angle of attack, reduce it with using endEffectFactor_
    // if endEffectFactor_ = 0, we should have the angle of attack that
    // gives zero lift coefficient
    angleOfAttackDiff_ = 0.0;
    if (endEffectsCorrectAngleOfAttack_ || endEffectsCorrectAngleOfRelWind_)
    {
        scalar zeroLiftAngleOfAttack = profileData_.zeroLiftAngleOfAttack();
        scalar angleOfAttackNew =
            (angleOfAttack_ - zeroLiftAngleOfAttack)*endEffectFactor_
            + zeroLiftAngleOfAttack;
        angleOfAttackDiff_ = degToRad(angleOfAttack_ - angleOfAttackNew);
        if (endEffectsCorrectAngleOfAttack_)
        {
            angleOfAttack_ = angleOfAttackNew;
        }
    }

    // Lookup lift and drag coefficients
    lookupCoefficients();

    if (debug)
    {
        Info<< "    inflowVelocity: "
            << inflowVelocity_[azimuthIndex_] << endl;
        Info<< "    relativeVelocity: "
            << relativeVelocity_[azimuthIndex_] << endl;
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
            mag(relativeVelocity_[azimuthIndex_]),
            angleOfAttack_,
            liftCoefficient_[azimuthIndex_],
            dragCoefficient_[azimuthIndex_],
            momentCoefficient_[azimuthIndex_]
        );
    }

    // Correct for added mass effects
    if (addedMassActive_)
    {
        addedMass_.correct
        (
            liftCoefficient_[azimuthIndex_],
            dragCoefficient_[azimuthIndex_],
            momentCoefficient_[azimuthIndex_],
            degToRad(angleOfAttack_),
            mag
            (
                chordDirection_[azimuthIndex_]
                    & relativeVelocity_[azimuthIndex_]
            ),
            mag(planformNormal_ & relativeVelocity_[azimuthIndex_])
        );
    }

    // Apply end effect correction factor to lift coefficient
    if (endEffectsCorrectAngleOfAttack_ == false)
    {
        liftCoefficient_[azimuthIndex_] *= endEffectFactor_;
    }

    // Calculate force per unit density
    scalar area = chordLength_ * spanLength_;
    scalar magSqrU = magSqr(relativeVelocity_[azimuthIndex_]);
    scalar lift = 0.5*area*liftCoefficient_[azimuthIndex_]*magSqrU;
    scalar drag = 0.5*area*dragCoefficient_[azimuthIndex_]*magSqrU;
    vector liftDirection =
        relativeVelocity_[azimuthIndex_] ^ spanDirection_[azimuthIndex_];
    liftDirection /= mag(liftDirection);
    vector dragDirection =
        relativeVelocity_[azimuthIndex_]/mag(relativeVelocity_[azimuthIndex_]);

    if (endEffectsCorrectAngleOfRelWind_)
    {
        rotateVector
        (
            liftDirection,
            vector::zero,
            spanDirection_[azimuthIndex_],
            angleOfAttackDiff_
        );
        rotateVector
        (
            dragDirection,
            vector::zero,
            spanDirection_[azimuthIndex_],
            angleOfAttackDiff_
        );
    }
    forceVector_[azimuthIndex_] = lift*liftDirection + drag*dragDirection;

    if (debug)
    {
        Info<< "    liftDirection: " << liftDirection << endl;
        Info<< "    dragDirection: " << dragDirection << endl;
        Info<< "    force (per unit density): "
            << forceVector_[azimuthIndex_] << endl;
    }
    calcProjectionEpsilon();

    // Write performance to file
    if (Pstream::master() && (writePerf_ || writePerfEnd_))
    {
        writePerf();
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
        Info<< "Initial position: " << position_[azimuthIndex_] << endl;
        Info<< "Initial chordDirection: "
            << chordDirection_[azimuthIndex_] << endl;
        Info<< "Initial spanDirection: "
            << spanDirection_[azimuthIndex_] << endl;
        Info<< "Initial velocity: " << velocity_[azimuthIndex_] << endl;
    }

    // Rotation matrices make a rotation about the origin, so need to subtract
    // rotation point off the point to be rotated.
    vector point = position_[azimuthIndex_];
    point -= rotationPoint;

    // Perform the rotation.
    point = RM & point;

    // Return the rotated point to its new location relative to the rotation
    // point
    point += rotationPoint;

    label nextIndex = (azimuthIndex_ + 1)% position_.size();
    // Set the position of the element
    position_[nextIndex] = point;

    // Rotate the span and chord vectors of the element
    chordDirection_[nextIndex] = RM & chordDirection_[azimuthIndex_];
    spanDirection_[nextIndex] = RM & spanDirection_[azimuthIndex_];

    // Rotate the element's velocity vector if specified
    if (rotateVelocity)
    {
        velocity_[nextIndex] = RM & velocity_[azimuthIndex_];
        chordRefDirection_[nextIndex] = RM & chordRefDirection_[azimuthIndex_];
    }

    if (debug)
    {
        Info<< "Final position: " << position_[nextIndex] << endl;
        Info<< "Final chordDirection: " << chordDirection_[nextIndex] << endl;
        Info<< "Final chordRefDirection: "
            << chordRefDirection_[nextIndex] << endl;
        Info<< "Final spanDirection: " << spanDirection_[nextIndex] << endl;
        Info<< "Final velocity: " << velocity_[nextIndex] << endl << endl;
    }
}


void Foam::fv::actuatorLineElement::pitch
(
    scalar radians,
    scalar chordFraction
)
{
    vector rotationPoint = position_[azimuthIndex_];
    rotationPoint +=
        chordDirection_[azimuthIndex_]*(chordMount_ - chordFraction);
    rotate(rotationPoint, spanDirection_[azimuthIndex_], radians, false);
}

label Foam::fv::actuatorLineElement::locationCount() const
{
    label n = position_.size();  // Central point

    if (velocitySampleRadius_ > 0.0)
    {
        if (previousRingLocation_.empty() == false)
        {
            n += previousRingLocation_.size() * previousRingLocation_[0].size();
        }
    }
    return n;
}

label Foam::fv::actuatorLineElement::epsilonCount() const
{
    return epsilon_.size();
}

const scalar &Foam::fv::actuatorLineElement::chordFactor() const
{
    return chordFactor_;
}

const scalar &Foam::fv::actuatorLineElement::dragFactor() const
{
    return dragFactor_;
}

const scalar &Foam::fv::actuatorLineElement::meshFactor() const
{
    return meshFactor_;
}

void Foam::fv::actuatorLineElement::collectLocationData
(
    List<label> &globalCellI,
    List<label> &globalProcI,
    List<point> &globalLocations,
    label &index,
    bool includeRing
) const
{
    const label nCenter = position_.size();

    // --- Center ---
    SubList<point>(globalLocations, nCenter, index) = position_;
    SubList<label>(globalCellI,     nCenter, index) = centerCellI_;
    SubList<label>(globalProcI,     nCenter, index) = centerProcI_;
    index += nCenter;

    // --- Rings ---
    if (velocitySampleRadius_ > 0.0 && includeRing)
    {
        forAll(previousRingLocation_, ringI)
        {
            const List<point>&  ringPos  = previousRingLocation_[ringI];
            const List<label>&  ringCell = ringCellI_[ringI];
            const List<label>&  ringProc = ringProcI_[ringI];

            const label nRing = ringPos.size();

            SubList<point>(globalLocations, nRing, index) = ringPos;
            SubList<label>(globalCellI,     nRing, index) = ringCell;
            SubList<label>(globalProcI,     nRing, index) = ringProc;

            index += nRing;
        }
    }
}

void Foam::fv::actuatorLineElement::distributeVelocityData
(
    const List<label> &globalCellI,
    const List<label> &globalProcI,
    const List<vector> &globalVelocities,
    label &index
)
{
    const label nCenter = position_.size();

    // --- Center ---
    centerCellI_ = SubList<label>(globalCellI, nCenter, index);
    centerProcI_ = SubList<label>(globalProcI, nCenter, index);
    SubList<vector>(inflowVelocity_, nCenter, 0) =
        SubList<vector>(globalVelocities, nCenter, index);

    index += nCenter;

    // --- Rings ---
    if (velocitySampleRadius_ > 0.0)
    {
        forAll(previousRingLocation_, ringI)
        {
            List<label>&  ringCell = ringCellI_[ringI];
            List<label>&  ringProc = ringProcI_[ringI];
            List<vector>& ringVel  = velocitiesRing_[ringI];

            const label nRing = ringCell.size();

            ringCell = SubList<label>(globalCellI,  nRing, index);
            ringProc = SubList<label>(globalProcI,  nRing, index);
            ringVel  = SubList<vector>(globalVelocities, nRing, index);

            index += nRing;
        }
    }
}

void Foam::fv::actuatorLineElement::distributeCenterCellI
(
    const List<label>& globalCellI,
    const List<label>& globalProcI,
    label& index
)
{
    const label nCenter = centerCellI_.size();
    centerCellI_ = SubList<label>(globalCellI, nCenter, index);
    centerProcI_ = SubList<label>(globalProcI, nCenter, index);
    index += nCenter;
}

void Foam::fv::actuatorLineElement::collectEpsilonData
(
    List<scalar> &globalEpsilon,
    label &index
) const
{
    const label nCenter = epsilon_.size();
    SubList<scalar>(globalEpsilon, nCenter, index) = epsilon_;
    index += nCenter;
}

void Foam::fv::actuatorLineElement::calcInfluenceEpsilon
(
    scalar dragCoefficient
)
{
    // Calculate projection width
    dragCoefficient_ = dragCoefficient;
    calcProjectionEpsilon();
}

void Foam::fv::actuatorLineElement::distributeEpsilonData
(
    List<scalar> &globalEpsilon,
    label &index
)
{
    const label nCenter = epsilon_.size();
    epsilon_ = SubList<scalar>(globalEpsilon, nCenter, index);
    index += nCenter;
}

void Foam::fv::actuatorLineElement::distributeRhoData
(
    List<scalar> &globalRho,
    label &index
)
{
    const label nCenter = localRho_.size();
    localRho_ = SubList<scalar>(globalRho, nCenter, index);
    index += nCenter;
}

void Foam::fv::actuatorLineElement::distributeMuData
(
    List<scalar> &globalMu,
    label &index
)
{
    const label nCenter = localMu_.size();
    localMu_ = SubList<scalar>(globalMu, nCenter, index);
    index += nCenter;
}

void Foam::fv::actuatorLineElement::translate(vector translationVector)
{
    position_[azimuthIndex_] += translationVector;
}


void Foam::fv::actuatorLineElement::setVelocity(vector velocity)
{
    if (debug)
    {
        Info<< "Changing velocity of " << name_ << " from "
            << velocity_[azimuthIndex_] << " to " << velocity << endl << endl;
    }
    velocity_[azimuthIndex_] = velocity;
}


void Foam::fv::actuatorLineElement::setSpeed(scalar speed)
{
    if (mag(velocity_[azimuthIndex_]) > 0)
    {
        velocity_[azimuthIndex_] /= mag(velocity_[azimuthIndex_]);
        velocity_[azimuthIndex_] *= speed;
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
        Info<< "    Initial velocity: " << velocity_[azimuthIndex_] << endl;
    }

    // First find radius from axis to element position -- formula from
    // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    vector point2 = point + axis;
    scalar radius = mag((position_[azimuthIndex_] - point) ^
                    (position_[azimuthIndex_] - point2))
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
        velocityLE_ = velocity_[azimuthIndex_]*radiusLE/radius;
        rotateVector
        (
            velocityLE_,
            vector::zero,
            spanDirection_[azimuthIndex_],
            angleLE
        );

        // Set velocity at trailing edge
        scalar radiusTE = sqrt(magSqr(0.75*chordLength_) + magSqr(radius));
        angleTE = atan2(-0.75*chordLength_, radius);
        velocityTE_ = velocity_[azimuthIndex_]*radiusTE/radius;
        rotateVector
        (
            velocityTE_,
            vector::zero,
            spanDirection_[azimuthIndex_],
            angleTE
        );
    }

    // Also set omega for flow curvature correction
    setOmega(omega);

    if (debug)
    {
        Info<< "    Radius: " << radius << endl;
        Info<< "    Final velocity: " << velocity_[azimuthIndex_] << endl;
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
    velocity_[azimuthIndex_] *= scale;
}


const Foam::vector& Foam::fv::actuatorLineElement::force()
{
    return forceVector_[azimuthIndex_];
}


Foam::vector Foam::fv::actuatorLineElement::moment(vector point)
{
    // Calculate radius vector
    vector radius = position_[azimuthIndex_] - point;
    vector moment = radius ^ forceVector_[azimuthIndex_];
    vector pitchingMoment = 0.5*chordLength_*chordLength_*spanLength_
                          * momentCoefficient_[azimuthIndex_]
                          * magSqr(relativeVelocity_[azimuthIndex_])
                          * spanDirection_[azimuthIndex_];
    return moment + pitchingMoment;
}


void Foam::fv::actuatorLineElement::addForce
(
    volVectorField& forceField,
    scalar scale,
    bool compressible
)
{
    applyForceField(forceField, scale);

    //if (compressible)
    //{
        // Multiply force vector by local density
    //    multiplyForceRho();
    //}
}

void Foam::fv::actuatorLineElement::addTurbulence
(
    fvMatrix<scalar>& eqn,
    word fieldName
)
{

    // Calculate projection radius
    //scalar epsilon = calcProjectionEpsilon();
    scalar projectionRadius =
        (epsilon_[azimuthIndex_]*Foam::sqrt(Foam::log(1.0/0.001)));

    // Calculate TKE injection rate
    scalar k = 0.1*mag(dragCoefficient_[azimuthIndex_]);

    // Add turbulence to the cells within the element's sphere of influence
    scalar sphereRadius = chordLength_ + projectionRadius;
    scalar sphereRadiusSqr = sphereRadius*sphereRadius;
    scalar invepsilonSqr =
        1.0/(epsilon_[azimuthIndex_]*epsilon_[azimuthIndex_]);
    scalar internalFactor = 1.0/(Foam::pow(epsilon_[azimuthIndex_], 3)
                          * Foam::pow(Foam::constant::mathematical::pi, 1.5));
                          
    const vectorField& C = mesh_.C();
    scalarField& src = eqn.source();
    const scalarField& V = mesh_.V();

    if (influenceCells_.empty())
    {
        forAll(mesh_.cells(), cellI)
        {
            scalar dis = magSqr(C[cellI] - position_[azimuthIndex_]);
            if (dis <= sphereRadiusSqr)
            {
                scalar factor = Foam::exp(-dis*invepsilonSqr)*internalFactor;
                if (fieldName == "k")
                {
                    //turbulence[cellI] += factor*k;
                    src[cellI] += factor*k * V[cellI];
                }
                else if (fieldName == "epsilon")
                {
                    //turbulence[cellI] += factor*Foam::pow(k, 1.5)
                    //              * 0.09/(chordLength_/10.0);
                    src[cellI] += factor*Foam::pow(k, 1.5)
                                  * 0.09/(chordLength_/10.0)
                                  * V[cellI];
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
            scalar dis = magSqr(C[cellI] - position_[azimuthIndex_]);
            if (dis <= sphereRadiusSqr)
            {
                scalar factor = Foam::exp(-dis*invepsilonSqr)*internalFactor;
                if (fieldName == "k")
                {
                    //turbulence[cellI] += factor*k;
                    src[cellI] += factor*k * V[cellI];
                }
                else if (fieldName == "epsilon")
                {
                    //turbulence[cellI] += factor*Foam::pow(k, 1.5)
                    //              * 0.09/(chordLength_/10.0);
                    src[cellI] += factor*Foam::pow(k, 1.5)
                                  * 0.09/(chordLength_/10.0)
                                  * V[cellI];
                }
            }
        }
    }

    //eqn += turbulence;
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


vector& Foam::fv::actuatorLineElement::inflowVelocity()
{
    return inflowVelocity_[azimuthIndex_];
}


const vector& Foam::fv::actuatorLineElement::inflowVelocity() const
{
    return inflowVelocity_[azimuthIndex_];
}


const scalar& Foam::fv::actuatorLineElement::epsilon() const
{
    return epsilon_[azimuthIndex_];
}
// ************************************************************************* //
