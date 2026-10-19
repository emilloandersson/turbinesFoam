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

#include "actuatorLineSource.H"
#include "interpolateUtils.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "vector.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "simpleMatrix.H"
#include "SVD.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorLineSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        actuatorLineSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::fv::actuatorLineSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // Look up information in dictionary
        coeffs_.lookup("elementProfiles") >> elementProfiles_;
        profileData_ = coeffs_.subDict("profileData");
        coeffs_.lookup("elementGeometry") >> elementGeometry_;
        coeffs_.lookup("nElements") >> nElements_;
        coeffs_.lookup("freeStreamVelocity") >> freeStreamVelocity_;
        freeStreamDirection_ = freeStreamVelocity_/mag(freeStreamVelocity_);
        endEffectsActive_ = coeffs_.lookupOrDefault("endEffects", false);
        maxNElementsEndEffects_ =
            coeffs_.lookupOrDefault<label>("maxNElementsEndEffects", 240);
        maxConditionNumberEndEffects_ = coeffs_.lookupOrDefault<scalar>
        (
            "maxConditionNumberEndEffects",
            1e10
        );
        maxResidualEndEffects_ =
            coeffs_.lookupOrDefault<scalar>("maxResidualEndEffects", 1e-8);
        normalizeEndEffects_ =
            coeffs_.lookupOrDefault<bool>("normalizeEndEffects", true);
        ignoreLowValuesEndEffects_ =
            coeffs_.lookupOrDefault<scalar>("ignoreLowValuesEndEffects", 0.2);
        ignoreHighValuesEndEffects_ =
            coeffs_.lookupOrDefault<scalar>("ignoreHighValuesEndEffects", 0.2);

        // Read harmonic pitching parameters if present
        dictionary pitchDict = coeffs_.subOrEmptyDict("harmonicPitching");
        harmonicPitchingActive_ = pitchDict.lookupOrDefault("active", false);
        reducedFreq_ = pitchDict.lookupOrDefault("reducedFreq", 0.0);
        pitchAmplitude_ = pitchDict.lookupOrDefault("amplitude", 0.0);

        if (debug)
        {
            Info<< "Debugging for actuatorLineSource on" << endl;
            printCoeffs();
        }

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::actuatorLineSource::createOutputFile()
{
    fileName dir;

    if (Pstream::parRun())
    {
        dir = mesh_.time().path()/"../postProcessing/actuatorLines"
            / mesh_.time().timeName();
    }
    else
    {
        dir = mesh_.time().path()/"postProcessing/actuatorLines"
            / mesh_.time().timeName();
    }

    if (not isDir(dir))
    {
        mkDir(dir);
    }

    outputFile_.open(dir/name_ + ".csv", std::ios::out);

    if (outputFile_.is_open())
    {
        outputFile_ << "time,x,y,z,rel_vel_mag,alpha_deg,alpha_geom_deg,cl,"
                    << "cd,cm" << std::endl;
    }
}


void Foam::fv::actuatorLineSource::createElements()
{
    elements_.setSize(nElements_);

    label nGeometryPoints = elementGeometry_.size();
    label nGeometrySegments = nGeometryPoints - 1;
    label nElementsPerSegment = nElements_/nGeometrySegments;
    if (nElements_ % nGeometrySegments)
    {
        // Need to have integer number of elements per geometry segment
        FatalErrorIn("void actuatorLineSource::createElements()")
            << "Number of actuator line elements must be multiple of the "
            << "number of actuator line geometry segments"
            << abort(FatalError);
    }
    List<vector> points(nGeometryPoints);
    List<vector> spanDirs(nGeometryPoints);
    List<scalar> chordLengths(nGeometryPoints);
    List<scalar> spanLengths(nGeometrySegments);
    List<vector> chordRefDirs(nGeometryPoints);
    List<scalar> pitches(nGeometryPoints);
    List<scalar> chordMounts(nGeometryPoints);
    List<scalar> coneAngles(nGeometryPoints);
    totalLength_ = 0.0;
    chordLength_ = 0.0;

    forAll(points, i)
    {
        // Extract geometry point
        scalar x = elementGeometry_[i][0][0];
        scalar y = elementGeometry_[i][0][1];
        scalar z = elementGeometry_[i][0][2];
        points[i] = vector(x, y, z);
        if (i > 0)
        {
            spanLengths[i - 1] = mag(points[i] - points[i-1]);
            totalLength_ += spanLengths[i - 1];
        }
        // Read span direction
        x = elementGeometry_[i][1][0];
        y = elementGeometry_[i][1][1];
        z = elementGeometry_[i][1][2];
        spanDirs[i] = vector(x, y, z);
        // Read chord length
        chordLengths[i] = elementGeometry_[i][2][0];
        chordLength_ += chordLengths[i];
        // Read chord ref dir
        x = elementGeometry_[i][3][0];
        y = elementGeometry_[i][3][1];
        z = elementGeometry_[i][3][2];
        chordRefDirs[i] = vector(x, y, z);
        // Read chord mount
        chordMounts[i] = elementGeometry_[i][4][0];
        // Read pitch
        pitches[i] = elementGeometry_[i][5][0];
        // coneAngle of 90 degrees means a horizontal axis turbine here
        coneAngles[i] =
            elementGeometry_[i].size() > 6 ? elementGeometry_[i][6][0] :
            Foam::constant::mathematical::pi/2;
    }

    // Store blade root and tip locations for distance calculations
    vector rootLocation = points[0];
    vector tipLocation = points[nGeometryPoints - 1];

    // Compute average chord length
    chordLength_ /= nGeometryPoints;

    // Compute aspect ratio
    aspectRatio_ = totalLength_/chordLength_;

    // Lookup initial element velocities if present
    List<vector> initialVelocities(nGeometryPoints, vector::zero);
    coeffs_.readIfPresent("initialVelocities", initialVelocities);

    if (debug)
    {
        Info<< "Total length: " << totalLength_ << endl;
        Info<< "Elements per geometry segment: " << nElementsPerSegment
            << endl;
        Info<< "Points:" << endl << points << endl;
        Info<< "Span directions:" << endl << spanDirs << endl;
        Info<< "Span lengths: " << endl << spanLengths << endl;
        Info<< "Chord lengths:" << endl << chordLengths << endl;
        Info<< "Pitches:" << endl << pitches << endl;
        Info<< "Root location: " << rootLocation << endl;
        Info<< "Tip location: " << tipLocation << endl;
    }

    forAll(elements_, i)
    {
        std::stringstream ss;
        ss << i;
        string str = ss.str();
        const word name = name_ + ".element" + str;

        // Actuator point geometry to be calculated from elementGeometry
        label geometrySegmentIndex = i/nElementsPerSegment;
        label pointIndex = i % nElementsPerSegment;
        label elementProfileIndex = i*elementProfiles_.size()/nElements_;
        word profileName = elementProfiles_[elementProfileIndex];
        vector position;
        scalar chordLength;
        vector chordDirection;
        vector chordRefDirection;
        scalar spanLength = spanLengths[geometrySegmentIndex];
        spanLength /= nElementsPerSegment;
        vector spanDirection;
        scalar pitch;
        scalar chordMount;
        scalar cone;
        vector initialVelocity;

        // Linearly interpolate position
        vector point1 = points[geometrySegmentIndex];
        vector point2 = points[geometrySegmentIndex + 1];
        vector segment = point2 - point1;
        position = point1
                 + segment/nElementsPerSegment*pointIndex
                 + segment/nElementsPerSegment/2;

        // Linearly interpolate chordLength
        scalar chordLength1 = chordLengths[geometrySegmentIndex];
        scalar chordLength2 = chordLengths[geometrySegmentIndex + 1];
        scalar deltaChordTotal = chordLength2 - chordLength1;
        chordLength = chordLength1
                    + deltaChordTotal/nElementsPerSegment*pointIndex
                    + deltaChordTotal/nElementsPerSegment/2;

        // Linearly interpolate spanDirection
        vector spanDir1 = spanDirs[geometrySegmentIndex];
        vector spanDir2 = spanDirs[geometrySegmentIndex + 1];
        vector deltaSpanTotal = spanDir2 - spanDir1;
        spanDirection = spanDir1
                      + deltaSpanTotal/nElementsPerSegment*pointIndex
                      + deltaSpanTotal/nElementsPerSegment/2;

        // Linearly interpolate section pitch
        scalar pitch1 = pitches[geometrySegmentIndex];
        scalar pitch2 = pitches[geometrySegmentIndex + 1];
        scalar deltaPitchTotal = pitch2 - pitch1;
        pitch = pitch1
              + deltaPitchTotal/nElementsPerSegment*pointIndex
              + deltaPitchTotal/nElementsPerSegment/2;

        // Linearly interpolate chord mount
        scalar cm1 = chordMounts[geometrySegmentIndex];
        scalar cm2 = chordMounts[geometrySegmentIndex + 1];
        scalar deltaCmTotal = cm2 - cm1;
        chordMount = cm1 + deltaCmTotal/nElementsPerSegment*pointIndex
                   + deltaCmTotal/nElementsPerSegment/2;

        // Linearly interpolate element velocity
        vector vel1 = initialVelocities[geometrySegmentIndex];
        vector vel2 = initialVelocities[geometrySegmentIndex + 1];
        vector deltaVelTotal = vel2 - vel1;
        initialVelocity = vel1
                        + deltaVelTotal/nElementsPerSegment*pointIndex
                        + deltaVelTotal/nElementsPerSegment/2;

        // Linearly interpolate chordDirection
        vector chordDir1 = chordRefDirs[geometrySegmentIndex];
        vector chordDir2 = chordRefDirs[geometrySegmentIndex + 1];
        vector deltaChordDirTotal = chordDir2 - chordDir1;
        chordDirection = chordDir1
                       + deltaChordDirTotal/nElementsPerSegment*pointIndex
                       + deltaChordDirTotal/nElementsPerSegment/2;
                       
        // Linearly interpolate cone angle
        scalar cone1 = coneAngles[geometrySegmentIndex];
        scalar cone2 = coneAngles[geometrySegmentIndex + 1];
        scalar deltaConeTotal = cone2 - cone1;
        cone = cone1 + deltaConeTotal/nElementsPerSegment*pointIndex
             + deltaConeTotal/nElementsPerSegment/2;

        // Chord reference direction (before pitching)
        chordRefDirection = chordDirection;
        
        // Calculate nondimensional root distance
        scalar rootDistance = mag(position - rootLocation)/totalLength_;

        // Create a dictionary for this actuatorLineElement
        dictionary dict;
        dict.add("position", position);
        dictionary profileDataDict = profileData_.subDict(profileName);
        dict.add("profileData", profileDataDict);
        dict.add("profileName", profileName);
        dict.add("chordLength", chordLength);
        dict.add("chordDirection", chordDirection);
        dict.add("chordRefDirection", chordRefDirection);
        dict.add("spanLength", spanLength);
        dict.add("spanDirection", spanDirection);
        dict.add("freeStreamVelocity", freeStreamVelocity_);
        dict.add("chordMount", chordMount);
        dict.add("rootDistance", rootDistance);
        dict.add("cone", cone);
        dict.add("addedMass", coeffs_.lookupOrDefault("addedMass", false));
        dict.add
        (
            "endEffectsCorrectAoA",
            coeffs_.lookupOrDefault("endEffectsCorrectAoA", false)
        );
        dict.add
        (
            "endEffectsCorrectWindAngle",
            coeffs_.lookupOrDefault("endEffectsCorrectWindAngle", false)
        );
        dict.add
        (
            "endEffectsCorrectAngleClamp",
            coeffs_.lookupOrDefault("endEffectsCorrectAngleClamp", -1)
        );
        dict.add
        (
            "velocitySampleRadius",
            coeffs_.lookupOrDefault("velocitySampleRadius", 0.0)
        );
        dict.add
        (
            "nVelocitySamples",
            coeffs_.lookupOrDefault("nVelocitySamples", 20)
        );
        if (coeffs_.found("dynamicStall"))
        {
            dictionary dsDict = coeffs_.subDict("dynamicStall");
            dsDict.add("chordLength", chordLength);
            dict.add("dynamicStall", dsDict);
        }
        if (coeffs_.found("meshFactor"))
        {
            scalar meshFactor = 2.0;
            coeffs_.lookup("meshFactor") >> meshFactor;
            dict.add("meshFactor", meshFactor);
        }
        if (coeffs_.found("dragFactor"))
        {
            scalar dragFactor = 1.0;
            coeffs_.lookup("dragFactor") >> dragFactor;
            dict.add("dragFactor", dragFactor);
        }
        if (coeffs_.found("chordFactor_"))
        {
            scalar chordFactor_ = 0.25;
            coeffs_.lookup("chordFactor_") >> chordFactor_;
            dict.add("chordFactor_", chordFactor_);
        }
        dictionary fcDict = coeffs_.subOrEmptyDict("flowCurvature");
        dict.add("flowCurvature", fcDict);
        bool writeElementPerf
        (
            coeffs_.lookupOrDefault("writeElementPerf", false)
        );
        bool writeElementPerfEnd
        (
            coeffs_.lookupOrDefault("writeElementPerfEnd", false)
        );
        dict.add("writePerf", writeElementPerf);
        dict.add("writePerfEnd", writeElementPerfEnd);

        if (debug)
        {
            Info<< "Creating actuatorLineElement: " << name << endl;
            Info<< "Geometry segment index: " << geometrySegmentIndex << endl;
            Info<< "Position: " << position << endl;
            Info<< "Chord length: " << chordLength << endl;
            Info<< "Chord direction (before pitching): " << chordDirection
                << endl;
            Info<< "Pitch (degrees): " << pitch << endl;
            Info<< "Span length: " << spanLength << endl;
            Info<< "Span direction: " << spanDirection << endl;
            Info<< "Profile name index: " << elementProfileIndex << endl;
            Info<< "Profile name: " << profileName << endl;
            Info<< "writePerf: " << writeElementPerf << endl;
            Info<< "writePerfEnd: " << writeElementPerfEnd << endl;
            Info<< "Root distance (nondimensional): " << rootDistance << endl;
        }

        actuatorLineElement* element = new actuatorLineElement
        (
            name, dict, mesh_
        );
        elements_.set(i, element);
        pitch = Foam::degToRad(pitch);
        elements_[i].pitch(pitch);
        elements_[i].setVelocity(initialVelocity);
    }
    if (endEffectsActive_)
    {
        elementChordLengths_.setSize(nElements_);
        elementRootDistances_.setSize(nElements_);

        forAll(elements_, i)
        {
            elementChordLengths_[i]   = elements_[i].chordLength();
            elementRootDistances_[i]  = elements_[i].rootDistance();
        }
    }
}


void Foam::fv::actuatorLineSource::writePerf()
{
    scalar time = mesh_.time().value();
    scalar totalArea = 0.0;
    scalar x = 0.0;
    scalar y = 0.0;
    scalar z = 0.0;
    scalar relVelMag = 0.0;
    scalar alphaDeg = 0.0;
    scalar alphaGeom = 0.0;
    scalar cl = 0.0;
    scalar cd = 0.0;
    scalar cm = 0.0;

    forAll(elements_, i)
    {
        scalar area = elements_[i].chordLength()*elements_[i].spanLength();
        totalArea += area;
        vector pos = elements_[i].position();
        x += pos[0]; y += pos[1]; z += pos[2];
        relVelMag += mag(elements_[i].relativeVelocity())*area;
        alphaDeg += elements_[i].angleOfAttack()*area;
        alphaGeom += elements_[i].angleOfAttackGeom()*area;
        cl += elements_[i].liftCoefficient()*area;
        cd += elements_[i].dragCoefficient()*area;
        cm += elements_[i].momentCoefficient()*area;
    }

    x /= nElements_; y /= nElements_; z /= nElements_;
    relVelMag /= totalArea;
    alphaDeg /= totalArea;
    alphaGeom /= totalArea;
    cl /= totalArea; cd /= totalArea; cm /= totalArea;

    // write time,x,y,z,rel_vel_mag,alpha_deg,alpha_geom_deg,cl,cd,cm
    stringBuffer_ << time << "," << x << "," << y << "," << z << ","
                << relVelMag << "," << alphaDeg << "," << alphaGeom << ","
                << cl << "," << cd << "," << cm << std::endl;

    // only write to file with writePerf_, writePerfEnd_ writes in destructor
    if (writePerf_ && outputFile_.is_open())
    {
        outputFile_ << stringBuffer_.str();
        stringBuffer_.str("");
        stringBuffer_.clear();
    }
}


void Foam::fv::actuatorLineSource::calcEndEffects()
{
    if (debug)
    {
        Info<< "Calculating end effects for " << name_ << endl;
    }

    // Decide how many stations to use for the linear system
    label nCalc = min(nElements_, maxNElementsEndEffects_);

    scalar pi = Foam::constant::mathematical::pi;
    List<scalar> c; // Chord lengths
    List<scalar> theta; // Span distance rescaled on [0, pi]
    List<scalar> A; // Fourier coefficients
    List<scalar> circulation;
    List<scalar> cl;
    List<scalar> factorsCalc;
    

    bool acceptable = false;

    while (nCalc > 4 && acceptable == false)   // never go below 4 stations
    {
        theta.setSize(nCalc);
        c.setSize(nCalc);
        circulation.setSize(nCalc);
        cl.setSize(nCalc);
        factorsCalc.setSize(nCalc);
        List<scalar> alpha(nCalc, 0.1);
        List<scalar> relVelMag(nCalc, 1.0);
        

        // Create lists from element parameters
        if (nCalc == nElements_)
        {
            forAll(elements_, n)
            {
                theta[n] = elements_[n].rootDistance()*pi;
                c[n] = elements_[n].chordLength();
                //~ alpha[n] = Foam::degToRad(elements_[n].angleOfAttackGeom());
                //~ relVelMag[n] = mag(elements_[n].relativeVelocityGeom());
            }
        }
        else
        {
            // Uniform spacing in root-distance, staying away from the tips
            const scalar dr = 1.0/(nCalc + 1);
            for (label n = 0; n < nCalc; ++n)
            {
                const scalar rootDist = (n + 1)*dr;
                theta[n] = rootDist*pi;

                // Linear interpolation of chord from the stored element data
                c[n] = interpolateUtils::interpolate1D
                (
                    rootDist,
                    elementRootDistances_,
                    elementChordLengths_
                );
            }
        }

        // resize D matrix
        simpleMatrix<scalar> D(nCalc, 0.0, 0.1);

        forAll(theta, i)
        {
            scalar n = i + 1;
            forAll(theta, m)
            {
                D[m][i] = 2.0*totalLength_/(pi*c[m])*sin(n*theta[m])
                        + n*sin(n*theta[m]) / sin(theta[m]);
            }
            D.source()[i] = alpha[i];
        }

        SVD svd(D, SMALL);                 // SMALL avoids zero singular values
        const scalarList& S = svd.S();
        const scalar cond = max(S)/(min(S) + VSMALL);

        if (debug)
        {
            Info<< "  nCalc = " << nCalc
                << ", estimated κ₂(D) = " << cond << endl;
        }

        if (cond < maxConditionNumberEndEffects_)
        {
            acceptable = true;
            A = D.solve();

            List<scalar> residual = D * A - D.source();
            scalar resNorm = Foam::sqrt(sum(residual*residual));
            if (resNorm > maxResidualEndEffects_)
            {
                acceptable = false;
                Info<< "End correction has residual " << resNorm
                    << ", rejecting solution with "
                    << nCalc << " elements" << endl;
            }

            forAll(theta, m)
            {
                scalar sumA = 0.0;
                forAll(theta, i)
                {
                    scalar n = i + 1;
                    sumA += A[i]*sin(n*theta[m]);
                }
                circulation[m] = 2*totalLength_*relVelMag[m]*sumA;
                cl[m] = circulation[m]/(0.5*c[m]*relVelMag[m]);
            }
            factorsCalc = cl/(2.0 * constant::mathematical::pi * alpha);
            if (normalizeEndEffects_)
            {
                // Root-distance coordinates of the collocation stations
                // (theta runs from 0 -> π, so rootDist = theta/π)
                List<scalar> rootDist(nCalc);
                forAll(theta, i)
                {
                    rootDist[i] = theta[i] / constant::mathematical::pi;
                }

                const scalar low  = ignoreLowValuesEndEffects_;
                const scalar high = 1.0 - ignoreHighValuesEndEffects_;

                // Collect indices inside the central interval [low, high]
                DynamicList<label> centralIdx;
                forAll(rootDist, i)
                {
                    if (rootDist[i] >= low && rootDist[i] <= high)
                    {
                        centralIdx.append(i);
                    }
                }

                scalar normValue = 0.0;

                if (centralIdx.size() > 0)
                {
                    // Normal case: take the maximum factor inside the interval
                    normValue = factorsCalc[centralIdx[0]];
                    forAll(centralIdx, j)
                    {
                        normValue = max(normValue, factorsCalc[centralIdx[j]]);
                    }
                }
                else
                {
                    // Interval too narrow – fall back to the station closest
                    // to the middle of the requested interval
                    const scalar mid = 0.5*(low + high);
                    label closest = 0;
                    scalar minDist = mag(rootDist[0] - mid);

                    for (label i = 1; i < nCalc; ++i)
                    {
                        const scalar d = mag(rootDist[i] - mid);
                        if (d < minDist)
                        {
                            minDist = d;
                            closest = i;
                        }
                    }
                    normValue = factorsCalc[closest];

                    if (debug)
                    {
                        Info<< "normalizeEndEffects: no station inside ["
                            << low << ", " << high
                            << "], using station " << closest
                            << " (rootDist = " << rootDist[closest]
                            << ")" << endl;
                    }
                }

                // Protect against a zero / negative normValue
                if (normValue < SMALL)
                {
                    WarningInFunction
                        << "Normalisation value for end-effect "
                        << "factors is too small (" << normValue
                        << ").  Skipping normalisation." << endl;
                }
                else
                {
                    factorsCalc = factorsCalc/normValue;
                }
            }
            // Optional safety: never allow factors > 1
            forAll(factorsCalc, i)
            {
                factorsCalc[i] = min(factorsCalc[i], 1.0);
            }
            // Additional checks to avoid incorrect solutions
            // No negative factors
            if (min(factorsCalc) < 0.0)
            {
                acceptable = false;
                Info<< "End correction has negative factors, rejecting "
                    << "solution with " << nCalc << " elements" << endl;
                    
            }
        }
        else
        {
            Info<< "End correction has condition number " << cond
                << ", rejecting solution with "
                << nCalc << " elements" << endl;
        }
        if (acceptable == false)
        {
            // Too ill-conditioned → try fewer stations
            nCalc = max(4, nCalc*3/4);   // reduce by ~25%
            Info<< "Reducing number of elements for end effect correction to "
                << nCalc << endl;
        }
        if (debug == 2)
        {
            Info<< "D.source: " << D.source() << endl;
            Info<< "D: " << D << endl;
        }
    }

    if (acceptable == false)
    {
        WarningInFunction
            << "Could not find a well-conditioned Prandtl matrix "
            << "even with nCalc = " << nCalc
            << ". Using the last attempt." << endl;
        // (A and factorsCalc are left from the last iteration)
    }

    // Set endEffectFactor for all elements
    
    if (nCalc == nElements_)
    {
        forAll(elements_, i)
        {
            elements_[i].setEndEffectFactor(factorsCalc[i]);
        }
    }
    else
    {
        // Interpolate onto the true element locations

        // Build extended tables that force zero factor at the tips
        List<scalar> rootDistCalc(nCalc + 2);
        List<scalar> factorCalc  (nCalc + 2);

        rootDistCalc[0]           = 0.0;
        factorCalc[0]             = 0.0;
        rootDistCalc[nCalc + 1]   = 1.0;
        factorCalc[nCalc + 1]     = 0.0;

        for (label i = 0; i < nCalc; i++)
        {
            rootDistCalc[i + 1] = theta[i]/pi;
            factorCalc[i + 1]   = factorsCalc[i];
        }

        // Now map onto every original element
        forAll(elements_, i)
        {
            const scalar f = interpolateUtils::interpolate1D
            (
                elementRootDistances_[i],
                rootDistCalc,
                factorCalc
            );
            elements_[i].setEndEffectFactor(f);
        }
    }

    if (debug == 2)
    {
        Info<< "Debug output from actuatorLineSource::calcEndEffects:" << endl;
        Info<< "theta: " << theta << endl;
        Info<< "A: " << A << endl;
        Info<< "c: " << c << endl;
        Info<< "cl: " << cl << endl;
        Info<< "factors:" << factorsCalc << endl;
    }
}


void Foam::fv::actuatorLineSource::harmonicPitching()
{
    // Pitch the actuator line if time has changed
    scalar t = mesh_.time().value();
    if (t != lastMotionTime_)
    {
        scalar omega = reducedFreq_*2*mag(freeStreamVelocity_)/chordLength_;
        scalar dt = mesh_.time().deltaT().value();
        scalar deltaPitch = degToRad(pitchAmplitude_)*(sin(omega*t)
                          - sin(omega*(t - dt)));
        pitch(deltaPitch);
        lastMotionTime_ = t;
    }
}

//- virtual function setupPositions is protected,
// but findCells needs to be public
void Foam::fv::actuatorLineSource::setupPositions(bool includeRing)
{
    // If harmonic pitching is active, do harmonic pitching
    if (harmonicPitchingActive_)
    {
        harmonicPitching();
    }
    findCells();
}

//- same as above, calculateForces is the virtual function
// that has to be implemented
void Foam::fv::actuatorLineSource::calculateForces()
{
    calculateElementForces();
}

//- same as above, need to implement allocateAL
void Foam::fv::actuatorLineSource::allocateAL()
{
    allocateInfluenceCells(1, false);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineSource::actuatorLineSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    actuatorModelBase(name, modelType, dict, mesh),
    force_(vector::zero),
    writePerf_(coeffs_.lookupOrDefault("writePerf", false)),
    writePerfEnd_(coeffs_.lookupOrDefault("writePerfEnd", false)),
    lastMotionTime_(mesh.time().value()),
    endEffectsActive_(false),
    elementChordLengths_(0),
    elementRootDistances_(0)
{
    read(dict_);
    createElements();
    if (writePerf_ || writePerfEnd_)
    {
        createOutputFile();
    }

    // Calculate end effects
    if (endEffectsActive_)
    {
        calcEndEffects();
    }
    int precision = 6;
    if (mesh_.time().controlDict().found("writePrecision"))
    {
        mesh_.time().controlDict().lookup("writePrecision") >> precision;
    }
    stringBuffer_.precision(precision);

    // an actuatorLineSource will only have itself in actuatorLines
    actuatorLines_.setSize(1);
    actuatorLines_[0] = this;
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fv::actuatorLineSource::~actuatorLineSource()
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

void Foam::fv::actuatorLineSource::printCoeffs() const
{
    // Print turbine properties
    Info<< "Actuator line properties:" << endl;
    Info<< "Profile data:" << endl;
    Info<< profileData_ << endl;
    Info<< "First item of element geometry:" << endl;
    Info<< elementGeometry_[0] << endl;
}


void Foam::fv::actuatorLineSource::rotate
(
    vector rotationPoint,
    vector axis,
    scalar radians
)
{
    forAll(elements_, i)
    {
        elements_[i].rotate(rotationPoint, axis, radians, true);
    }
}


void Foam::fv::actuatorLineSource::pitch(scalar radians)
{
    forAll(elements_, i)
    {
        elements_[i].pitch(radians);
    }
}


void Foam::fv::actuatorLineSource::pitch(scalar radians, scalar chordFraction)
{
    forAll(elements_, i)
    {
        elements_[i].pitch(radians, chordFraction);
    }
}


void Foam::fv::actuatorLineSource::translate(vector translationVector)
{
    forAll(elements_, i)
    {
        elements_[i].translate(translationVector);
    }
}


void Foam::fv::actuatorLineSource::setSpeed
(
    vector point,
    vector axis,
    scalar omega
)
{
    forAll(elements_, i)
    {
        elements_[i].setSpeed(point, axis, omega);
    }
}


void Foam::fv::actuatorLineSource::scaleVelocity(scalar scale)
{
    forAll(elements_, i)
    {
        elements_[i].scaleVelocity(scale);
    }
}


void Foam::fv::actuatorLineSource::setOmega(scalar omega)
{
    forAll(elements_, i)
    {
        elements_[i].setOmega(omega);
    }
}


void Foam::fv::actuatorLineSource::setCustomTime
(
    scalar time,
    scalar deltaT,
    bool useCustomTime
)
{
    forAll(elements_, i)
    {
        elements_[i].setCustomTime(time, deltaT, useCustomTime);
    }
}


void Foam::fv::actuatorLineSource::calcInfluenceEpsilon
(
    scalar dragCoefficient
)
{
    forAll(elements_, i)
    {
        elements_[i].calcInfluenceEpsilon(dragCoefficient);
    }
}

void Foam::fv::actuatorLineSource::allocateInfluenceCells
(
    label count,
    bool cacheInteractions_
)
{
    forAll(elements_, i)
    {
        elements_[i].allocateInfluenceCells(count, cacheInteractions_);
    }
}

void Foam::fv::actuatorLineSource::constructInfluenceCellList
(
    label azimuthIndex,
    labelList& globalToLocal,
    label& nActive
)
{
    forAll(elements_, i)
    {
        elements_[i].constructInfluenceCellList
        (
            azimuthIndex,
            globalToLocal,
            nActive
        );
    }
}

void Foam::fv::actuatorLineSource::setCompactFields
(
    vectorField& activePositions,
    vectorField& activeForceField
)
{
    forAll(elements_, i)
    {
        elements_[i].setCompactFields
        (
            activePositions,
            activeForceField
        );
    }
}

void Foam::fv::actuatorLineSource::setAzimuthIndex
(
    label azimuthIndex,
    bool clearBuffer
)
{
    forAll(elements_, i)
    {
        elements_[i].setAzimuthIndex(azimuthIndex, clearBuffer);
    }
    if (clearBuffer &&  azimuthIndex == 0)
    {
        stringBuffer_.str("");
        stringBuffer_.clear();
    }
}

void Foam::fv::actuatorLineSource::findCells(bool includeRing)
{
    forAll(elements_, i)
    {
        elements_[i].findCells(includeRing);
    }
}

void Foam::fv::actuatorLineSource::calculateElementForces()
{
    forAll(elements_, i)
    {
        elements_[i].calculateForce();
    }
    // Write performance to file
    if (Pstream::master() && (writePerf_ || writePerfEnd_))
    {
        writePerf();
    }
}

const Foam::vector Foam::fv::actuatorLineSource::force()
{
    // recalculate instead of using force_ to get correct azimuthIndex
    Foam::vector force = vector::zero;
    
    forAll(elements_, i)
    {
        force += elements_[i].force();
    }
    return force;
}

PtrList<Foam::fv::actuatorLineElement>& Foam::fv::actuatorLineSource::elements()
{
    return elements_;
}


Foam::vector Foam::fv::actuatorLineSource::moment(vector point)
{
    vector moment(vector::zero);
    forAll(elements_, i)
    {
        moment += elements_[i].moment(point);
    }

    if (debug)
    {
        Info<< "Moment on " << name_ << " about " << point << ": " << moment
            << endl;
    }

    return moment;
}


void Foam::fv::actuatorLineSource::addForce
(
    volVectorField& forceField,
    scalar scale,
    bool compressible
)
{
    // Zero the total force vector
    force_ = vector::zero;
    
    forAll(elements_, i)
    {
        elements_[i].addForce(forceField, scale, compressible);
        force_ += elements_[i].force();
    }

    if (printPerf_)
    {
        Info<< "Force (per unit density) on " << name_ << ": "
            << endl << force_ << endl << endl;
    }
}

void Foam::fv::actuatorLineSource::addTurbulence
(
    fvMatrix<scalar>& eqn,
    const word fieldName
)
{
    forAll(elements_, i)
    {
        elements_[i].addTurbulence(eqn, fieldName);
    }
}

// ************************************************************************* //
