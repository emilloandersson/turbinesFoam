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

#include "axialFlowTurbineADSource.H"
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
    defineTypeNameAndDebug(axialFlowTurbineADSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        axialFlowTurbineADSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::axialFlowTurbineADSource::axialFlowTurbineADSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    axialFlowTurbineALSource(name, modelType, dict, mesh)
{
    read(dict);
    customTime_ = mesh.time().value();
    rotateAD(true);

    // override the nBlades value for the end effects calculation 
    // if bladeMultiplier is used
    effectiveNBlades_ = bladeMultiplier_*nBlades_; 
    forAll(blades_, i)
    {
        blades_[i].setApplyForce(false);
    }
    
    if (hasHub_)
    {
        hub_->setApplyForce(false);
    }
    
    if (hasTower_)
    {
        tower_->setApplyForce(false);
    }

    if (hasNacelle_)
    {
        nacelle_->setApplyForce(false);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::axialFlowTurbineADSource::~axialFlowTurbineADSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::axialFlowTurbineADSource::rotateAD(bool updateOnly)
{
    
    scalar radians = 2*mathematical::pi/divisions_;
    scalar deltaT = radians/omega_;

    //updateOnly is intended for first step only to set custom speed
    if (updateOnly == false)
    {
        customTime_ += deltaT;
        rotate(radians);
        angleDeg_ += radToDeg(radians);
        //lastRotationTime_ = time_.value();
    }
    updateTSROmega();
    
    //Info << "rotateAD called for time " << time_.value() 
    //     << " custom time: " << customTime_ << endl;
    forAll(blades_, i)
    {
        blades_[i].setCustomTime(customTime_, deltaT);
    }
    
    if (hasHub_)
    {
        hub_->setCustomTime(customTime_, deltaT);
    }
    
    if (hasTower_)
    {
        tower_->setCustomTime(customTime_, deltaT);
    }

    if (hasNacelle_)
    {
        nacelle_->setCustomTime(customTime_, deltaT);
    }
}

void Foam::fv::axialFlowTurbineADSource::rotate(scalar radians)
{
    if (debug)
    {
        Info<< "Rotating " << name_ << " " << radians << " radians"
            << endl << endl;
    }

    forAll(blades_, i)
    {
        blades_[i].rotate(origin_, axis_, radians);
        blades_[i].setSpeed(origin_, axis_, omega_);
    }

    if (hasHub_)
    {
        hub_->rotate(origin_, axis_, radians);
        hub_->setSpeed(origin_, axis_, omega_);
    }
}


void Foam::fv::axialFlowTurbineADSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // tower and nacelle are not rotating,
    // so we only need to calculate the force field once
    if (hasTower_)
    {
        // Add source for tower actuator line
        tower_->addSup(eqn, fieldI);
    }

    if (hasNacelle_)
    {
        // Add source for tower actuator line
        nacelle_->addSup(eqn, fieldI);
    }
    // code can run extra revolutions to make dynamic stall converge
    for (int currentLoop = 0; currentLoop < dynStallLoop_; currentLoop++)
    {
        // forceField_ should be the average during one revolution here
        forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
        for (int innerStep = 0; innerStep < divisions_; innerStep++)
        {
            // Zero out force vector and field
            force_ *= 0;

            // Check dimensions of force field and correct if necessary
            if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
            {
                forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
            }

            // Create local moment vector
            vector moment(vector::zero);
            
            if (endEffectsActive_ and endEffectsModel_ != "liftingLine")
            {
                // Calculate end effects based on current velocity field
                calcEndEffects();
            }

            // Add source for blade actuator lines
            forAll(blades_, i)
            {
                blades_[i].addSup(eqn, fieldI);
                forceField_ +=
                    (bladeMultiplier_/divisions_)*blades_[i].forceField();
                //Info<< "Added blade" << endl;
                force_ += bladeMultiplier_*blades_[i].force();
                bladeMoments_[i] = blades_[i].moment(origin_);
                moment += bladeMultiplier_*bladeMoments_[i];
            }

            if (hasHub_)
            {
                // Add source for hub actuator line
                hub_->addSup(eqn, fieldI);
                forceField_ += (1.0/divisions_)*hub_->forceField();
                force_ += hub_->force();
                moment += hub_->moment(origin_);
            }

            if (hasTower_)
            {
                // Add source for tower actuator line
                //tower_->addSup(eqn, fieldI);
                forceField_ += (1.0/divisions_)*tower_->forceField();
                if (includeTowerDrag_)
                {
                    force_ += tower_->force();
                }
            }

            if (hasNacelle_)
            {
                // Add source for tower actuator line
                //nacelle_->addSup(eqn, fieldI);
                forceField_ += (1.0/divisions_)*nacelle_->forceField();
                if (includeNacelleDrag_)
                {
                    force_ += nacelle_->force();
                }
            }

            // only needed to to do calculations for the actual loop,
            // not the dummy loops that are used to make dynamic stall converge
            if (currentLoop == dynStallLoop_ - 1)
            {
                // Torque is the projection of the moment from
                // all blades on the axis
                torque_ = moment & axis_;

                torqueCoefficient_ =
                    torque_/(0.5*frontalArea_*rotorRadius_
                    * magSqr(freeStreamVelocity_));
                powerCoefficient_ = torqueCoefficient_*tipSpeedRatio_;
                dragCoefficient_ = 
                    force_ & freeStreamDirection_
                    / (0.5*frontalArea_*magSqr(freeStreamVelocity_));


                // Print performance to terminal
                //printPerf();

                // Write performance data
                // Note this will write multiples if there are
                // multiple PIMPLE loops
                if (Pstream::master())
                {
                    writePerf();
                }
            }
            rotateAD();
        }
    }
    eqn += forceField_;
}


void Foam::fv::axialFlowTurbineADSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // code can run extra revolutions to make dynamic stall converge
    for (int currentLoop = 0; currentLoop < dynStallLoop_; currentLoop++)
    {
        // forceField_ should be the average during one revolution here
        forceField_ *=
            dimensionedScalar("zero", forceField_.dimensions(), 0.0);
        for (int innerStep = 0; innerStep < divisions_; innerStep++)
        {
            // Zero out force vector and field
            force_ *= 0;

            // Create local moment vector
            vector moment(vector::zero);

            // Check dimensions of force field and correct if necessary
            if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
            {
                forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
            }
            
            if (endEffectsActive_ and endEffectsModel_ != "liftingLine")
            {
                // Calculate end effects based on current velocity field
                calcEndEffects();
            }

            // Add source for blade actuator lines
            forAll(blades_, i)
            {
                blades_[i].addSup(rho, eqn, fieldI);
                forceField_ +=
                    (bladeMultiplier_/divisions_)*blades_[i].forceField();
                force_ += bladeMultiplier_*blades_[i].force();
                bladeMoments_[i] = blades_[i].moment(origin_);
                moment += bladeMultiplier_*bladeMoments_[i];
            }

            if (hasHub_)
            {
                // Add source for hub actuator line
                hub_->addSup(rho, eqn, fieldI);
                forceField_ += (1.0/divisions_)*hub_->forceField();
                force_ += hub_->force();
                moment += hub_->moment(origin_);
            }

            if (hasTower_)
            {
                // Add source for tower actuator line
                tower_->addSup(rho, eqn, fieldI);
                forceField_ += (1.0/divisions_)*tower_->forceField();
                if (includeTowerDrag_)
                {
                    force_ += tower_->force();
                }
            }

            if (hasNacelle_)
            {
                // Add source for tower actuator line
                nacelle_->addSup(rho, eqn, fieldI);
                forceField_ += (1.0/divisions_)*nacelle_->forceField();
                if (includeNacelleDrag_)
                {
                    force_ += nacelle_->force();
                }
            }
            
            if (currentLoop == dynStallLoop_ - 1)
            {
                // Torque is the projection of the moment from all blades
                // on the axis
                torque_ = moment & axis_;

                scalar rhoRef;
                coeffs_.lookup("rhoRef") >> rhoRef;
                torqueCoefficient_ = 
                    torque_/(0.5*rhoRef*frontalArea_*rotorRadius_
                    * magSqr(freeStreamVelocity_));
                powerCoefficient_ = torqueCoefficient_*tipSpeedRatio_;
                dragCoefficient_ = 
                    force_ & freeStreamDirection_
                    / (0.5*rhoRef*frontalArea_*magSqr(freeStreamVelocity_));

                // Print performance to terminal
                //printPerf();

                // Write performance data -- note this will write multiples
                // if there are multiple PIMPLE loops
                if (Pstream::master())
                {
                    writePerf();
                }
            }
            rotateAD();
        }
    }
    eqn += forceField_;
}


void Foam::fv::axialFlowTurbineADSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // code can run extra revolutions to make dynamic stall converge
    for (int currentLoop = 0; currentLoop < dynStallLoop_; currentLoop++)
    {
        // forceField_ should be the average during one revolution here
        fvMatrix<scalar> kField(eqn.psi(), eqn.dimensions());
        kField *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
        for (int innerStep = 0; innerStep < divisions_; innerStep++)
        {
            if (endEffectsActive_ and endEffectsModel_ != "liftingLine")
            {
                // Calculate end effects based on current velocity field
                calcEndEffects();
            }
            
            // Add scalar source term from blades
            forAll(blades_, i)
            {
                blades_[i].addSup(kField, fieldI);
            }

            if (hasHub_)
            {
                // Add source for hub actuator line
                hub_->addSup(kField, fieldI);
            }

            if (hasTower_)
            {
                // Add source for tower actuator line
                tower_->addSup(kField, fieldI);
            }

            if (hasNacelle_)
            {
                // Add source for nacelle actuator line
                nacelle_->addSup(kField, fieldI);
            }
            rotateAD();
        }
        eqn += (bladeMultiplier_/divisions_)*kField;
    }
}


bool Foam::fv::axialFlowTurbineADSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        //crossFlowTurbineALSource::read(dict);

        // Get number of divisions
        divisions_ = coeffs_.lookupOrDefault("divisions", 180);
        
        // Get number of divisions
        dynStallLoop_ = coeffs_.lookupOrDefault("dynStallLoop", 1);
        
        // Get blade multiplier
        bladeMultiplier_ = coeffs_.lookupOrDefault("bladeMultiplier", 1);
        
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
