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

Class
    Foam::fv::actuatorLineElement

Description
    An object, e.g. an airfoil, that imparts a force on the flow field
    calculated via lift and drag coefficients

SourceFiles
    actuatorLineElement.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorLineElement_H
#define actuatorLineElement_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "vector.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "dynamicStallModel.H"
#include "interpolationCellPoint.H"
#include "profileData.H"
#include "addedMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class actuatorLineElement Declaration
\*---------------------------------------------------------------------------*/

class actuatorLineElement
{

protected:

    // Protected data

        //- Coefficients dictionary
        const dictionary dict_;

        //- Name of element
        const word name_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Mesh bounding box
        boundBox meshBoundBox_;

        //- Chord direction -- pointing from trailing to leading edge
        List<vector> chordDirection_;

        //- Chord length
        scalar chordLength_;

        //- Span direction -- direction about which element will pitch
        //  Can be obtained by sweeping chord direction into planform normal
        List<vector> spanDirection_;

        //- Element span length, for calculating force
        scalar spanLength_;

        //- Element planform normal vector -- direction when the opposite of
        //  the chord direction is swept into the span direction
        vector planformNormal_;

        //- Element chord reference direction
        List<vector> chordRefDirection_;

        //- Location of element's quarter-chord
        List<vector> position_;

        //- Velocity of the element
        List<vector> velocity_;

        //- Free stream velocity
        vector freeStreamVelocity_;

        //- Free stream direction
        vector freeStreamDirection_;

        //- Fluid force on element -- vector that is perpendicular to
        //  spanDirection,
        List<vector> forceVector_;

        //- Reference density for incompressible case
        scalar rhoRef_;

        //- Density at element location for compressible case
        List<scalar> localRho_;

        //- Density at element location for compressible case
        List<scalar> localMu_;

        //- Inflow velocity
        List<vector> inflowVelocity_;

        //- Epsilon for force projection width
        List<scalar> epsilon_;

        //- Previously found cell, to prevent unnecessary lookups
        List<point> previousLocation_;

        //- Previously found ring cell (actuator disc only)
        List<List<point>> previousRingLocation_;

        //- for acutatorDisc, have previous positions stored for fast lookup
        List<label> centerCellI_;

        //- for acutatorDisc, have previous positions stored for fast lookup
        List<label> centerProcI_;

        //- for actuatorDisc, ring sampling
        List<List<label>> ringCellI_;

        //- for actuatorDisc, ring sampling
        List<List<label>> ringProcI_;

        //- corresponding velocities to ringCellI_
        List<List<vector>> velocitiesRing_;

        //- Relative flow velocity
        List<vector> relativeVelocity_;

        //- Geometric relative flow velcity
        vector relativeVelocityGeom_;

        //- Angle of attack (degrees)
        scalar angleOfAttack_;

        //- Geometric angle of attack (degrees, no flow curvature correction)
        scalar angleOfAttackGeom_;

        //- Lift coefficient
        List<scalar> liftCoefficient_;

        //- Drag coefficient;
        List<scalar> dragCoefficient_;

        //- 1/4 chord moment coefficient
        List<scalar> momentCoefficient_;

        //- Profile name
        word profileName_;

        //- Profile data (encapsulates coefficients)
        profileData profileData_;

        //- Dynamic stall model
        autoPtr<dynamicStallModel> dynamicStall_;

        //- Switch for applying dynamic stall model
        bool dynamicStallActive_;

        //- Kinematic viscosity (for calculating Re)
        scalar nu_;

        //- Chord Reynolds number based on relative velocity
        scalar Re_;

        //- Angular velocity for flow curvature correction
        scalar omega_;

        //- Nondimensional chord mount distance
        scalar chordMount_;

        //- Switch for applying flow curvature correction
        bool flowCurvatureActive_;

        //- Flow curvature model name
        word flowCurvatureModelName_;

        //- Leading edge velocity vector
        vector velocityLE_;

        //- Trailing edge velocity vector
        vector velocityTE_;

        //- Switch for writing performance
        bool writePerf_;

        //- Switch for writing performance only at final step
        bool writePerfEnd_;

        //- Output file stream (Foam::OFStream does not work with destructor)
        std::ofstream outputFile_;

        //- For WritePerfEnd, save to stream instead
        std::stringstream stringBuffer_;

        //- Nondimensional distance from actuator line root
        scalar rootDistance_;
        
        //- cone angle of the blade for use by curvature correction
        scalar cone_;

        //- End effect correction factor [0, 1]
        scalar endEffectFactor_;

        //- Correct angle of attack instead of CL
        bool endEffectsCorrectAngleOfAttack_;

        //- Correct angle of relative wind
        bool endEffectsCorrectAngleOfRelWind_;

        //- Limit how much the angle of the relative wind can be corrected
        scalar endEffectsCorrectAngleClamp_;

        //- The angle to correct the wind with
        scalar angleOfAttackDiff_;

        //- Switch for added mass correction
        bool addedMassActive_;

        //- Added mass model
        addedMassModel addedMass_;

        //- Radius of circle used to sample velocities
        scalar velocitySampleRadius_;

        //- Number of elements used to sample velocities
        label nVelocitySamples_;

        //- List of cell indices for the actuatorDisc loop
        List<List<label>> influenceCells_;

        //- List of cell indices for the actuatorDisc loop
        List<label> tmpInfluenceCells_;

        //- ActuatorDisc - pointer to the compact field of the owner class
        vectorField *activePositionsPtr_;

        //- ActuatorDisc - pointer to the compact field of the owner class
        vectorField *activeForceFieldPtr_;

        //- The current azimuthIndex to use for the influenceCells
        label azimuthIndex_;

        //- Factor for determining epsilon based on chord lengt
        scalar chordFactor_;

        //- Factor for determining epsilon based on drag coefficient
        scalar dragFactor_;

        //- Factor for determining epsilon based on cell volume
        scalar meshFactor_;


    // Protected Member Functions

        //- Rotate a vector
        void rotateVector
        (
            vector& vectorToRotate,
            vector rotationPoint,
            vector axis,
            scalar radians
        );

        //- Find cell containing location
        //label findCell(const point& location);

        //- Lookup force coefficients
        void lookupCoefficients();

        //- Calculate projection width epsilon
        void calcProjectionEpsilon();

        //- Correct for flow curvatue
        void correctFlowCurvature(scalar& angleOfAttackRad);

        //- Multiply force vector by local density
        void multiplyForceRho();

        //- Apply force field based on force vector
        void applyForceField
        (
            volVectorField& forceField,
            scalar scale
        );


        //- Get inflow velocity
        void calculateRingVelocity();

        //- Create the performance output file
        virtual void createOutputFile();

        //- Write performance to CSV
        void writePerf();


public:

    //- Runtime type information
    TypeName("actuatorLineElement");


        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            actuatorLineElement,
            dictionary,
            (
                const dictionary& dict,
                const word& modelName
            ),
            (dict, modelName)
        );


    // Selectors

        //- Return a reference
        static autoPtr<actuatorLineElement> New(const dictionary& dict);


    //- Constructor
    actuatorLineElement
    (
        const word& name,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~actuatorLineElement();


    // Member functions

        // Access

            //- Return const access to the element name
            const word& name() const;

            //- Return the element chord length
            const scalar& chordLength() const;

            //- Return the element span length
            const scalar& spanLength();

            //- Return the element location
            const vector& position();

            //- Return the element velocity
            const vector& velocity();

            //- Return force vector from element
            const vector& force();

            //- Return relative velocity
            const vector& relativeVelocity();

            //- Return geometric relative velocity
            const vector& relativeVelocityGeom();

            //- Return element angle of attack in degrees
            const scalar& angleOfAttack();

            //- Return element geometric angle of attack in degrees
            const scalar& angleOfAttackGeom();

            //- Return lift coefficient
            const scalar& liftCoefficient();

            //- Return drag coefficient
            const scalar& dragCoefficient();

            //- Return moment coefficient
            const scalar& momentCoefficient();

            //- Return nondimensional distance from actuator line root
            const scalar& rootDistance();

            //- Return tangential reference coefficient
            scalar tangentialRefCoefficient();
            
            //- Return normal reference coefficient
            scalar normalRefCoefficient();

            //- Return tangential reference force (per unit span)
            scalar tangentialRefForce();
            
            //- Return normal reference force (per unit span)
            scalar normalRefForce();

            //- Return the inflow velocity angle w.r.t. the chord reference
            //  direction (degrees)
            //  E.g., for an axial-flow turbine this would be the angle of
            //  attack plus twist
            scalar inflowRefAngle();

            //- Return inflow velocity
            vector& inflowVelocity();

            //- Return inflow velocity
            const vector& inflowVelocity() const;

            //- Return epsilon for force projection width
            const scalar& epsilon() const;

            //- Return the number of locations for velocity sampling
            label locationCount() const;

            //- Return the number of epsilon values stored
            label epsilonCount() const;

            //- Return chord factor for determining epsilon
            const scalar& chordFactor() const;

            //- Return chord factor for determining epsilon
            const scalar& dragFactor() const;

            //- Return chord factor for determining epsilon
            const scalar& meshFactor() const;
            
            //- collect the locations of the actuator line element
            // for velocity sampling
            void collectLocationData
            (
                List<label>& globalCellI,
                List<label>& globalProcI,
                List<point>& globalLocations,
                label& index,
                bool includeRing = true
            ) const;

            //- distribute the locations of the actuator line element
            // and sampled velocity
            void distributeVelocityData
            (
                const List<label>& globalCellI,
                const List<label>& globalProcI,
                const List<vector>& globalVelocities,
                label& index
            );

            //- Distribute only centerCellI for the inital setup
            void distributeCenterCellI
            (
                const List<label>& globalCellI,
                const List<label>& globalProcI,
                label& index
            );

            //- collect epsilon of the actuator line elements
            void collectEpsilonData
            (
                List<scalar>& globalEpsilon,
                label& index
            ) const;

            //- Calculate projection width epsilon for influence cells
            void calcInfluenceEpsilon(scalar dragCoefficient);


            //- distribute epsilon of the actuator line elements
            void distributeEpsilonData
            (
                List<scalar>& globalEpsilon,
                label& index
            );

            //- distribute rho of the actuator line elements
            void distributeRhoData
            (
                List<scalar>& globalRho,
                label& index
            );

            //- distribute rho of the actuator line elements
            void distributeMuData
            (
                List<scalar>& globalMu,
                label& index
            );
        // Manipulation

            //- Translate the element's location
            void translate(vector translationVector);

            //- Rotate the element about some axis
            void rotate
            (
                vector rotationPoint,
                vector axis,
                scalar radians,
                bool rotateVelocity
            );

            //- Pitch the element about some fraction of the chord
            void pitch(scalar radians, scalar chordFraction=0.25);

            //- Set velocity of element
            void setVelocity(vector velocity);

            //- Set element speed (velocity magnitude)
            void setSpeed(scalar speed);

            //- Set element speed based on rotation
            void setSpeed(vector point, vector axis, scalar omega);

            //- Scale the velocity
            void scaleVelocity(scalar scale);

            //- Set dynamic stall active
            void setDynamicStallActive(bool active);

            //- Set omega for flow curvature correction
            void setOmega(scalar omega);

            //- Set end effect correction factor
            void setEndEffectFactor(scalar factor);

            //- Set velocity sample radius
            void setVelocitySampleRadius(scalar radius);

            //- Set number of velocity samples
            void setNVelocitySamples(label nSamples);
            
            //- Set custom time value for actuator disc model, where
            //  the internal time is not the same as the simulation time
            void setCustomTime
            (
                scalar time,
                scalar deltaT,
                bool useCustomTime = true
            );

            //- allocate the list of influence cells (for actuator disc)
            void allocateInfluenceCells(label count, bool cacheInteractions);

            //- Create the list of influence cells for one angle
            void constructInfluenceCellList
            (
                label azimuthIndex,
                labelList& globalToLocal,
                label& nActive
            );

            //- For actuatorDisc, set the compact position and force fields
            void setCompactFields
            (
                vectorField& activePositions,
                vectorField& activeForceField
            );

            //- allocate the list of influence cells (for actuator disc)
            void setAzimuthIndex(label azimuthIndex, bool clearBuffer);

            //- Check if one can find the cell in
            // the nearby positions of the last cell
            label findNearbyCell(
                const point& location,
                const label previousCell
            );

            //- Find cells that can be determined without parallel interaction
            void findCells(bool includeRing=true);


        // Evaluation

            //- Calculate forces
            void calculateForce();

            //- Read coefficient data
            void read();

            //- Calculate and return moment about specified point
            vector moment(vector point);


        // Source term addition

            //- Source term to momentum equation
            void addForce
            (
                volVectorField& force,
                scalar scale,
                bool compressible
            );

            //- Add source term to turbulence quantity
            virtual void addTurbulence(fvMatrix<scalar>& eqn, word fieldName);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
