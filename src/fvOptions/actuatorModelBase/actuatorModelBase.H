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
    Foam::fv::actuatorModelBase

Description
    Base class for all actuator models, with focus on handling parallel
    computing efficiently

SourceFiles
    actuatorModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorModelBase_H
#define actuatorModelBase_H

#include "cellSetOption.H"
//#include "NamedEnum.H"
#include "vector.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "interpolationCellPoint.H"
#include "OFstream.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

class actuatorLineSource;   // forward declaration

/*---------------------------------------------------------------------------*\
                      Class actuatorModelBase  Declaration
\*---------------------------------------------------------------------------*/

class actuatorModelBase
:
    public cellSetOption
{

protected:

    // Protected data

        //- Mesh bounding box
        boundBox meshBoundBox_;

        //- Save force field from this particular object
        bool writeForceField_;

        //- Pointer to force field (if used) (initialized to null)
        autoPtr<volVectorField> forceFieldPtr_;

        //- for mesh searches
        autoPtr<meshSearch> searchEnginePtr;

        //- Sample velocities
        List<vector> velocities_;

        //- point locations
        List<vector> locations_;

        //- List of cell indices
        List<label> cellI_;

        //. List of processor ownership
        List<label> procI_;
        
        //- List of epsilon values for each actuator line element
        List<scalar> epsilon_;

        //- List of density values for each actuator line element
        List<scalar> rho_;

        //- List of mu values for each actuator line element (compressible)
        List<scalar> mu_;
        
        //- All actuator lines belonging to this object (for fast iterations)
        List<actuatorLineSource*> actuatorLines_;

        //- Pointer to the density field (for compressible flow)
        mutable const volScalarField* rhoPtr_;

        //- Pointer to the dynamic viscocity field (for compressible flow)
        mutable const volScalarField* muPtr_;

        //- Flag to indicate if the object has been initialized
        bool initialized_;

        //- Last time the inflow velocities were calculated
        scalar lastTime_;

        //- First field to be called from addSup,
        // to avoid unnecessary calculations
        label firstField_;

        //- Print status to terminal
        bool printPerf_;

        //- name for force field
        word forceFieldName_;

        //- Flag to write fields at final time step
        bool writeFieldsAtEnd_;

        //- Make sure that final fields are written once
        bool finalFieldsWritten_;

    // Protected Member Functions

        //- Gather data, calculate velocity and distribute
        void calculateALData(label fieldI);

        // Gather epsilon, reduce and distribute to all processors
        void distributeEpsilon();

        //- Create local force field if needed
        void createForceField();

        //- Write fields on the final time step
        void writeFieldsAtEnd();

public:

    //- Runtime type information
    TypeName("actuatorModelBase");


    // Constructors

        //- Construct from components
        actuatorModelBase
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~actuatorModelBase();


    // Member Functions

        // Routines for actuatorModelBase for parallel efficiency

        //- Alloction part of the initial setup
        virtual void allocateAL();

        //- Work to be done first time for calculateALData
        virtual void initializeAL();

        //- Set up positions by calling findCells on all actuator line elements
        virtual void setupPositions(bool includeRing);

        //- Calculate the force once velocity is sampled
        virtual void calculateForces();


        // Access

            // to easier collect actuatorLines from several sources
            const List<actuatorLineSource*>& actuatorLines() const;


        // Source term addition

            virtual void addForce
            (
                volVectorField& forceField,
                scalar scale,
                bool compressible
            );

            //- Source term to turbulence scalars
            virtual void addTurbulence
            (
                fvMatrix<scalar>& eqn,
                const word fieldName
            );

            // select proper force field and call addForce
            void addForceFromChild
            (
                volVectorField& forceField,
                scalar scale,
                bool compressible
            );

            //- Add source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add source term to turbulence model equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldI
            );


        // I-O

            //- Read source dictionary
            virtual bool read(const dictionary& dict);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
