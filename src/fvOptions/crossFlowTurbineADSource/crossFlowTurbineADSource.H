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
    Foam::fv::crossFlowTurbineADSource

Description
    Cell based momentum source that represents a cross-flow turbine

SourceFiles
    crossFlowTurbineADSource.C

\*---------------------------------------------------------------------------*/

#ifndef crossFlowTurbineADSource_H
#define crossFlowTurbineADSource_H

#include "crossFlowTurbineALSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class crossFlowTurbineALSource Declaration
\*---------------------------------------------------------------------------*/

class crossFlowTurbineADSource
:
    public crossFlowTurbineALSource
{

protected:

    // Protected data

        //- number of divisions per revolution
        label divisions_;
        
        //- Run the force calculation loop additional times to allow
        // dynamic stall models with hysteresis develop the solution
        // 1 additional loop should be sufficient
        label dynStallLoop_;
        
        //- Multiply the force with this value (basically allows to set up
        // one bladed turbines that act like 3 bladed turbines)
        scalar bladeMultiplier_;

        //- internal time used to perform a loop within the openfoam-step
        scalar customTime_;
        
        
    // Protected Member Functions

        //- Create the co-ordinate system
        // virtual void createCoordinateSystem();

        //- Use turbine base class rotate method
        //using turbineALSource::rotate;
        void rotateAD(bool updateOnly = false);

        //- Set up a list of neighbor cells for faster force projection
        void buildInfluenceCells();

        //- Rotate the turbine a specified angle about its axis
        virtual void rotate(scalar radians);


public:

    //- Runtime type information
    TypeName("crossFlowTurbineADSource");


    // Constructors

        //- Construct from components
        crossFlowTurbineADSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~crossFlowTurbineADSource();


    // Member Functions

        // Source term addition

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

            //- Print dictionary values
            //virtual void printCoeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
