/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

Class
    Foam::solidBodyMotionSolverTrim

Description
    Solid-body motion of the mesh specified by a run-time selectable
    motion function.

SourceFiles
    solidBodyMotionSolverTrim.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyMotionSolverTrim_H
#define solidBodyMotionSolverTrim_H

#include "points0MotionSolver.H"
#include "solidBodyMotionFunctionTrim.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class solidBodyMotionSolverTrim Declaration
\*---------------------------------------------------------------------------*/

class solidBodyMotionSolverTrim
:
    public points0MotionSolver
{
    // Private data

        //- The motion control function
        autoPtr<solidBodyMotionFunctionTrim> SBMFPtr_;

        //- Points to move when cell zone is supplied
        labelList pointIDs_;

        //- Flag to indicate whether all cells should move
        bool moveAllCells_;


    // Private Member Functions

        //- No copy construct
        solidBodyMotionSolverTrim(const solidBodyMotionSolverTrim&) = delete;

        //- No copy assignment
        void operator=(const solidBodyMotionSolverTrim&) = delete;


public:

    //- Runtime type information
    TypeName("solidBodyTrim");


    // Constructors

        //- Construct from mesh and dictionary
        solidBodyMotionSolverTrim
        (
            const polyMesh&,
            const IOdictionary&
        );


    //- Destructor
    ~solidBodyMotionSolverTrim();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
