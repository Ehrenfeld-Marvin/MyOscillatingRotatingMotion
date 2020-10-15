/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::solidBodyMotionFunctions::MyoscillatingRotatingMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Oscillating rotation.

SourceFiles
    MyoscillatingRotatingMotion.C

\*---------------------------------------------------------------------------*/

#ifndef MyoscillatingRotatingMotion_H
#define MyoscillatingRotatingMotion_H

#include "solidBodyMotionFunctionTrim.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class MyoscillatingRotatingMotion Declaration
\*---------------------------------------------------------------------------*/

class MyoscillatingRotatingMotion
:
    public solidBodyMotionFunctionTrim
{
    // Private Data

        //- Centre of gravity
        point origin_;

        //- Amplitude
        mutable vector amplitude_;

        //- Radial velocity
        scalar omega_;
        
	mutable scalar Ampli = 0.1;

    // Private Member Functions

        //- No copy construct
        MyoscillatingRotatingMotion(const MyoscillatingRotatingMotion&) = delete;

        //- No copy assignment
        void operator=(const MyoscillatingRotatingMotion&) = delete;


public:

    //- Runtime type information
    TypeName("MyoscillatingRotatingMotion");


    // Constructors

        //- Construct from components
        MyoscillatingRotatingMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunctionTrim> clone() const
        {
            return autoPtr<solidBodyMotionFunctionTrim>
            (
                new MyoscillatingRotatingMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~MyoscillatingRotatingMotion() = default;


    // Member Functions

	int Amplitude() const;

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
