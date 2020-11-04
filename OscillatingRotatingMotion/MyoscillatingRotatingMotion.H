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

#include "solidBodyMotionFunction.H"
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
    public solidBodyMotionFunction
{
    // Private Data

        //- Centre of gravity
        point origin_;

        //- Amplitude
        mutable vector amplitude_;

        //- Radial velocity
        scalar omega_;
        
        int NuOfOsc;
        
	mutable float changeTime = 0;		//Intervall when the amplitude is supposed to change


	mutable double Target;			//Zielwert der Force
	mutable float K_P;			//PID Koeffizienten
	mutable float K_I;
	mutable float K_D;
	mutable double IST;
	mutable double IST_alt=0;
	mutable double IST_SUMME=0;
	mutable double e_t;
	mutable double y_t;


	
	mutable double Average_X=0;
	mutable double Average_Y=0;
	mutable double Average_Z=0;

//	mutable int Sinus_Iterator=0;
	

	mutable int I_iterator=5;
	mutable int I_intervall=5;
	mutable vector I_vector;
        
        
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
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
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

	vector Amplitude() const;

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
