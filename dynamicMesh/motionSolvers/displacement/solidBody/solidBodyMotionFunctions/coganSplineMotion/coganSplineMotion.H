/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::solidBodyMotionFunctions::coganSplineMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Oscillating rotation.

SourceFiles
    coganSplineMotion.C

\*---------------------------------------------------------------------------*/

#ifndef coganSplineMotion_H
#define coganSplineMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "Function1.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class coganSplineMotion Declaration
\*---------------------------------------------------------------------------*/

class coganSplineMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of gravity
        point origin_;


        //- Axis vector
        const vector axis_;


        //- initial offset (only used to accomodate old offseted mesh)
        const scalar initialOffset_;

        //- Angular velocity (rad/sec)
        autoPtr<Function1<scalar>> omega_;
/*
        //- splineMultiplier
        autoPtr<Function1<scalar>> splineMultiplier_;
	*/

	//- Spline parameters:
    scalar a[8] = {-0.055182006826239, 0.035337914131215, -0.029636063996805, -0.143299538627938, 0.473373979170938, -0.26193409891157, 0.002528778311838, -0.02118896325144};
    scalar b[8] = {-0.000929954772379, -0.130949495214121, -0.047686496643261, -0.11751482734352, -0.455156410705772, 0.660204750817218, 0.043037070168281, 0.048995363693552};
    scalar c[8] = {0.29432821473092, 0.190750336921657, 0.050449957000189, -0.079298859449855, -0.529073798044311, -0.368029208313086, 0.184295626313202, 0.256577730841283};
    scalar d[8] = {-0.032963, 0.170894, 0.257053, 0.252903, 0.048708, -0.418252, -0.426954, -0.254436};

    // Private Member Functions

        //- No copy construct
        coganSplineMotion(const coganSplineMotion&) = delete;

        //- No copy assignment
        void operator=(const coganSplineMotion&) = delete;


public:

    //- Runtime type information
    TypeName("coganSplineMotion");


    // Constructors

        //- Construct from components
        coganSplineMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new coganSplineMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~coganSplineMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);

//	void wrapAngle(scalar& angle);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
