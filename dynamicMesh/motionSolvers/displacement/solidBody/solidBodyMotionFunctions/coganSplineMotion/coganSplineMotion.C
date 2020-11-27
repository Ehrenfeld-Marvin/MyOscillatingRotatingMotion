/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "coganSplineMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(coganSplineMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        coganSplineMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::coganSplineMotion::
coganSplineMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    origin_(SBMFCoeffs_.get<vector>("origin")),
    axis_(SBMFCoeffs_.get<vector>("axis")),
    initialOffset_(SBMFCoeffs_.get<scalar>("initialOffset")),
    omega_(Function1<scalar>::New("omega", SBMFCoeffs_))
//    splineMultiplier_(Function1<scalar>::New("splineMultiplier", SBMFCoeffs_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::coganSplineMotion::
~coganSplineMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::coganSplineMotion::
transformation() const
{
    scalar t = time_.value();
    scalar omega = omega_->value(t);
    // Rotation around axis
//    scalar angle = omega_->integrate(0, t);


// OLD:    vector eulerAngles = amplitude_ * sign(omega) * ( sin(fabs(omega*t) + initialOffset_*pi ) - sin( initialOffset_*pi ) ); /// the second sin is to ensure that the value is null at time = 0 /// WARNING: the t inside the sine that the function1 will not really be respected, as the angle will be calculated at each time step as if it had for the whole time the ang. vel. of the current timestep


    // Convert the rotational motion from deg to rad
// OLD:    eulerAngles *= degToRad();

//    int Nblade = 0; // temp for testing
//    scalar offset = - Nblade / 6 * 2*pi; // taking first blade as 0 and increamenting clockwise
    scalar Psi = omega*t + initialOffset_*pi; // offset for cycloidal rotor blade position
//    scalar twoPi = 2.0 * pi;
    Psi = Psi - 2.0*pi * floor( Psi / (2.0*pi) );
//    scalar Psi2 = 3.12;
//    wrapAngle(Psi2);

    int i = floor(Psi / (pi/4)); // in which part of the spline are we?
    

    scalar angle = fmax(-t*20.,-1.0) * ( // using fmax to make everything negative (Cogan theta is computed in opposite direction)
	    a[i] * pow( (Psi - pi/4*i) , 3 )
	    + b[i] * pow( (Psi - pi/4*i) , 2 )
	    + c[i] * (Psi - pi/4*i)
	    + d[i] 
	    - 25*degToRad()*sin( initialOffset_*pi ) // correction for initially offsetted mesh
	    );
    		

    quaternion R(axis_, angle);
    septernion TR(septernion(-origin_)*R*septernion(origin_));

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::coganSplineMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

//    SBMFCoeffs_.lookup("omega") >> omega_;
     omega_.reset
     (
         Function1<scalar>::New("omega", SBMFCoeffs_)
     );
     /*
     splineMultiplier_.reset
     (
         Function1<scalar>::New("splineMultiplier", SBMFCoeffs_)
     );
*/
    return true;
}
/*
void Foam::solidBodyMotionFunctions::coganSplineMotion::wrapAngle
(
    scalar& angle
)
{
    scalar twoPi = 2.0 * 3.141592865358979;
    angle = angle - twoPi * floor( angle / twoPi );
}
*/


// ************************************************************************* //
