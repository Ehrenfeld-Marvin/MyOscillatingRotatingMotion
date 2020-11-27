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

#include "oscillatingRotatingFunction1Motion.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(oscillatingRotatingFunction1Motion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        oscillatingRotatingFunction1Motion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::oscillatingRotatingFunction1Motion::
oscillatingRotatingFunction1Motion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    origin_(SBMFCoeffs_.lookup("origin")),
    amplitude_(SBMFCoeffs_.lookup("amplitude")),
    omega_(Function1<scalar>::New("omega", SBMFCoeffs_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::oscillatingRotatingFunction1Motion::
~oscillatingRotatingFunction1Motion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::oscillatingRotatingFunction1Motion::
transformation() const
{
    scalar t = time_.value();
    scalar omega = omega_->value(t);

    vector eulerAngles = amplitude_*sin(omega*t); /// WARNING: this implies that the function1 will not really be respected, as the angle will be calculated at each time step as if it had for the whole time the ang. vel. of the current timestep

    // Convert the rotational motion from deg to rad
    eulerAngles *= degToRad();

    quaternion R(quaternion::XYZ, eulerAngles);
    septernion TR(septernion(-origin_)*R*septernion(origin_));

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::oscillatingRotatingFunction1Motion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    omega_.reset
    (
        Function1<scalar>::New("omega", SBMFCoeffs_).ptr()
    );

    return true;
}


// ************************************************************************* //
