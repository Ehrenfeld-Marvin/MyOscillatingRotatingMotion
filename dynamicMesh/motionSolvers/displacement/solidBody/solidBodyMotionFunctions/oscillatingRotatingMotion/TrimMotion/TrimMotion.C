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

#include "TrimMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"
#include <fstream>
#include <iostream>

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(TrimMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        TrimMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::TrimMotion::
TrimMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
	read(SBMFCoeffs);
	#include "Constructor.H"
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::TrimMotion::
~TrimMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

#include "WriteFilePID.H"
#include "Amplitude_Member_Function.H"
#include "PID_BackUp_Data.H"


Foam::septernion
Foam::solidBodyMotionFunctions::TrimMotion::
transformation() const
{
    scalar t = time_.value();
    scalar omega = omega_->value(t);
    
    // Rotation around axis
//    scalar angle = omega_->integrate(0, t);


//    vector eulerAngles = amplitude_ * sign(omega) * ( sin(fabs(omega*t) + initialOffset_*pi ) - sin( initialOffset_*pi ) ); /// the second sin is to ensure that the value is null at time = 0 /// WARNING: the t inside the sine that the function1 will not really be respected, as the angle will be calculated at each time step as if it had for the whole time the ang. vel. of the current timestep

	vector eulerAngles = Amplitude(omega) * sign(omega) * ( sin(fabs(omega*t) + initialOffset_*pi ) - sin( initialOffset_*pi ) );

    // Convert the rotational motion from deg to rad
    eulerAngles *= degToRad();

    quaternion R(quaternion::XYZ, eulerAngles);
    septernion TR(septernion(-origin_)*R*septernion(origin_));

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::TrimMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

	SBMFCoeffs_.lookup("origin") >> origin_;
	SBMFCoeffs_.lookup("omega_value") >> omega_PID;
	SBMFCoeffs_.lookup("amplitude") >> amplitude_;
//	SBMFCoeffs_.lookup("omega") >> omega_;
	SBMFCoeffs_.lookup("initialOffset") >> initialOffset_;
	SBMFCoeffs_.lookup("K_P") >> K_P;
	SBMFCoeffs_.lookup("K_I") >> K_I;
	SBMFCoeffs_.lookup("K_D") >> K_D;
	SBMFCoeffs_.lookup("Target") >> Target;
	SBMFCoeffs_.lookup("amplitude_begin") >> amplitude_begin;
		
    omega_.reset
    (
        Function1<scalar>::New("omega", SBMFCoeffs_).ptr()
    );

    return true;
}


// ************************************************************************* //
