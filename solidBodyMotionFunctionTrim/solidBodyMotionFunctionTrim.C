/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "solidBodyMotionFunctionTrim.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyMotionFunctionTrim, 0);
    defineRunTimeSelectionTable(solidBodyMotionFunctionTrim, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctionTrim::solidBodyMotionFunctionTrim
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    SBMFCoeffs_
    (
        SBMFCoeffs.optionalSubDict
        (
            SBMFCoeffs.get<word>("solidBodyMotionFunctionTrim") + "Coeffs"
        )
    ),
    time_(runTime)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctionTrim::~solidBodyMotionFunctionTrim()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::solidBodyMotionFunctionTrim::read(const dictionary& SBMFCoeffs)
{
    SBMFCoeffs_ = SBMFCoeffs.optionalSubDict(type() + "Coeffs");

    return true;
}


void Foam::solidBodyMotionFunctionTrim::writeData(Ostream& os) const
{
    os << SBMFCoeffs_;
}


// ************************************************************************* //
