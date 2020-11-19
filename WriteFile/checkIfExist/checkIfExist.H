/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
#ifndef functionObjects_checkIfExist_H
#define functionObjects_checkIfExist_H

#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace writeFile
{

/*---------------------------------------------------------------------------*\
                 Class functionObjects::writeFile Declaration
\*---------------------------------------------------------------------------*/

class checkIfExist
{
protected:

    // Protected data

        //- Reference to the region objectRegistry
        const objectRegistry& fileObr_;

        //- Prefix
        const fileName prefix_;

        //- Name of file
        word fileName_;

        //- File pointer
        autoPtr<OFstream> filePtr_;

        //- Write precision
        label writePrecision_;

        //- Flag to enable/disable writing to file
        bool writeToFile_;

        //- Flag to update the header, e.g. on mesh changes.
        //- Default is true.
        bool updateHeader_;

        //- Flag to identify whether the header has been written
        bool writtenHeader_;

        //- Flag to use the specified user time, e.g. CA deg instead
        //- of seconds.  Default = true
        bool useUserTime_;

        //- Start time value
        scalar startTime_;


    // Protected Member Functions

        virtual autoPtr<OFstream> createFile
        (
            const word& name,
            scalar timeValue
        ) const;

        //- Return autoPtr to a new file using the simulation start time
        virtual autoPtr<OFstream> createFile
        (
            const word& name
        ) const;



public:



    // Constructors

        //- Construct from objectRegistry, prefix, fileName
        writeFile
        (
            const objectRegistry& obr,
            const fileName& prefix,
            const word& name = "undefined",
            const bool writeToFile = true
        );

        //- Construct from objectRegistry, prefix, fileName
        //- and read options from dictionary
        writeFile
        (
            const objectRegistry& obr,
            const fileName& prefix,
            const word& name,
            const dictionary& dict,
            const bool writeToFile = true
        );

        //- Construct copy
        writeFile(const writeFile& wf);


    //- Destructor
    virtual ~writeFile() = default;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
