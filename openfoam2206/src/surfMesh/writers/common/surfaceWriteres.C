/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "surfaceWriteres.H"

#include "MeshedSurfaceProxy.H"
#include "proxySurfaceWriter.H"

#include "HashTable.H"
#include "word.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceWriteres, 0);
    defineRunTimeSelectionTable(surfaceWriteres, word);
    defineRunTimeSelectionTable(surfaceWriteres, wordDict);
    addNamedToRunTimeSelectionTable
    (
        surfaceWriteres,
        surfaceWriteres,
        word,
        null
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceWriteres>
Foam::surfaceWriteres::New(const word& writeType)
{
    // Constructors without dictionary options
    auto* ctorPtr = wordConstructorTable(writeType);

    if (!ctorPtr)
    {
        if (MeshedSurfaceProxy<face>::canWriteType(writeType))
        {
            // Generally unknown, but handle via 'proxy' handler
            return autoPtr<surfaceWriteres>
            (
                new surfaceWriters::proxyWriter(writeType)
            );
        }

        FatalErrorInFunction
            << "Unknown write type \"" << writeType << "\"\n\n"
            << "Valid write types : "
            << flatOutput(wordConstructorTablePtr_->sortedToc()) << nl
            << "Valid proxy types : "
            << MeshedSurfaceProxy<face>::writeTypes() << endl
            << exit(FatalError);
    }

    return autoPtr<surfaceWriteres>(ctorPtr());
}


Foam::autoPtr<Foam::surfaceWriters>
Foam::surfaceWriters::New(const word& writeType, const dictionary& optDict)
{
    // Constructors with dictionary options
    {
        auto* ctorPtr = wordDictConstructorTable(writeType);

        if (ctorPtr)
        {
            return autoPtr<surfaceWriteres>(ctorPtr(writeOpts));
        }
    }


    // Constructors without dictionary options
    auto* ctorPtr = wordConstructorTable(writeType);

    if (!ctorPtr)
    {
        if (MeshedSurfaceProxy<face>::canWriteType(writeType))
        {
            // Generally unknown, but handle via 'proxy' handler
            return autoPtr<surfaceWriteres>
            (
                new surfaceWriters::proxyWriter(writeType, writeOpts)
            );
        }

        FatalErrorInFunction
            << "Unknown write type \"" << writeType << "\"\n\n"
            << "Valid write types : "
            << wordConstructorTablePtr_->sortedToc() << nl
            << "Valid proxy types : "
            << MeshedSurfaceProxy<face>::writeTypes() << endl
            << exit(FatalError);
    }

    return autoPtr<surfaceWriteres>(ctorPtr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriteres::surfaceWriteres()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceWriteres::~surfaceWriteres()
{}


// ************************************************************************* //
