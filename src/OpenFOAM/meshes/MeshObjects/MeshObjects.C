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

#include "MeshObjects.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mesh, class Type>
Foam::MeshObjects<Mesh, Type>::MeshObjects(const Mesh& mesh)
:
    regIOobject
    (
        IOobject
        (
            Type::typeName,
            mesh.thisDb().instance(),
            mesh.thisDb()
        )
    ),
    mesh_(mesh)
{
    if (Mesh::debug)
    {
        InfoIn("MeshObjects<Mesh, Type>::MeshObjects(const Mesh& mesh)")
            << "Creating meshObjects for type " << Type::typeName << endl;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Mesh, class Type>
const Type& Foam::MeshObjects<Mesh, Type>::New
(
    const Mesh& mesh
)
{
    if (!mesh.thisDb().objectRegistry::template foundObject<Type>(Type::typeName))
    {
        return store(new Type(mesh));
    }
    else
    {
        return mesh.thisDb().objectRegistry::template lookupObject<Type>
        (
            Type::typeName
        );
    }
}


template<class Mesh, class Type>
template<class Data1>
const Type& Foam::MeshObjects<Mesh, Type>::New
(
    const Mesh& mesh,
    const Data1& d
)
{
    if (!mesh.thisDb().objectRegistry::template foundObject<Type>(Type::typeName))
    {
        return store(new Type(mesh, d));
    }
    else
    {
        return mesh.thisDb().objectRegistry::template lookupObject<Type>
        (
            Type::typeName
        );
    }
}


template<class Mesh, class Type>
template<class Data1, class Data2>
const Type& Foam::MeshObjects<Mesh, Type>::New
(
    const Mesh& mesh,
    const Data1& d1,
    const Data2& d2
)
{
    if (!mesh.thisDb().objectRegistry::template foundObject<Type>(Type::typeName))
    {
        return store(new Type(mesh, d1, d2));
    }
    else
    {
        return mesh.thisDb().objectRegistry::template lookupObject<Type>
        (
            Type::typeName
        );
    }
}


template<class Mesh, class Type>
template<class Data1, class Data2, class Data3>
const Type& Foam::MeshObjects<Mesh, Type>::New
(
    const Mesh& mesh,
    const Data1& d1,
    const Data2& d2,
    const Data3& d3
)
{
    if (!mesh.thisDb().objectRegistry::template foundObject<Type>(Type::typeName))
    {
        return store(new Type(mesh, d1, d2, d3));
    }
    else
    {
        return mesh.thisDb().objectRegistry::template lookupObject<Type>
        (
            Type::typeName
        );
    }
}


template<class Mesh, class Type>
template<class Data1, class Data2, class Data3, class Data4>
const Type& Foam::MeshObjects<Mesh, Type>::New
(
    const Mesh& mesh,
    const Data1& d1,
    const Data2& d2,
    const Data3& d3,
    const Data4& d4
)
{
    if (!mesh.thisDb().objectRegistry::template foundObject<Type>(Type::typeName))
    {
        return store(new Type(mesh, d1, d2, d3, d4));
    }
    else
    {
        return mesh.thisDb().objectRegistry::template lookupObject<Type>
        (
            Type::typeName
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Mesh, class Type>
bool Foam::MeshObjects<Mesh, Type>::Delete(const Mesh& mesh)
{
    if (mesh.thisDb().objectRegistry::template foundObject<Type>(Type::typeName))
    {
        return mesh.thisDb().checkOut
        (
            const_cast<Type&>
            (
                mesh.thisDb().objectRegistry::template lookupObject<Type>
                (
                    Type::typeName
                )
            )
        );
    }
    else
    {
        return false;
    }
}


template<class Mesh, class Type>
Foam::MeshObjects<Mesh, Type>::~MeshObjects()
{
    release();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Mesh, class Type>
bool Foam::MeshObjects<Mesh, Type>::deleteObject() const
{
    return mesh().thisDb().checkOut
    (
        const_cast<MeshObjects<Mesh, Type>&>(*this)
    );
}



// ************************************************************************* //
