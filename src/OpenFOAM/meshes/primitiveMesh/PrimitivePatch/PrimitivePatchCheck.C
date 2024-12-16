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

Description
    Checks topology of the patch.

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"
#include "Map.H"
#include "ListOps.H"
#include "OSstream.H"
#include "OFstream.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::visitPointRegion
(
    const label pointi,
    const labelList& pFaces,
    const label startFacei,
    const label startEdgeI,
    boolList& pFacesHad
) const
{
    label index = pFaces.find(startFacei);

    if (!pFacesHad[index])
    {
        // Mark face as been visited.
        pFacesHad[index] = true;

        // Step to next edge on face which is still using pointi
        const labelList& fEdges = faceEdges()[startFacei];

        label nextEdgeI = -1;

        forAll(fEdges, i)
        {
            label edgeI = fEdges[i];

            const edge& e = edges()[edgeI];

            if (edgeI != startEdgeI && (e[0] == pointi || e[1] == pointi))
            {
                nextEdgeI = edgeI;

                break;
            }
        }

        if (nextEdgeI == -1)
        {
            FatalErrorInFunction
                << "Problem: cannot find edge out of " << fEdges
                << "on face " << startFacei << " that uses point " << pointi
                << " and is not edge " << startEdgeI << abort(FatalError);
        }

        // Walk to next face(s) across edge.
        const labelList& eFaces = edgeFaces()[nextEdgeI];

        forAll(eFaces, i)
        {
            if (eFaces[i] != startFacei)
            {
                visitPointRegion
                (
                    pointi,
                    pFaces,
                    eFaces[i],
                    nextEdgeI,
                    pFacesHad
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
typename Foam::PrimitivePatch<FaceList, PointField>::surfaceTopo
Foam::PrimitivePatch<FaceList, PointField>::surfaceType() const
{
    DebugInFunction << "Calculating patch topology" << nl;

    const labelListList& edgeFcs = edgeFaces();

    surfaceTopo pType = MANIFOLD;

    forAll(edgeFcs, edgeI)
    {
        label nNbrs = edgeFcs[edgeI].size();

        if (nNbrs < 1 || nNbrs > 2)
        {
            pType = ILLEGAL;

            // Can exit now. Surface is illegal.
            return pType;
        }
        else if (nNbrs == 1)
        {
            // Surface might be open or illegal so keep looping.
            pType = OPEN;
        }
    }

    DebugInFunction << "Calculated patch topology" << nl;

    return pType;
}


template<class FaceList, class PointField>
bool
Foam::PrimitivePatch<FaceList, PointField>::checkTopology
(
    const bool report,
    labelHashSet* setPtr
) const
{
    DebugInFunction << "Checking patch topology" << nl;

    // Check edgeFaces

    const labelListList& edgeFcs = edgeFaces();

    bool illegalTopo = false;

    forAll(edgeFcs, edgeI)
    {
        label nNbrs = edgeFcs[edgeI].size();

        if (nNbrs < 1 || nNbrs > 2)
        {
            illegalTopo = true;

            if (report)
            {
                Info<< "Edge " << edgeI << " with vertices:" << edges()[edgeI]
                    << " has " << nNbrs << " face neighbours"
                    << endl;
            }

            if (setPtr)
            {
                const edge& e = edges()[edgeI];

                setPtr->insert(meshPoints()[e.start()]);
                setPtr->insert(meshPoints()[e.end()]);
            }
        }
    }

    DebugInFunction << "Checked patch topology" << nl;

    return illegalTopo;
}


template<class FaceList, class PointField>
bool
Foam::PrimitivePatch<FaceList, PointField>::checkPointManifold
(
    const bool report,
    labelHashSet* setPtr
) const
{
    const labelListList& pf = pointFaces();
    const labelListList& pe = pointEdges();
    const labelListList& ef = edgeFaces();
    const labelList& mp = meshPoints();

    bool foundError = false;

    forAll(pf, pointi)
    {
        const labelList& pFaces = pf[pointi];

        // Visited faces (as indices into pFaces)
        boolList pFacesHad(pFaces.size(), false);

        // Starting edge
        const labelList& pEdges = pe[pointi];
        label startEdgeI = pEdges[0];

        const labelList& eFaces = ef[startEdgeI];

        forAll(eFaces, i)
        {
            // Visit all faces using pointi, starting from eFaces[i] and
            // startEdgeI. Mark off all faces visited in pFacesHad.
            this->visitPointRegion
            (
                pointi,
                pFaces,
                eFaces[i],  // starting face for walk
                startEdgeI, // starting edge for walk
                pFacesHad
            );
        }

        // After this all faces using pointi should have been visited and
        // marked off in pFacesHad.

        label unset = pFacesHad.find(false);

        if (unset != -1)
        {
            foundError = true;

            label meshPointi = mp[pointi];

            if (setPtr)
            {
                setPtr->insert(meshPointi);
            }

            if (report)
            {
                Info<< "Point " << meshPointi
                    << " uses faces which are not connected through an edge"
                    << nl
                    << "This means that the surface formed by this patched"
                    << " is multiply connected at this point" << nl
                    << "Connected (patch) faces:" << nl;

                forAll(pFacesHad, i)
                {
                    if (pFacesHad[i])
                    {
                        Info<< "    " << pFaces[i] << endl;
                    }
                }

                Info<< nl << "Unconnected (patch) faces:" << nl;
                forAll(pFacesHad, i)
                {
                    if (!pFacesHad[i])
                    {
                        Info<< "    " << pFaces[i] << endl;
                    }
                }
            }
        }
    }

    return foundError;
}

//---------------------------------------------------------------------------sajad

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::
writeVTK
(
    const fileName& name,
    const FaceListType& faces,
    const Field<point_type>& points
)
{
    // Write patch and points into VTK
    OFstream mps(name + ".vtk");

    mps << "# vtk DataFile Version 2.0" << nl
        << name << ".vtk" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << points.size() << " float" << nl;

    // Write points
    List<float> mlpBuffer(3*points.size());

    label counter = 0;
    forAll (points, i)
    {
        mlpBuffer[counter++] = float(points[i].x());
        mlpBuffer[counter++] = float(points[i].y());
        mlpBuffer[counter++] = float(points[i].z());
    }

    forAll (mlpBuffer, i)
    {
        mps << mlpBuffer[i] << ' ';

        if (i > 0 && (i % 10) == 0)
        {
            mps << nl;
        }
    }

    // Write faces
    label nFaceVerts = 0;

    forAll (faces, faceI)
    {
        nFaceVerts += faces[faceI].size() + 1;
    }
    labelList mlfBuffer(nFaceVerts);

    counter = 0;
    forAll (faces, faceI)
    {
     //   const face& f = faces[faceI];sajad

        mlfBuffer[counter++] = faces[faceI].size();

        forAll (faces[faceI], fpI)
        {
            mlfBuffer[counter++] = faces[faceI][fpI];
        }
    }
    mps << nl;

    mps << "POLYGONS " << faces.size() << ' ' << nFaceVerts << endl;

    forAll (mlfBuffer, i)
    {
        mps << mlfBuffer[i] << ' ';

        if (i > 0 && (i % 10) == 0)
        {
            mps << nl;
        }
    }
    mps << nl;
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::
writeVTKNormals
(
    const fileName& name,
    const FaceListType& faces,
    const Field<point_type>& points
)
{
    // Write patch and points into VTK
    OFstream mps(name + ".vtk");

    mps << "# vtk DataFile Version 2.0" << nl
        << name << ".vtk" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << faces.size() << " float" << nl;

    // Write points
    List<float> mlPointBuffer(3*faces.size());

    label counter = 0;
    forAll (faces, i)
    {
        const vector c = faces[i].centre(points);

        mlPointBuffer[counter++] = float(c.x());
        mlPointBuffer[counter++] = float(c.y());
        mlPointBuffer[counter++] = float(c.z());
    }

    forAll (mlPointBuffer, i)
    {
        mps << mlPointBuffer[i] << ' ';

        if (i > 0 && (i % 10) == 0)
        {
            mps << nl;
        }
    }
    mps << nl;

    // Write normals
    mps << "POINT_DATA " << faces.size() << nl
        << "FIELD attributes " << 1 << nl
        << "normals" << " 3 "
        << faces.size() << " float" << nl;

    List<float> mlNormalBuffer(3*faces.size());

    counter = 0;
    forAll (faces, i)
    {
        const vector n = faces[i].normal(points);

        mlNormalBuffer[counter++] = float(n.x());
        mlNormalBuffer[counter++] = float(n.y());
        mlNormalBuffer[counter++] = float(n.z());
    }

    forAll (mlNormalBuffer, i)
    {
        mps << mlNormalBuffer[i] << ' ';

        if (i > 0 && (i % 10) == 0)
        {
            mps << nl;
        }
    }
    mps << nl;
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::
writeVTK
(
    const fileName& name
) const
{
    PrimitivePatch<FaceList, PointField>::writeVTK
    (
        name,
        this->localFaces(),
        this->localPoints()
    );
}

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::
writeVTKNormals
(
    const fileName& name
) const
{
    PrimitivePatch<FaceList, PointField>::writeVTKNormals
    (
        name,
        this->localFaces(),
        this->localPoints()
    );
}

// ************************************************************************* //
