    forAll (mesh.boundary(), patchI)
    {
        if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {
            const immersedBoundaryFvPatch& ibPatch =
                refCast<const immersedBoundaryFvPatch>
                (
                    mesh.boundary()[patchI]
                );
                const labelList& ibCells = ibPatch.ibPolyPatch().ibCells();
                const scalarField& ibCellVolumes=ibPatch.ibPolyPatch().ibCellVolumes();
                const vectorField& ibFaceAreas=ibPatch.ibPolyPatch().ibPatch().faceNormals();
                //Info<<"ibCells:"<<ibCells<<endl;
                Info<<"ibCellVolumes:"<<ibCellVolumes<<endl;
                //Info<<"ibFaceAreas:"<<ibFaceAreas<<endl;


        }
    }
