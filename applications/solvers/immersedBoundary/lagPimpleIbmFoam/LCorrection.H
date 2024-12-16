//typedef typename Foam::KinematicParcel<Foam::Particle> parcelType ;

forAllIter( basicKinematicCloud , kinematicCloud, iter)
{
   // parcelType& p = iter(); //paticle
    if(gamma[iter().cell()]==0)
        {
                Info<<"MassCell:"<<iter().mass()<<endl;
                ParticleMassInIbPatchCell[iter().cell()]+=iter().mass();
                particleCounter++;
		kinematicCloud.deleteParticle(iter());
	}
}
//-----------------------------------------------------------------------------------
    // forAll(mesh.C(),II)
    // {
		//~ Info<<II<<"\t"<<mesh.C()[II]<<endl;
    // }

//-----------------------------------------------------------------------------------
   // forAll (mesh.boundary(), patchI)
   // {
    //    if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
   //     {
     //       const immersedBoundaryFvPatch& ibPatch =
     //           refCast<const immersedBoundaryFvPatch>
     //           (
      //              mesh.boundary()[patchI]
     //           );
     //           const labelList& ibCells = ibPatch.ibPolyPatch().ibCells();
      //         // Info<<"ibCells:"<<ibCells<<endl;
     //   }
   // }

