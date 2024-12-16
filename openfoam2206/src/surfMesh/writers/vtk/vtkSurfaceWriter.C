/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "vtkSurfaceWriter.H"
#include "foamVtkSurfaceWriter.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceWriter.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(vtkWriter);
    addToRunTimeSelectionTable(surfaceWriter, vtkWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, vtkWriter, wordDict);

    // Accept vtp ending as well
    addNamedToRunTimeSelectionTable
    (
        surfaceWriter,
        vtkWriter,
        word,
        vtp
    );
    addNamedToRunTimeSelectionTable
    (
        surfaceWriter,
        vtkWriter,
        wordDict,
        vtp
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::vtkWriter::vtkWriter()
:
    surfaceWriter(),
    fmtType_(static_cast<unsigned>(vtk::formatType::INLINE_BASE64)),
    precision_(IOstream::defaultPrecision()),
    writeNormal_(false),
    writer_(nullptr)
{}


Foam::surfaceWriters::vtkWriter::vtkWriter
(
    const vtk::outputOptions& opts
)
:
    surfaceWriter(),
    fmtType_(static_cast<unsigned>(opts.fmt())),
    precision_(opts.precision()),
    writeNormal_(false),
    writer_(nullptr)
{}


Foam::surfaceWriters::vtkWriter::vtkWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    fmtType_(static_cast<unsigned>(vtk::formatType::INLINE_BASE64)),
    precision_
    (
        options.getOrDefault("precision", IOstream::defaultPrecision())
    ),
    writeNormal_(options.getOrDefault("normal", false)),
    writer_(nullptr)
{
    // format: ascii | binary
    // legacy: true | false

    vtk::outputOptions opts(vtk::formatType::INLINE_BASE64);

    opts.ascii
    (
        IOstreamOption::ASCII
     == IOstreamOption::formatEnum("format", options, IOstreamOption::BINARY)
    );

    opts.legacy(options.getOrDefault("legacy", false));

    // Convert back to raw data type
    fmtType_ = static_cast<unsigned>(opts.fmt());
}


Foam::surfaceWriters::vtkWriter::vtkWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    vtkWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::vtkWriter::vtkWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    vtkWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceWriters::vtkWriter::~vtkWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceWriters::vtkWriter::close()
{
    writer_.clear();
    surfaceWriter::close();
}


void Foam::surfaceWriters::vtkWriter::beginTime(const Time& t)
{
    writer_.clear();
    surfaceWriter::beginTime(t);
}


void Foam::surfaceWriters::vtkWriter::beginTime(const instant& inst)
{
    writer_.clear();
    surfaceWriter::beginTime(inst);
}


void Foam::surfaceWriters::vtkWriter::endTime()
{
    writer_.clear();
    surfaceWriter::endTime();
}


Foam::fileName Foam::surfaceWriters::vtkWriter::write()
{
    checkOpen();

    if (needsUpdate())
    {
        writer_.clear();
    }
    merge();

    // From raw unsigned values to vtk::outputOptions
    vtk::outputOptions opts(static_cast<vtk::formatType>(fmtType_), precision_);


    // Geometry:  rootdir/<TIME>/surfaceName.{vtk|vtp}

    fileName outputFile = outputPath_;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile = outputPath_.path() / timeName() / outputPath_.name();
    }
    outputFile.ext(vtk::surfaceWriter::ext(opts));

    if (verbose_)
    {
        Info<< "Writing geometry to " << outputFile << endl;
    }

    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (!writer_ && (Pstream::master() || !parallel_))
    {
        writer_.reset
        (
            new vtk::surfaceWriter
            (
                surf.points(),
                surf.faces(),
                opts,
                outputFile,
                false  // serial!
            )
        );

        if (this->hasTime())
        {
            // Time name in title
            writer_->setTime(currTime_);
            writer_->writeTimeValue();
        }
        else
        {
            // Surface name in title
            writer_->beginFile(outputPath_.nameLessExt());
        }

        writer_->writeGeometry();

        if (writeNormal_)
        {
            const faceList& fcs = surf.faces();
            const pointField& pts = surf.points();

            Field<vector> normals(fcs.size());
            forAll(fcs, facei)
            {
                normals[facei] = fcs[facei].areaNormal(pts);
            }

            label nCellData = 1;

            if (!this->isPointData())
            {
                // Ill-defined with legacy() if nFields_ not properly set...
                nCellData += nFields_;
            }

            writer_->beginCellData(nCellData);
            writer_->write("area-normal", normals);
        }
    }

    wroteGeom_ = true;
    return outputFile;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::surfaceWriters::vtkWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    // Field:  rootdir/<TIME>/surfaceName.{vtk|vtp}

    // Open file, writing geometry (if required)
    fileName outputFile = this->write();

    // Implicit geometry merge()
    tmp<Field<Type>> tfield = adjustField(fieldName, mergeField(localValues));

    if (verbose_)
    {
        Info<< " to " << outputFile << endl;
    }


    if (Pstream::master() || !parallel_)
    {
        if (!nFields_ && writer_->legacy())
        {
            // Emit error message, but attempt to recover anyhow
            nFields_ = 1;

            FatalErrorInFunction
                << "Using VTK legacy format, but did not define nFields!"
                << nl
                << "Assuming nFields=1 (may be incorrect) and continuing..."
                << nl
                << "    Field " << fieldName << " to " << outputFile << nl;

            Info<< FatalError;
            Info<< endl;
        }

        if (this->isPointData())
        {
            writer_->beginPointData(nFields_);
        }
        else
        {
            writer_->beginCellData(nFields_);
        }

        writer_->write(fieldName, tfield());
    }

    wroteGeom_ = true;
    return outputFile;
}

    void Foam::surfaceWriters::vtkWriter::write
    (                                                                         
        const fileName& outputDir,                                            
        const fileName& surfaceName,                                          
        const pointField& points,                                             
        const faceList& faces,                                                
        const word& fieldName,                                                
        const Field<scalar>& values,                                       
        const bool isNodeValues,                                              
        const bool verbose                                                    
    ) const                                                                   
    {                                                                         
        writeTemplate                                                        
        (                                                                     
            outputDir,                                                        
            surfaceName,                                                      
            points,                                                           
            faces,                                                            
            fieldName,                                                        
            values,                                                           
            isNodeValues,                                                     
            verbose                                                           
        );                                                                    
    }
    void Foam::surfaceWriters::vtkWriter::write
    (                                                                         
        const fileName& outputDir,                                            
        const fileName& surfaceName,                                          
        const pointField& points,                                             
        const faceList& faces,                                                
        const word& fieldName,                                                
        const Field<vector>& values,                                       
        const bool isNodeValues,                                              
        const bool verbose                                                    
    ) const                                                                   
    {                                                                         
        writeTemplate                                                        
        (                                                                     
            outputDir,                                                        
            surfaceName,                                                      
            points,                                                           
            faces,                                                            
            fieldName,                                                        
            values,                                                           
            isNodeValues,                                                     
            verbose                                                           
        );                                                                    
    }
    void Foam::surfaceWriters::vtkWriter::write
    (                                                                         
        const fileName& outputDir,                                            
        const fileName& surfaceName,                                          
        const pointField& points,                                             
        const faceList& faces,                                                
        const word& fieldName,                                                
        const Field<sphericalTensor>& values,                                       
        const bool isNodeValues,                                              
        const bool verbose                                                    
    ) const                                                                   
    {                                                                         
        writeTemplate                                                        
        (                                                                     
            outputDir,                                                        
            surfaceName,                                                      
            points,                                                           
            faces,                                                            
            fieldName,                                                        
            values,                                                           
            isNodeValues,                                                     
            verbose                                                           
        );                                                                    
    }
    void Foam::surfaceWriters::vtkWriter::write
    (                                                                         
        const fileName& outputDir,                                            
        const fileName& surfaceName,                                          
        const pointField& points,                                             
        const faceList& faces,                                                
        const word& fieldName,                                                
        const Field<symmTensor>& values,                                       
        const bool isNodeValues,                                              
        const bool verbose                                                    
    ) const                                                                   
    {                                                                         
        writeTemplate                                                        
        (                                                                     
            outputDir,                                                        
            surfaceName,                                                      
            points,                                                           
            faces,                                                            
            fieldName,                                                        
            values,                                                           
            isNodeValues,                                                     
            verbose                                                           
        );                                                                    
    }
    void Foam::surfaceWriters::vtkWriter::write
    (                                                                         
        const fileName& outputDir,                                            
        const fileName& surfaceName,                                          
        const pointField& points,                                             
        const faceList& faces,                                                
        const word& fieldName,                                                
        const Field<diagTensor>& values,                                       
        const bool isNodeValues,                                              
        const bool verbose                                                    
    ) const                                                                   
    {                                                                         
        writeTemplate                                                        
        (                                                                     
            outputDir,                                                        
            surfaceName,                                                      
            points,                                                           
            faces,                                                            
            fieldName,                                                        
            values,                                                           
            isNodeValues,                                                     
            verbose                                                           
        );                                                                    
    }
    void Foam::surfaceWriters::vtkWriter::write
    (                                                                         
        const fileName& outputDir,                                            
        const fileName& surfaceName,                                          
        const pointField& points,                                             
        const faceList& faces,                                                
        const word& fieldName,                                                
        const Field<tensor>& values,                                       
        const bool isNodeValues,                                              
        const bool verbose                                                    
    ) const                                                                   
    {                                                                         
        writeTemplate                                                        
        (                                                                     
            outputDir,                                                        
            surfaceName,                                                      
            points,                                                           
            faces,                                                            
            fieldName,                                                        
            values,                                                           
            isNodeValues,                                                     
            verbose                                                           
        );                                                                    
    }

template<class Type>
void Foam::surfaceWriters::vtkWriter::writeTemplate
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os(outputDir/fieldName + '_' + surfaceName + ".vtk");

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << os.name() << endl;
    }

    writeGeometry(os, points, faces);

    // start writing data
    if (isNodeValues)
    {
        os  << "POINT_DATA ";
    }
    else
    {
        os  << "CELL_DATA ";
    }

    os  << values.size() << nl
        << "FIELD attributes 1" << nl
        << fieldName << " ";

    // Write data
    writeData(os, values);
}

void Foam::surfaceWriters::vtkWriter::writeGeometry
(
    Ostream& os,
    const pointField& points,
    const faceList& faces
)
{
    // header
    os
        << "# vtk DataFile Version 2.0" << nl
        << "sampleSurface" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl;

    // Write vertex coords
    os  << "POINTS " << points.size() << " double" << nl;
    forAll(points, pointI)
    {
        const point& pt = points[pointI];
        os  << float(pt.x()) << ' '
            << float(pt.y()) << ' '
            << float(pt.z()) << nl;
    }
    os  << nl;


    // Write faces
    label nNodes = 0;
    forAll(faces, faceI)
    {
        nNodes += faces[faceI].size();
    }

    os  << "POLYGONS " << faces.size() << ' '
        << faces.size() + nNodes << nl;

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        os  << f.size();
        forAll(f, fp)
        {
            os  << ' ' << f[fp];
        }
        os  << nl;
    }
}


namespace Foam
{

    template<>
    void Foam::surfaceWriters::vtkWriter::writeData
    (
        Ostream& os,
        const Field<scalar>& values
    )
    {
        os  << "1 " << values.size() << " double" << nl;

        forAll(values, elemI)
        {
            if (elemI)
            {
                if (elemI % 10)
                {
                    os  << ' ';
                }
                else
                {
                    os  << nl;
                }
            }

            os  << float(values[elemI]);
        }
        os  << nl;
    }


    template<>
    void Foam::surfaceWriters::vtkWriter::writeData
    (
        Ostream& os,
        const Field<vector>& values
    )
    {
        os  << "3 " << values.size() << " double" << nl;

        forAll(values, elemI)
        {
            const vector& v = values[elemI];
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << nl;
        }
    }


    template<>
    void Foam::surfaceWriters::vtkWriter::writeData
    (
        Ostream& os,
        const Field<sphericalTensor>& values
    )
    {
        os  << "1 " << values.size() << " double" << nl;

        forAll(values, elemI)
        {
            const sphericalTensor& v = values[elemI];
            os  << float(v[0]) << nl;
        }
    }


    template<>
    void Foam::surfaceWriters::vtkWriter::writeData
    (
        Ostream& os,
        const Field<symmTensor>& values
    )
    {
        os  << "6 " << values.size() << " double" << nl;

        forAll(values, elemI)
        {
            const symmTensor& v = values[elemI];
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << ' '
                << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
                << nl;

        }
    }


    template<>
    void Foam::surfaceWriters::vtkWriter::writeData
    (
        Ostream& os,
        const Field<tensor>& values
    )
    {
        os  << "9 " << values.size() << " double" << nl;

        forAll(values, elemI)
        {
            const tensor& v = values[elemI];
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << ' '
                << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
                << ' '
                << float(v[6]) << ' ' << float(v[7]) << ' ' << float(v[8])
                << nl;
        }
    }

}
// Write generic field in vtk format
template<class Type>
void Foam::surfaceWriters::vtkWriter::writeData
(
    Ostream& os,
    const Field<Type>& values
)
{
    os  << "1 " << values.size() << " double" << nl;

    forAll(values, elemI)
    {
        os  << float(0) << nl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::vtkWriter);


// ************************************************************************* //
