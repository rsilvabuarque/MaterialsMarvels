'use client';

import React, { useEffect, useRef } from 'react';
import * as NGL from 'ngl';

const NGLViewerComponent = ({ mol2Url, trajUrl }) => {
    const viewerRef = useRef();

    useEffect(() => {
        const stage = new NGL.Stage(viewerRef.current);

        // Load the .mol2 file
//         var stringBlob = new Blob( [ pdbData ], { type: 'text/plain'} );
// stage.loadFile( stringBlob, { ext: "pdb" } );

        // stage.loadFile(mol2Url).then((o) => {
        stage.loadFile(mol2Url).then((o) => {
            o.addRepresentation('cartoon');
            o.autoView();
        }).catch(error => {
            console.error('Error loading the .mol2 file:', error);
        });
        
        // // Load the trajectory
        // stage.loadFile(trajUrl, { ext: 'lammpstrj' }).then((o) => {
        //     o.addTrajectory();
        // }).catch(error => {
        //     console.error('Error loading the trajectory file:', error);
        // });

        return () => {
            stage.dispose();
        };
    }, [mol2Url, trajUrl]);

    return <div ref={viewerRef} style={{ width: '600px', height: '400px' }}></div>;
};

export default NGLViewerComponent;