'use client'

import NGLViewerComponent from '../visuals/NGLViewerComponent';

const VisualsPage = () => {
    return (
        <div>
            <h1>Self-Assembly Visualizer</h1>
            <NGLViewerComponent 
              mol2Url="/api/getfile/LCO_pristine.mol2" 
              trajUrl="/api/getfile/LCO_mlff.heat.lammpstrj"
          />
        </div>
    );
  }
  
  export default VisualsPage;