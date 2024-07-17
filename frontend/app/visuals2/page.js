'use client';

import React from 'react'
import ProteinViewer, {ProteinStage, ViewerStage} from "@jowillianto/ngl-viewer/dist"
import { 
  ComponentUIDataT 
} from "@jowillianto/ngl-viewer/dist/ngl-viewer/user-interface/component-data";

const App = () => {
  const component = {
    type : "file",
    props : {
      file : "http://localhost:8000/api/getfile/LCO_pristine.mol2",
      fileSettings : {},
      viewSettings : [{
        type : 'cartoon', params : {}
      }]
    },
    config : {}
  }
  return (
    <ProteinViewer initialComponents = {[component]}>
      <ViewerStage height = "800px" width = "800px" />
    </ProteinViewer>
  )
}

export default App;