// "use client";

// import Head from 'next/head';

// export default function VisualsPage() {
//   return (
//   // <div>
//   // <Head>
//   //   <script type="text/javascript" src={"/static/script.js"}></script>
//   //   <script src="https://unpkg.com/ngl"></script>
//   // </Head>
//   // <div id="viewport" style={{"width":"400px", "height":"300px"}}></div>
//   // </div>
//   <p>Hello World</p>
//   );
// }

import NGLViewerComponent from './NGLViewerComponent';

const HomePage = () => {
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

export default HomePage;