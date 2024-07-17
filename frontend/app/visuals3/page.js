'use client'

import Head from 'next/head';
import Script from 'next/script';
import { useEffect, useRef } from 'react';

const NGLView = () => {
  const viewportRef = useRef(null);

  useEffect(() => {
    const initializeNGL = () => {
      if (window.NGL && viewportRef.current) {
        const stage = new window.NGL.Stage(viewportRef.current);
        stage.loadFile('rcsb://1crn.mmtf', { defaultRepresentation: true });

        return () => {
          stage.dispose(); // Cleanup function to dispose of resources
        };
      }
    };

    if (document.readyState === "complete") {
      initializeNGL();
    } else {
      window.addEventListener("load", initializeNGL);
      return () => window.removeEventListener("load", initializeNGL);
    }
  }, []);

  return (
    <div>
      <Head>
        <title>NGL - embedded</title>
        <meta charSet="utf-8" />
        <meta name="viewport" content="width=device-width, initial-scale=1, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0" />
        <Script src="ngl.js" strategy="lazyOnload" />
      </Head>
      <div id="viewport" ref={viewportRef} style={{ width: '800px', height: '800px' }}></div>
    </div>
  );
};

export default NGLView;



// import Head from 'next/head';
// import { useEffect, useRef } from 'react';

// const NGLView = () => {
//   const viewportRef = useRef(null);

//   useEffect(() => {
//     const NGL = require('ngl');  // Import NGL if available as a module
//     let stage;

//     if (viewportRef.current) {
//       stage = new NGL.Stage(viewportRef.current);
//       stage.loadFile('rcsb://1crn.mmtf', { defaultRepresentation: true });
//     }

//     return () => {
//       // Cleanup function to dispose of resources when component unmounts
//       if (stage) {
//         stage.dispose();
//       }
//     };
//   }, []);

//   return (
//     <div>
//       <Head>
//         <title>NGL - embedded</title>
//         <meta charSet="utf-8" />
//         <meta name="viewport" content="width=device-width, initial-scale=1, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0" />
//       </Head>
//       <div id="viewport" ref={viewportRef} style={{ width: '800px', height: '800px' }}></div>
//     </div>
//   );
// };

// export default NGLView;