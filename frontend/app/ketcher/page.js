"use client";

import 'ketcher-react/dist/index.css';
import './page.module.css';
import dynamic from 'next/dynamic';
import { useEffect, useState } from 'react';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import LoadingButton from './LoadingBtn';


export default function EditorPage() {
  const [isClient, setIsClient] = useState(false);

  useEffect(() => {
    // Ensure this code runs only on the client
    setIsClient(true);
  }, []);

  if (!isClient || typeof window === 'undefined') {
    return <div style={{ position: 'relative', height: '100vh', overflow: 'hidden' }}>Loading...</div>;
  }

  // Dynamically import the Editor component from ketcher-react
  const Editor = dynamic(() => import('ketcher-react').then(mod => mod.Editor), {
    ssr: false,
  });

  const structServiceProvider = new StandaloneStructServiceProvider();

  return (
    <>
      <div style={{ position: 'relative', height: '80vh', overflow: 'hidden' }}>
        <Editor
          staticResourcesUrl={process.env.PUBLIC_URL}
          structServiceProvider={structServiceProvider}
          style={{ width: '100%', height: '100%' }}
          onInit={(ketcher) => {
            // If window defined, attach the ketcher instance to it
            if (typeof window !== 'undefined'){
              window.ketcher = ketcher;
            } else {
              // Try again later
              setTimeout(() => {
                window.ketcher = ketcher;
              }, 1000);
            }
          }}
        />
      </div>
      <LoadingButton>Visualize Structure</LoadingButton>
    </>
  );
}
