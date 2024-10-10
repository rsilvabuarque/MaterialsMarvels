"use client";

import 'ketcher-react/dist/index.css';
import './page.module.css';
import { useEffect, useState } from 'react';
import Card from 'react-bootstrap/Card';
import Button from 'react-bootstrap/Button';
import LoadingButton from './LoadingBtn';
import Navigation from '../components/Navigation';
import dynamic from "next/dynamic";
// import StructureEditor from './StructureEditor';

export default function EditorPage() {
  const [isClient, setIsClient] = useState(false);

  useEffect(() => {
    // Ensure this code runs only on the client
    setIsClient(true);
  }, []);

  if (!isClient || typeof window === 'undefined') {
    return <div style={{ position: 'relative', height: '100vh', overflow: 'hidden' }}>Loading...</div>;
  }

  const StructureEditor = dynamic(() =>
    import("./StructureEditor"), {   ssr: false });

  return (
    <>
      <Navigation />
      <center>
      <Card bg={"dark"} text={"white"} style={{ width: '40rem' }}>
      <Card.Header as="h5">Visualizing the Self-Assembly of Ions</Card.Header>
      <Card.Body>
        <Card.Title>What does this tool do?</Card.Title>
        <Card.Text>
          Explanation of this tool.
          <br />
          Explanation of this tool.
        </Card.Text>
        <Card.Title>How do I use this tool?</Card.Title>
        <Card.Text>
          Instructions on how to create a structure.
          <br />
          Instructions on how to create a structure.
        </Card.Text>
        <Button variant="primary">Watch the demo video</Button>
      </Card.Body>
    </Card>
    </center>
    <br />
      <StructureEditor />
      <LoadingButton>Visualize Self-Assembly</LoadingButton>
    </>
  );
}
