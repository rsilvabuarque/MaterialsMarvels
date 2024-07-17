'use client'

import React, { useEffect, useState } from 'react';
import styles from './page.module.css';

const HomePage = () => {
  const [isLoading, setLoading] = useState(false);

  useEffect(() => {
    if (isLoading) {
      plotEnergyVsStep('/Users/lauracharria/Downloads/test_ice_melting3.log');
      setLoading(false);
    }
  }, [isLoading]);

  const handleClick = () => setLoading(true);

  return (
    <div className={styles.container}>
      <h1 className={styles.header}>Self-Assembly Visualizer</h1>
      <div className={styles.visualizerSection}>
        <button 
          className={styles.button} 
          disabled={isLoading} 
          onClick={!isLoading ? handleClick : null}
        >
          {isLoading ? 'Loading…' : 'View Visualization'}
        </button>
        <div className={styles.visualizer}>
          {/*Visualization*/}
        </div>
      </div>
      <div className={styles.textSection}>
        <p>
          {/*Explanation*/}
        </p>
      </div>
      <div className={styles.plotsSection}>
        <canvas id="energyPlot" className={styles.plot}></canvas>
      </div>
    </div>
  );
}

export default HomePage;

async function plotEnergyVsStep(logFile) {
  let steps = [];
  let totalEnergy = [];
  let insideData = false;

  // Simulating file read with mock data
  let data = `Step CPU TotEng KinEng PotEng E_vdwl E_coul E_long Temp Press Volume\n0 0 10\n1 1 9\n2 2 8\nLoop time of 0.1 on 1 procs for 100 steps with 100 atoms`; // Replace with actual file read

  data.split('\n').forEach(line => {
    if (line.includes("Step CPU TotEng KinEng PotEng E_vdwl E_coul E_long Temp Press Volume")) {
      insideData = true;
      return;
    }
    if (insideData && line.includes("Loop time of")) {
      insideData = false;
      return;
    }
    if (insideData) {
      let columns = line.trim().split(/\s+/);
      let step = parseInt(columns[0]);
      let energy = parseFloat(columns[2]); 
      steps.push(step);
      totalEnergy.push(energy);
    }
  });

  
  const canvas = document.getElementById('energyPlot');
  const ctx = canvas.getContext('2d');
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  ctx.beginPath();
  ctx.moveTo(0, canvas.height - totalEnergy[0] * 10); 
  for (let i = 1; i < steps.length; i++) {
    ctx.lineTo(steps[i] * (canvas.width / steps.length), canvas.height - totalEnergy[i] * 10); 
  }
  ctx.stroke();
}
