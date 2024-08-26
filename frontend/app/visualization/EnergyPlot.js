'use client'

import React from "react";
import { useEffect, useState } from 'react';
import Chart from 'chart.js/auto';
import { Line } from "react-chartjs-2";
import regression from 'regression';

const options = {
    plugins: {
        title: {
            text: "Total Energy",
            display: true,
            font: {
                size: 20,
            },
            padding: {
                top: 10,
                bottom: 30
            }
        },
        legend: {
            display: true,
            position: 'top',
        },
    },
    scales: {
        x: {
            title: {
                display: true,
                text: 'Step',
                font: {
                    size: 14,
                }
            }
        },
        y: {
            title: {
                display: true,
                text: 'Total Energy',
                font: {
                    size: 14,
                }
            }
        }
    }
};

// Chart.register(annotationPlugin);
// Chart.register('chartjs-plugin-regression');

export default function EnergyPlot({ visualId } ) {
    const [log, setLog] = useState('');

    let steps = [];
    let totalEnergy = [];
    let insideData = false;

    useEffect(() => {
        const fetchData = async () => {
          try {
            const response = await fetch('/api/getfile/logfile'+ visualId +'.log');
            const data = await response.json();
            const content = data.content;
            setLog(content);
          } catch (error) {
            console.error('Error fetching data:', error);
          }
        };
    
        fetchData();
      }, []);
    
    log.split('\n').forEach(line => {
        if (line.includes("Step          CPU")) {
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

    const coords = steps.map((el, index)=> [el, totalEnergy[index]]);

    const clean_data = coords
        .filter(({ x, y }) => {
        return (
            typeof x === typeof y &&  // filter out one string & one number
            !isNaN(x) &&              // filter out `NaN`
            !isNaN(y) &&
            Math.abs(x) !== Infinity && 
            Math.abs(y) !== Infinity
        );
        })
        .map(({ x, y }) => {
        return [x, y];             // we need a list of [[x1, y1], [x2, y2], ...]
        });

    const my_regression = regression.linear(
    clean_data
    );

    const useful_points = my_regression.points.map(([x, y]) => {
    return y;    
    // or {x, y}, depending on whether you just want y-coords for a 'linear' plot
    // or x&y for a 'scatterplot'
    })

    const data = {
        labels: steps,
        datasets: [
            {
                label: "Step vs Total Energy",
                backgroundColor: "rgba(54, 162, 235, 0.5)",
                borderColor: "rgba(54, 162, 235, 1)",
                data: totalEnergy,
                tension: 0.4,
            }
        ]
    };
    return (
        <>
            <Line data={data} options={options} />
        </>
    );
};
