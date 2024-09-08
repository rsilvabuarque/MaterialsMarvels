'use client'

import React, { useEffect, useState } from 'react';
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

export default function EnergyPlot({ visualId }) {
    const [log, setLog] = useState('');
    const [sliderValue, setSliderValue] = useState(0);  
    const [maxSteps, setMaxSteps] = useState(100);      

    let steps = [];
    let totalEnergy = [];
    let insideData = false;

    useEffect(() => {
        const fetchData = async () => {
            try {
                const response = await fetch('/api/getfile/logfile' + visualId + '.log');
                const data = await response.json();
                const content = data.content;
                setLog(content);
            } catch (error) {
                console.error('Error fetching data:', error);
            }
        };

        fetchData();
    }, [visualId]);

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

    useEffect(() => {
        setMaxSteps(steps.length);
    }, [steps]);

    const visibleEnergy = totalEnergy.slice(0, sliderValue);
    const visibleSteps = steps.slice(0, sliderValue);

    const coords = visibleSteps.map((el, index) => [el, visibleEnergy[index]]);
    const polynomialRegression = regression.polynomial(coords, { order: 4, precision: 6 });
    const polynomialFitData = polynomialRegression.points.map(([x, y]) => ({ x, y }));

    const data = {
        labels: visibleSteps,
        datasets: [
            {
                label: "Step vs Total Energy",
                backgroundColor: "rgba(54, 162, 235, 0.5)",
                borderColor: "rgba(54, 162, 235, 1)",
                data: visibleEnergy,
                tension: 0.4,
                order: 1
            },
            {
                label: "Polynomial Fit",
                backgroundColor: "rgba(255, 99, 132, 0.5)",
                borderColor: "rgba(255, 99, 132, 1)",
                data: polynomialFitData.map(point => point.y), // Extract y values for the fit line
                borderDash: [5, 5], // Dashed line for distinction
                fill: false,
                pointRadius: 0, // No points, just the line
                tension: 0.4,
                order: 0
            }
        ]
    };

    return (
        <div>
            <input
                type="range"
                min="0"
                max={maxSteps}
                value={sliderValue}
                onChange={(e) => setSliderValue(parseInt(e.target.value))}
                style={{ width: '100%' }}
            />
            <p>Current Step: {sliderValue}</p>
            <Line data={data} options={options} />
        </div>
    );
};
