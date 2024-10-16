'use client'

import styles from './page.module.css';
import VariablePlot from './VariablePlot';
import VideoVisual from '../components/VideoVisual';
import Container from 'react-bootstrap/Container';
import Navigation from '../components/Navigation';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';
import { useEffect, useState } from 'react';
import { useSearchParams } from 'next/navigation'

export default function Page() {
    const [htmlContent, setHtmlContent] = useState('Loading...');
    const [sliderValue, setSliderValue] = useState(0);

    const searchParams = useSearchParams();
    const visualId = searchParams.get('visualId');

    const handleSliderChange = (value) => {
        setSliderValue(value); // Update sliderValue when it changes in VideoVisual
    };

    return (
        <>
            <Navigation />
            <Container className={styles.pageContainer}>
                <Row>
                    <Col className={styles.visualizationCol}>
                        <VideoVisual visualId={visualId} onProgressChange={handleSliderChange} />
                    </Col>
                    <Col className={styles.textCol}>
                        <div className={styles.explanationText}>Explanation text <br /> Explanation text</div>
                    </Col>
                </Row>
                <Row>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={2} variableName="Total Energy" variableUnit="eV" />
                    </Col>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={18} variableName="Temperature" variableUnit="K" />
                    </Col>
                </Row>
                <Row>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={4} variableName="Coulomb Energy" variableUnit="eV" />
                    </Col>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={5} variableName="Pair Energy" variableUnit="eV" />
                    </Col>
                </Row>
                <Row>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={6} variableName="Bond Energy" variableUnit="eV" />
                    </Col>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={7} variableName="Angle Energy" variableUnit="eV" />
                    </Col>
                </Row>
                <Row>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={8} variableName="Dihedral Energy" variableUnit="eV" />
                    </Col>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={9} variableName="Improper Energy" variableUnit="eV" />
                    </Col>
                </Row>
                <Row>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={10} variableName="Molecular Energy" variableUnit="eV" />
                    </Col>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={11} variableName="Long-Range Energy" variableUnit="eV" />
                    </Col>
                </Row>
                <Row>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={12} variableName="Tail Correction Energy" variableUnit="eV" />
                    </Col>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={3} variableName="Van der Waals Energy" variableUnit="eV" />
                    </Col>
                </Row>
                <Row>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={14} variableName="Energy Coupling" variableUnit="eV" />
                    </Col>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={15} variableName="Energy Conservation" variableUnit="eV" />
                    </Col>
                </Row>
                <Row>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={16} variableName="Kinetic Energy" variableUnit="eV" />
                    </Col>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={17} variableName="Potential Energy" variableUnit="eV" />
                    </Col>
                </Row>
                <Row>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={19} variableName="Pressure" variableUnit="atm" />
                    </Col>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={20} variableName="Volume" variableUnit="A^3" />
                    </Col>
                </Row>
                <Row>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={21} variableName="Density" variableUnit="g/cm^3" />
                    </Col>
                    <Col className={styles.plotCol}>
                        <VariablePlot visualId={visualId} sliderValue={sliderValue} variableIndex={13} variableName="Enthalpy" variableUnit="eV" />
                    </Col>
                </Row>

            </Container>
        </>
    );
}
