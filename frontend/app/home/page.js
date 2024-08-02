'use client'

import styles from './page.module.css';
import EnergyPlot from './EnergyPlot';
import Visualization from './Visualization';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';
import { useEffect, useState } from 'react';

const HomePage = () => {
  const [htmlContent, setHtmlContent] = useState('');

  useEffect(() => {
    const fetchHtml = async () => {
      try {
        const response = await fetch('/api/getfile/visualization.html');
        if (response.ok) {
          const resp = await response.json();
          setHtmlContent(resp.content);
        } else {
          console.error('Failed to fetch HTML content');
        }
      } catch (error) {
        console.error('Error fetching HTML content:', error);
      }
    };

    fetchHtml();
  }, []);

  return (
    <>
      <Container>
        <Row>
          <Col style={{ "borderColor": "navy", "borderStyle": "solid" }}><Visualization markup={htmlContent}></Visualization></Col>
          <Col style={{ "borderColor": "navy", "borderStyle": "solid" }}><EnergyPlot /></Col>
        </Row>
        <Row>
          <Col style={{ "borderColor": "navy", "borderStyle": "solid" }}><div>Explanation text <br /> Explanation text</div></Col>
          <Col style={{ "borderColor": "navy", "borderStyle": "solid" }}><EnergyPlot /></Col>
        </Row>
      </Container>
    </>
  )
}

export default HomePage;