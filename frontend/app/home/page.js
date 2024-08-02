'use client'

import styles from './page.module.css';
import EnergyPlot from './EnergyPlot';
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
                <Col style={{"borderColor": "navy", "border-style": "solid"}}><div dangerouslySetInnerHTML={{ __html: htmlContent }} /></Col>
                {/* <Col><iframe src="/api/getfile/visualization.html" frameborder="0"></iframe></Col> */}
                <Col style={{"borderColor": "navy", "border-style": "solid"}}><EnergyPlot /></Col>
            </Row>
            <Row>
                <Col style={{"borderColor": "navy", "border-style": "solid"}}><div>Explanation text <br /> Explanation text</div></Col>
                <Col style={{"borderColor": "navy", "border-style": "solid"}}><EnergyPlot /></Col>
            </Row>
            {/* <div className={styles.topPane}>test</div>
            <div className={styles.bottomPane}><EnergyPlot /></div> */}
        </Container>
        {/* <div>
        test 2
        </div> */}
        </>
    )
}

export default HomePage;