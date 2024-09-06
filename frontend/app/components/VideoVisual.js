import React, { useRef, useState } from 'react';
import FormRange from 'react-bootstrap/FormRange';
import Button from 'react-bootstrap/Button';

export default function VideoVisual({ filename }) {
    const vidRef = useRef(null);
    const [progress, setProgress] = useState(0);
    const [isPlaying, setIsPlaying] = useState(false);

    const handlePlayPause = () => {
        if (isPlaying) {
            vidRef.current.pause();
        } else {
            vidRef.current.play();
        }
        setIsPlaying(!isPlaying);
    };

    const handlePlayVideo = (event) => {
        const value = event.target.value;
        const duration = vidRef.current.duration;
        
        // Set the video's current time based on the slider value
        vidRef.current.currentTime = (value / 100) * duration;
        setProgress(value);
    };

    const handleTimeUpdate = () => {
        const currentTime = vidRef.current.currentTime;
        const duration = vidRef.current.duration;
        const currentProgress = (currentTime / duration) * 100;
        setProgress(currentProgress);
    };

    return (
        <>
            <video 
                ref={vidRef} 
                height="80%" 
                width="100%" 
                controls
                onTimeUpdate={handleTimeUpdate}
            >
                <source src={filename} type="video/mp4" />
            </video>

            <div>
                <span style={{ textAlign: 'center' }}>Current Step: {Math.round(progress)}</span>
            </div>
            
            {/* Flexbox container to align play button and slider */}
            <div style={{ display: 'flex', alignItems: 'center', marginTop: '10px' }}>
                <Button onClick={handlePlayPause} style={{ marginRight: '10px' }}>
                    {isPlaying ? 'Pause' : 'Play'}
                </Button>
                <FormRange 
                    onChange={handlePlayVideo} 
                    value={progress} 
                    min="0" 
                    max="100" 
                    step="1"
                    style={{ flexGrow: 1 }} // Make the slider grow to fill available space
                />
            </div>
        </>
    );
}
