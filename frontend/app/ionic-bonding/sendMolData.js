function sendMolData() {
  window.ketcher.getMolfile('v2000').then((molfile) => {
    // Send the molfile to the Flask API and retrieve the output after computation
    fetch('/api/visualize', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ 'molfile': molfile }),
    }).then((response) => {
      return response.json();
    }).then((data) => {
      // Simulation has been completed - redirect to the visualization page
      window.location.href = '/visualization?visualId=' + data.visualId;
    })
  });
  return new Promise((resolve) => setTimeout(resolve, 1000));
}

export default sendMolData;