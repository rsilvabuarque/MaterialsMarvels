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
      // Received the data from the Flask API
      window.alert("Received the data from the Flask API");
      console.log(data.output);
    })
  });
  return new Promise((resolve) => setTimeout(resolve, 1000));
}

export default sendMolData;