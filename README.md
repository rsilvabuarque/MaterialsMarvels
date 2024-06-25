# Self-Assembly Simulation Project

## Installation

- [Install](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm) `Node.js` and `npm`
- Install Python 3.10
- Install LAMMPS on the machine

### Install and Activate Venv

#### Windows
```bash
cd api
py -3 -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

#### MacOS/Linux
```bash
cd api
python3 -m venv .venv
. .venv/bin/activate
pip3 install -r requirements.txt
```

### Install Node.js packages
```bash
cd frontend
npm install
```

## Running the App
- Run the `api.py` file when in the `api` directory: `python3 api.py`
- Run `npm run dev` when in the `frontend` directory