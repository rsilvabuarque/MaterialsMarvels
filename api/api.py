from flask import Flask
from flask_restful import Resource, Api, reqparse
import os

app = Flask(__name__)
api = Api(app)

# For CORS
@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE')
    return response

class HelloWorld(Resource):
    def get(self):
        return {'hello': 'world'}

class Visualize(Resource):
    def post(self):
        # Parse the input data
        parser = reqparse.RequestParser()
        parser.add_argument('molfile', type=str, help='The molfile input (v3000)')
        args = parser.parse_args()

        # Find an unused filename from molfile1.mol, molfile2.mol, ...
        if not os.path.exists('temp'):
            os.makedirs('temp')
        
        i = 1
        while os.path.exists(f'temp/molfile{i}.mol'):
            i += 1

        # Save the molfile to a new file in the temp directory
        with open(f'temp/molfile{i}.mol', 'w') as f:
            f.write(args['molfile'])

        # Run the scripts to generate output files

        # Read output files and return the data
        with open(f'temp/molfile{i}.mol', 'r') as f:
            output = f.read()

        # Delete all the temp files created in the process
        os.remove(f'temp/molfile{i}.mol')

        return {'output': output}

api.add_resource(HelloWorld, '/api/')
api.add_resource(Visualize, '/api/visualize')

if __name__ == '__main__':
    app.run(port=7000, debug=True)
