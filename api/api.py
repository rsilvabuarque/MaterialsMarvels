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

class TempFileHandler(Resource):
    def get(self, filename):
        # parser = reqparse.RequestParser()
        # parser.add_argument('filename', type=str, help='The file name that you wish to retrieve from the temp directory')
        # args = parser.parse_args()

        # TODO - Add error handling for file not found
        # TODO - Prevent directory traversal attacks by changing how this is done
        with open(f'temp/{filename}', 'r') as f:
            output = f.read()
        
        return {"content": output}

class Visualize(Resource):
    def post(self):
        # Parse the input data
        parser = reqparse.RequestParser()
        parser.add_argument('molfile', type=str, help='The molfile input (v2000)')
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

class Calculate(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('num1', type=int, help='Number 1')
        parser.add_argument('num2', type=int, help='Number 2')
        args = parser.parse_args()
        result = args['num1'] + args['num2']
        return "The result is " + result


# test
api.add_resource(HelloWorld, '/api/')
api.add_resource(TempFileHandler, '/api/getfile/<string:filename>')
api.add_resource(Visualize, '/api/visualize')
api.add_resource(Calculate, '/api/calculate/')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000)
