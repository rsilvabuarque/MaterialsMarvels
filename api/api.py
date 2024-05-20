from flask import Flask
from flask_restful import Resource, Api, reqparse

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

class Test1(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('name1', type=int, help='a test input', location='args')
        args = parser.parse_args()
        return {'the input for name1': args['name1']}

api.add_resource(HelloWorld, '/api/')
api.add_resource(Test1, '/api/test1')

if __name__ == '__main__':
    app.run(port=7000, debug=True)
