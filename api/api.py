from flask import Flask
from flask_restful import Resource, Api, reqparse

app = Flask(__name__)
api = Api(app)

class HelloWorld(Resource):
    def get(self):
        return {'hello': 'world'}

class Test1(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('name1', type=int, help='a test input', location='args')
        args = parser.parse_args()
        return {'the input for name1': args['name1']}

api.add_resource(HelloWorld, '/')
api.add_resource(Test1, '/test1')

if __name__ == '__main__':
    app.run(debug=True)
