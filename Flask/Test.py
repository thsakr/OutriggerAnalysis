from flask import Flask, request, render_template

import math     as mt
import numpy	as np

import datetime as dt


class   Outrigger:
    
    def __init__(self, data):
        self.sc = 3


    def Out_Analysis   (self):

        return self.sc


app = Flask(__name__)


@app.route  ('/')

def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])


def analyze():
    My_Data = request.form.to_dict()
    #data = request.json
    outrigger = Outrigger(My_Data)
    result = outrigger.Out_Analysis()
    return  str( result)

if __name__ == '__main__':
     app.run(debug=True)

