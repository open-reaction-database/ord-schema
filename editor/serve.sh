#!/bin/bash

export PYTHONPATH=py:gen/py
export FLASK_APP=serve.py
export FLASK_ENV=development

python -m flask run
# Caution.
# python -m flask run --host=0.0.0.0 --port=80
