FROM python:3
WORKDIR /usr/src/app

RUN pip install Flask protobuf

# First run "make package".
COPY package/ord.tgz .
RUN tar zxf ord.tgz

WORKDIR /usr/src/app/ord
EXPOSE 5000
CMD ./serve.sh
