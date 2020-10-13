# Reaction Editor

A web app for browsing, modifying, and creating submissions to the Open
Reaction Database.

## Getting Started

From the ord-schema project root, 
```
$ export ORD_EDITOR_POSTGRES_PASSWORD=########
$ export ORD_EDITOR_MOUNT=$HOME/ord-editor-postgres
$ docker build -f editor/Dockerfile -t openreactiondatabase/ord-editor .
$ docker-compose -f editor/docker-compose.yml up
```
That builds the editor, launches it in a container, and points it at a second
container running Postgres where the editor keeps its state. You can see the
editor at [http://localhost:5000/](http://localhost:5000).

The first time you do this runs, you must initialize the Postgres schema.
```
psql -p 5432 -h localhost -U postgres -f schema.sql
```
You may also want to import some data to get started. The migration script
slurps the contents of the db/ directory.
```
export PGPASSWORD=########
./py/migrate.py
```

## How it Works

When you load a reaction in the editor, the editor reads the entire dataset
from the server and maps the content of the chosen reaction onto the DOM. When
you save, it reverses the process. Javascript functions called "load..." go
from reaction to DOM, and functions called "unload..." do the reverse.

Some fields like images are too big to go in the DOM. These fields have type
"bytes". When a bytes field is loaded, the load function creates a random token
for it, puts the token in the DOM element, and saves the bytes data in a table.
The unload function reads back the token and uses it to restore the bytes data.

You modify a bytes field by uploading a file. Since Javascript does not have
access to files, the model in this case is different. When you save, the bytes
field in the reaction is assigned its token and the file is uploaded in a
separate connection with the same token. The server merges the reaction and its
upload after both are available.

Really large datasets need to pass both ways on the network every time you save
them, including all their images.

## Development

You can get lightweight iterations by compiling and running the editor outside
Docker.

First install and run Postgres, or start its Docker image with the
"docker-compose" command above. Then set up the editor dependencies below. When
that's done, you can start the editor in a debug mode where modifications to
python or JS can be deployed quickly.
```
$ make
$ export POSTGRES_PASSWORD=########
$ export FLASK_ENV=development
$ ./serve.sh --port=5001
```
This starts service at [http://localhost:5001/](http://localhost:5001/).

### Dependencies

There are several.

#### ORD-Schema Environment

Activate Conda and install the setup.py modules by following the top-level
[instructions](https://github.com/Open-Reaction-Database/ord-schema/blob/main/README.md).

#### Closure JS Library

Unpack closure-library in this directory so that Closure can find it too.
```
$ url=https://github.com/google/closure-library/archive/v20200517.tar.gz
$ wget $url -O - | tar zxf -
```
The editor has been tested with [Closure
v20200517](https://github.com/google/closure-library/releases/).

#### Protobuf JS Runtime

Unpack the archive in the editor directory so the Closure compiler can find it.
```
$ url=https://github.com/protocolbuffers/protobuf/releases/download/v3.13.0/protobuf-js-3.13.0.tar.gz
$ wget $url -O - | tar zxf -
```
At the time of writing, [version
3.13.0](https://github.com/protocolbuffers/protobuf/releases/tag/v3.13.0)
matches the protoc compiler in pip's protoc-wheel-0 distribution referenced
from ../requirements.txt.

#### Ketcher

Ketcher draws molecule diagrams. First install Node.js and npm
([instructions](https://nodejs.org/en/download/)). Then, in this directory,
```
$ git clone git@github.com:Open-Reaction-Database/ketcher.git
$ cd ketcher
$ npm install && npm run build
```

## Testing and Validation

There is an end-to-end test that reads a densely populated dataset into the
editor, writes it back out again, and asserts that the proto is unchanged.

Start the editor as described above and then
```
$ make test
```

The test depends on [Node.js](https://nodejs.org/en/download/) and the
Puppeteer headless Chrome package.
```
$ npm i puppeteer
```

There is also an Abseil test. Start the editor as described above and then
```
$ export POSTGRES_PASSWORD=########
python ./py/serve_test.py
```
Note that this test exercises python code in the test process, not the editor's
Docker container. This test may behave differently from the editor as it
appears in a browser if the two become out of sync.
