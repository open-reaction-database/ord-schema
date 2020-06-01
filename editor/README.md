# A Reaction Editor

A web app for browsing, modifying, and creating submissions to the
Open Reactions Databse.

## Getting Started

```
$ make
$ ./serve.sh
```

This starts service at [http://localhost:5000/](http://localhost:5000/).

## Dependencies

The build needs:
* the `protoc` protobuf compiler;
* the protobuf runtime libraries for python and Javascript; and
* the Closure Library for Javascript.

Serving depends on:
* the protobuf python package; and
* the Flask python web framework.

And everything requires Python 3.

The editor has been tested with [protobuf
3.12.2](https://github.com/protocolbuffers/protobuf/releases) and [Closure
v20200517](https://github.com/google/closure-library/releases/). Unpack these
in this directory so that make can find them, and compile protobuf to get
the protoc compiler.

```
cd protobuf-3.12.2
./configure && make
```

(As an alternative, you can edit the Makefile to use a protoc binary you
already have. But make still needs the protobuf distribution for its JS common
files.)

To install the python packages,

```
$ pip install flask
$ pip install protobuf
```

The Closure compiler fails to resolve dependencies for test files in the
protobuf distribution. Feel free to delete any protobuf files Closure doesn't
like, e.g.

```
$ rm -rf protobuf-3.12.2/js/*test*
$ rm -rf protobuf-3.12.2/js/binary/*test*
$ rm -rf protobuf-3.12.2/js/compatibility_tests
$ rm -rf protobuf-3.12.2/js/experimental
```

## Testing and Validation

There is an end-to-end test that reads a densely populated Dataset into the
editor, writes it back out again, and asserts that the proto is unchanged.

```
$ ./serve.sh
$ make test
```

The test depends on [Node.js](https://nodejs.org/en/download/) and the
Puppeteer headless Chrome package.

```
$ npm i puppeteer
```

## Packaging and Deployment

There is a make target that bundles all required runtime resources for
deployment to a real server. The output appears at package/ord.tgz.

```
$ make package
```

The output includes the serve.sh script, but you probably want to use a better
web server instead.

## How it Works

Datasets to be edited are stored in the `db/` directory as `.pbtxt` files (protobuf
text format).

When you load a Reaction in the editor, the editor reads the entire Dataset
from the server and maps the content of the chosen Reaction onto the DOM. When
you save, it reverses the process. Javascript functions called "load..." go
from Reaction to DOM, and functions called "unload..." do the reverse.

Some fields like images are too big to go in the DOM. These fields have type
"bytes". When a bytes field is loaded, the load function creates a random token
for it, puts the token in the DOM element, and saves the bytes data in a table.
The unload function reads back the token and uses it to restore the bytes data.

You modify a bytes field by uploading a file. Since Javascript does not have
access to files, the model in this case is different. When you save, the bytes
field in the Reaction is assigned its token and the file is uploaded in a
separate connection with the same token. The server merges the Reaction and its
upload after both are available.

Since the backing store is a file system, writes happen in multiple steps, and
there are no incremental updates, users must take care.

* There is no database. If you kill the server while it's writing you will corrupt the entire Dataset. Forget about concurrent access.

* Really large Datasets need to pass both ways on the network every time you edit them, including all their images.

* Binary data are added by uploading files, and Javascript is not allowed to access files. This means that if you attach a new image to a Reaction, the two must be uploaded separately and merged on the server. This is another opportunity for corruption if the process is interrupted.
