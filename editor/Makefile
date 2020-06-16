# Use HEAD: git clone git@github.com:protocolbuffers/protobuf.git
protobuf = protobuf

# cd protobuf && ./autogen.sh && ./configure && make
# protoc = protobuf/src/protoc

# Instead of compiling, use these prebuilt protoc binaries available at
# https://storage.googleapis.com/ord-editor-test/editor_test_protobuf_1dae8fdd.tar
# The hex identifies the GitHub commit of the build.
UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
  protoc = protobuf/bin/protoc.elf
else ifeq ($(UNAME), Darwin)
  protoc = protobuf/bin/protoc.mach-o
else
  $(error protoc not found)
endif

# https://github.com/google/closure-library/releases/
closure = closure-library-20200517

all : prep py js

prep :
	# Ensure output directories exist.
	mkdir -p gen/js/ord gen/js/proto/ord
	# Prevent Closure from compiling test code.
	rm -rf $(protobuf)/js/binary/*test* $(protobuf)/js/*test* $(protobuf)/js/experimental

js : gen/js/proto/ord/dataset.js \
	   gen/js/proto/ord/reaction.js \
	   gen/js/ord/dataset.js \
	   gen/js/ord/reaction.js \

py : ../build/lib/ord_schema/proto/dataset_pb2.py \
	   ../build/lib/ord_schema/proto/reaction_pb2.py

gen/js/ord/dataset.js : js/dataset.js gen/js/proto/ord/dataset.js
	$(closure)/closure/bin/build/closurebuilder.py \
	  --root=$(closure) \
	  --root=$(protobuf)/js \
	  --root=gen/js/proto/ord \
	  --root=js \
	  --namespace=ord.dataset \
	  --output_mode=script \
	  --output_file=gen/js/ord/dataset.js \
	  || (rm -f gen/js/ord/dataset.js && false)

gen/js/ord/reaction.js : js/*.js gen/js/proto/ord/reaction.js
	$(closure)/closure/bin/build/closurebuilder.py \
	  --root=$(closure) \
	  --root=$(protobuf)/js \
	  --root=gen/js/proto/ord \
	  --root=js \
	  --namespace=ord.reaction \
	  --output_mode=script \
	  --output_file=gen/js/ord/reaction.js \
	  || (rm -f gen/js/ord/reaction.js && false)

gen/js/proto/ord/%.js : ../ord_schema/proto/%.proto
	$(protoc) \
	  --experimental_allow_proto3_optional \
	  --js_out=binary:gen/js/proto/ord \
	  --proto_path=../.. \
	  $<

../build/lib/ord_schema/proto/%_pb2.py : ../ord_schema/proto/%.proto
	$(error $< not found--run setup?)

package : all
	rm -rf package
	mkdir -p package/ord

	cp -r \
	  css \
	  gen \
	  html \
	  py \
	  serve.sh \
	  package/ord

	cp -r ../build/lib/ord_schema package/ord/py

	tar --create --directory=package --gzip --file=package/ord.tgz ord

test : all
	node js/test.js

clean :
	rm -rf gen
	rm -rf package
