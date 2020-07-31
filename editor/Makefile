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

all : prep py js ketcher-all

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

ketcher-all: ketcher ketcher-pull ketcher/node_modules ketcher/dist

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
	  --proto_path=.. \
	  $<

../build/lib/ord_schema/proto/%_pb2.py : ../ord_schema/proto/%.proto
	$(error $< not found--run setup?)

ketcher-update: ketcher-pull ketcher/node_modules ketcher/dist

ketcher-pull: ketcher
	cd ketcher && git pull

ketcher:
	git clone https://github.com/Open-Reaction-Database/ketcher.git

# Installs on changes to package.json or package-lock.json, not node_modules.
# (corruptions to node-modules, for example, will not be fixed; use make clean)
# TODO the existence of ketcher/node_modules does NOT imply a proper install,
# maybe use touch to create a file to signal and check for proper install?
ketcher/node_modules: ketcher/package.json ketcher/package-lock.json
	# Install Node.js modules.
	# Might have to sudo npm install, for phantom.

	# If make is run as root, then the ketcher folder has root perms.
	# Then, to install phantomjs properly, npm requires to be run as root.
	# If make is not run as root, then npm can be run normally.
	# Either way, serve.sh can be run as a regular user, since we only need to read protected files.
	@if [ $$USER = "root" ] ; then \
		cd ketcher && sudo npm install; \
	else \
		cd ketcher && npm install; \
	fi;

# Builds on changes to anything in ketcher/, that aren't in ketcher/node_modules or /ketcher_dist
# (so essentially, changes to any editable source code)
# TODO ketcher/dist existing doesnt imply finished build
ketcher/dist: @find ./ketcher/* ! -path "./ketcher/node_modules*" ! -path "./ketcher/dist*"
	# Build Ketcher from source.
	cd ketcher && npm run build

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

	mkdir package/ord/ketcher
	cp -r ketcher/dist package/ord/ketcher

	mkdir package/ord/db

	tar --create --directory=package --gzip --file=package/ord.tgz ord

test : all
	node js/test.js

clean :
	rm -rf gen
	rm -rf ketcher/node_modules
	rm -rf ketcher/dist
	rm -rf package
