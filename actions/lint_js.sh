# Copyright 2020 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Run the closure compiler, having it lint code and produce warnings during compilation.
# Takes in one argument: the entry point for the compiler.
# Assumes that all prerequisite Javascript has been generated, 
# and that the Google Closure compiler has been installed as a Node package,
# all in the directory that this script is run from.

# In order to compile properly, we must provide all the source files necessary, 
# and make sure the compiler skips over files that aren't needed as dependencies.
# Because we are only checking code, we skip generating output and running optimizations.
# We also produce warnings only for our written source code 
# (by suppressing warnings for third-party libraries or machine generated code).

set -ex

if [ "$#" -ne 1 ]; then
  echo "entry point required"
  exit 1
fi
ENTRY_POINT="$1"

java -jar node_modules/google-closure-compiler-java/compiler.jar \
    --entry_point "${ENTRY_POINT}" \
    --js 'closure-library-20200517/**.js' \
    --js 'protobuf-3.13.0/js/**.js' \
    --js 'gen/js/proto/ord/**.js'  \
    --js 'js/**.js' \
    --dependency_mode PRUNE \
    --checks_only \
    --jscomp_error lintChecks \
    --hide_warnings_for closure-library-20200517 \
    --hide_warnings_for protobuf-3.13.0/js \
    --hide_warnings_for gen/js/proto/ord
