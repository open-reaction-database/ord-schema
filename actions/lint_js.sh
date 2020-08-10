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

if [ "$#" -ne 1 ]; then
  echo "entry point required"
  exit 1
fi
ENTRY_POINT="$1"

java -jar node_modules/google-closure-compiler-java/compiler.jar \
    --entry_point $ENTRY_POINT \
    --js 'closure-library-20200517/**.js' \
    --js 'protobuf/js/**.js' \
    --js 'gen/js/proto/ord/**.js'  \
    --js 'js/**.js' \
    --dependency_mode PRUNE \
    --checks_only \
    --jscomp_warning lintChecks \
    --hide_warnings_for closure-library-20200517 \
    --hide_warnings_for protobuf/js \
    --hide_warnings_for gen/js/proto/ord
