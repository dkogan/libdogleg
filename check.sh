#!/bin/zsh

# This is a rudimentary test. Not exhaustive. Makes sure the simple problem
# defined in sample.c can be solved. Both with sparse and dense callbacks


GREEN="\x1b[32m"
RED="\x1b[31m"
COLOR_RESET="\x1b[0m"

CMDS=('./sample --check sparse'
      './sample --check dense'
      './sample --check dense-products-packed-upper'
      './sample --check dense-products-unpacked'
      './test-misc')

for cmd ($CMDS) {
    echo "Running:   $cmd";
    ${=cmd} || anyfailed=1
}

if [[ -n "$anyfailed" ]]; then
    echo -e $RED"Some tests failed"$COLOR_RESET
    exit 1;
fi

exit 0
