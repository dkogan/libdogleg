#!/bin/bash

# This is a rudimentary test. Not exhaustive. Makes sure the simple problem
# defined in sample.c can be solved. Both with sparse and dense callbacks


GREEN="\x1b[32m"
RED="\x1b[31m"
COLOR_RESET="\x1b[0m"


./sample --check sparse                       || anyfailed=1
./sample --check dense                        || anyfailed=1
./sample --check dense-products-packed-upper  || anyfailed=1
./sample --check dense-products-unpacked      || anyfailed=1

./test-misc || anyfailed=1

if [[ -n "$anyfailed" ]]; then
    echo -e $RED"Some tests failed"$COLOR_RESET
    exit 1;
fi

exit 0
