#!/bin/bash

TEST_DIR=${PWD}/test/src

if [[ ! -d ${TEST_DIR} ]]; then
	echo "Error: script must be executed from the fast_indexer_dev root dir"
	exit 1
fi

if [[ ! -d ${BUILD_DIR} ]]; then
        mkdir test-build
fi

cd ${TEST_DIR}
make
