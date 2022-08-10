#!/bin/bash

ROOT=${ROOT:-.}
EXE=(${EXE:-${ROOT}/test-build/xgandalf})

if [[ ! -x ${EXE[0]} ]]; then
	echo "$EXE not found, execute from root dir or set ROOT variable"
	exit 1
fi

if [[ -z "$LD_LIBRARY_PATH" ]]; then
	LD_LIBRARY_PATH=${ROOT}/xgandalf-build
else
	LD_LIBRARY_PATH+=:${ROOT}/xgandalf-build
fi

eval ${EXE[@]} "$@"
