#!/bin/bash

ROOT=${ROOT:-.}
EXE=${EXE:-${ROOT}/test-build/xgandalf}
LIBDIR=$(readlink -f $(dirname $(find ${ROOT}/xgandalf-install/lib -name libxgandalf.so)))

if [[ ! -x ${EXE[0]} ]]; then
	echo "$EXE not found, execute from root dir or set ROOT variable"
	exit 1
fi

if ! printenv LD_LIBRARY_PATH | grep -s xgandalf &>/dev/null; then
        if [[ -d ${LIBDIR} ]]; then
                if [[ -z "$LD_LIBRARY_PATH" ]]; then
                        LD_LIBRARY_PATH=${LIBDIR}
                else
                        LD_LIBRARY_PATH+=:${LIBDIR}
                fi
        fi
fi

export LD_LIBRARY_PATH

eval ${EXE[@]} "$@"
