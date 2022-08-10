#!/bin/bash

XGANDALF_DIR=${PWD}/xgandalf

if [[ ! -d ${XGANDALF_DIR} ]]; then
	echo "Error: script must be executed from the fast_indexer_dev root dir"
	exit 1
fi

usage () {
	echo "Use: $0 (con|com|cc|i|ci|cci) [meson args...]"
}

error () {
	usage
	echo "Error: $1"
	exit 1
}

if (($# < 1)); then
	error "no args"
fi

CONF=0
COMP=0
INST=0

case "$1" in
	cci) CONF=1; COMP=1; INST=1;;
	ci) COMP=1; INST=1;;
	i) INST=1;;
	cc) CONF=1; COMP=1;;
	com) COMP=1;;
	con) CONF=1;;
	*) error "unknown command <$1>"
esac

shift

BUILD_DIR=${PWD}/xgandalf-build
INSTALL_DIR=${PWD}/xgandalf-install

if ((CONF == 1)); then
	if [[ -d ${BUILD_DIR} ]]; then
		echo "cleaning ${BUILD_DIR} ..."
		pushd ${BUILD_DIR}
		rm -rf *
		rm .ninja*
		popd
	fi	
	pushd ${XGANDALF_DIR}
	echo "configure args=<$@> ..."
	if [[ -d ${BUILD_DIR}/meson-private ]]; then
		RECONFIGURE_OPTION=--reconfigure
	fi
	meson setup ${BUILD_DIR} --prefix=${INSTALL_DIR} ${RECONFIGURE_OPTION} "$@"
	popd
fi

pushd ${BUILD_DIR}
	if ((COMP == 1)); then
		echo "compile..."
		ninja -v
	fi

	if ((INST == 1)); then
		echo "install..."
		ninja install
		pushd ${INSTALL_DIR}
		if [[ ! -d lib/x86_64-linux-gnu && -d lib64 ]]; then
			mkdir -p lib
			pushd lib
			ln -sf ../lib64 x86_64-linux-gnu
			popd
		fi
		if [[ ! -d include/Eigen && -d include/eigen3 ]]; then
			ln -sf eigen3/Eigen include/Eigen
		fi
		popd
	fi
popd
