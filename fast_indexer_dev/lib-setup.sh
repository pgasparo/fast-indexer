if [[ -d ${PWD}/xgandalf-build ]]; then
	export LD_LIBRARY_PATH=${PWD}/xgandalf-build${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
else
	echo "Error: run this from the fast_indexer_dev root directory"
fi
