BASE_DIR = ../..
INSTALL_DIR ?= $(BASE_DIR)/test-build
XGANDALF_DIR = $(BASE_DIR)/xgandalf-install
PKGCONF_FILE = $(XGANDALF_DIR)/lib/x86_64-linux-gnu/pkgconfig/xgandalf.pc

xgandalf:
	$(CXX) $(CXXFLAGS) $(shell pkg-config --cflags $(PKGCONF_FILE)) -o $(INSTALL_DIR)/$@ XGandalfPerfTest.cpp $(shell pkg-config --libs $(PKGCONF_FILE))

.phony: xgandalf
