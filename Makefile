CXXFLAGS+=-O3 -ffast-math -fPIC -DPIC `pkg-config lv2core lv2-plugin fftw3f --cflags --libs`
LDFLAGS +=-shared `pkg-config lv2core lv2-plugin fftw3f --libs`

BUNDLE = p0ss3ss0r.lv2
INSTALL_DIR = $(DESTDIR)/usr/lib/lv2


all: $(BUNDLE)

$(BUNDLE): manifest.ttl p0ss3ss0r.ttl libp0ss3ss0r.so
	rm -rf $(BUNDLE)
	mkdir $(BUNDLE)
	cp manifest.ttl p0ss3ss0r.ttl libp0ss3ss0r.so $(BUNDLE)
	
p0ss3ss0r6.o: p0ss3ss0r.cpp p0ss3ss0r.hpp p0ss3ss0r.peg
	$(CXX) -c $(CXXFLAGS) p0ss3ss0r.cpp -o p0ss3ss0r.o
	
p0ss3ss0r.peg: p0ss3ss0r.ttl
	lv2peg p0ss3ss0r.ttl p0ss3ss0r.peg

install: $(BUNDLE)
	mkdir -p $(INSTALL_DIR)
	rm -rf $(INSTALL_DIR)/$(BUNDLE)
	cp -R $(BUNDLE) $(INSTALL_DIR)

uninstall:
	rm -rf $(INSTALL_DIR)/$(BUNDLE)
	
libp0ss3ss0r.so: p0ss3ss0r.o
	$(CXX) p0ss3ss0r.o $(LDFLAGS) -o libp0ss3ss0r.so
	
clean:
	rm rm -rf $(BUNDLE) *.so *.o *.peg
	
.phony: clean all install
