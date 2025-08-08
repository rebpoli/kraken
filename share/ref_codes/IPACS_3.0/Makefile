

build: .force cmake
	$(MAKE) -C $@

cmake:
	mkdir -p build
	cd build/$* && cmake ..

clean:
	rm -rf build

run:
	build/IPACS-V3.0

.force:
