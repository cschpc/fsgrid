MAKEFLAGS += --silent

project_name:=fsgrid
root_dir:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
source_dir:=$(root_dir)
build_dir_prefix:=$(root_dir)/build

.PHONY: all test clean

.ONESHELL:
all:
	build_type=Release
	build_dir=$(build_dir_prefix)/$${build_type}
	cmake \
		-B $${build_dir} \
		-S $(source_dir) \
		-DCMAKE_BUILD_TYPE:STRING=$${build_type} \
		-Dproject_name=$(project_name) \
		-DMPI_TEST_NUM_PROCS=16 \
		-DCMAKE_EXPORT_COMPILE_COMMANDS=ON
	cmake \
		--build $${build_dir} \
		-j8
	cp $${build_dir}/compile_commands.json $(source_dir)

clean:
	rm -rf $(build_dir_prefix) $(project_name)

test: all
	ctest --test-dir $(build_dir_prefix)/Release --output-on-failure
