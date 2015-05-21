DX = doxygen
DOC = docs/Doxyfile

default: release

all: debug release

debug:
	@echo [MAKE] $@
	@$(MAKE) -C debug

release:
	@echo [MAKE] $@
	@$(MAKE) -C release

docs:
	@echo [DX] generating documentation
	@$(DX) $(DOC) > /dev/null

clean:
	@echo [MAKE] clean
	@$(MAKE) -C debug clean
	@$(MAKE) -C release clean


.PHONY: default all install debug release docs clean uninstall
