DX = doxygen
DOC = docs/Doxyfile

NAME = cpppoa

default: release

all: debug release

debug:
	@echo [MAKE] $(NAME) $@
	@$(MAKE) -C debug

release:
	@echo [MAKE] $(NAME) $@
	@$(MAKE) -C release

docs:
	@echo [DX] generating documentation
	@$(DX) $(DOC) > /dev/null

clean:
	@echo [MAKE] clean
	@$(MAKE) -C debug clean
	@$(MAKE) -C release clean


.PHONY: default all install debug release docs clean uninstall
