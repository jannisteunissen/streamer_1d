SRC_DIRS	:= src
CREATE_DIRS	:= output

.PHONY:	all clean

all: | $(CREATE_DIRS)
	$(MAKE) -C src

clean:
	$(MAKE) -C src clean

$(CREATE_DIRS):
		mkdir -p $@
