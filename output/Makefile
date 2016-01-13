dataFiles = $(wildcard *.dat)

define MAKEPLOT
$(1).png: $(1) draw-image.plt
	gnuplot \
		-e "filename='$(strip $1)'; filetitle='$(subst -, ,$(basename $(strip $1)))'" \
		./draw-image.plt
endef

all: $(addsuffix .png, $(dataFiles))

$(foreach d, $(dataFiles), $(eval $(call MAKEPLOT, $(d))))