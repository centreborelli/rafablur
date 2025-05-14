BIN = image_gauss image_gaussian_sampler

LDLIBS = -lm

all: $(BIN)

jirafa_blur.o: CPPFLAGS += -DIMAGE_GAUSS_OMIT_MAIN
jirafa_blur.o: image_gauss.c ; $(COMPILE.c) $^ -o $@

clean: ; $(RM) $(BIN) *.o
