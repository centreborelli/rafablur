BIN = image_gauss image_gaussian_sampler

LDLIBS = -lm

all: $(BIN)

clean: ; $(RM) $(BIN)
