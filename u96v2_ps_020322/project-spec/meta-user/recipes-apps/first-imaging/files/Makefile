APP = first-imaging

# Add any other object files to this list below
APP_OBJS = first-imaging.o
APP_OBJS += DataHandler.o 
APP_OBJS += ImageProcessor.o
APP_OBJS += NeonMath.o

#FFLAGS = -O3
#LDFLAGS = $(FFLAGS)

all: build

BLT:= 1.0.122122111111
build: $(APP)

$(APP): $(APP_OBJS)  #<Target> : <Dependencies> ,<Recipe>-Tab Indent $@:Automatic variables(Target name), CXX:g++(macro)
	$(CXX) -o $@ $(APP_OBJS) $(LDFLAGS) $(LDLIBS) -lm
clean:
	rm -f $(APP) *.o   