APP = first-imaging

# Add any other object files to this list below
APP_OBJS = first-imaging.o
APP_OBJS += DataHandler.o 
APP_OBJS += ImageProcessor.o
APP_OBJS += NeonMath.o

all: build

build: $(APP)

$(APP): $(APP_OBJS)  #<Target> : <Dependencies> ,<Recipe>-Tab Indent $@:Automatic variables(Target name), CXX:g++(macro)
	$(CXX) -o $@ $(APP_OBJS) $(LDFLAGS) $(LDLIBS) -lm
clean:
	rm -f $(APP) *.o   