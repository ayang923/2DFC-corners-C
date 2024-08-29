ACCELERATE_LIB = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers

CC = clang
CFLAGS = -Wall -Wextra -Iinclude/ -I$(ACCELERATE_LIB)

TARGETS = main
LIBS = q_patch_lib num_linalg_lib

OBJS = $(addprefix out/,$(LIBS:=.o)) $(addprefix out/,$(TARGETS:=.o)) 
BINS = $(addprefix bin/,$(TARGETS))

# find <dir> is the command to find files in a directory
# ! -name .gitignore tells find to ignore the .gitignore
# -type f only finds files
# -delete deletes all the files found
CLEAN_COMMAND = find out/ ! -name .gitignore -type f -delete && \
find bin/ ! -name .gitignore -type f -delete

all: $(BINS)

$(BINS): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) -framework Accelerate

out/%.o: src/%.c
	$(CC) $(CFLAGS) -c $^ -o $@ -framework Accelerate

clean:
	$(CLEAN_COMMAND)

# This special rule tells Make that "all", "clean", and "test" are rules
# that don't build a file.
.PHONY: all clean
# Tells Make not to delete the .o files after the executable is built
.PRECIOUS: out/%.o