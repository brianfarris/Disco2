APP      = disco

SRCEXT   = c
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
H55 = /home/install/app/hdf5-1.6_intel_mpi
#H55 = /share/apps/hdf5/1.8.2/openmpi/intel
endif
ifeq ($(UNAME),Darwin)
#H55 = /usr/local/hdf5
H55 = $(HOME)/Coding/Libraries/hdf5-1.8.9-par
endif

SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)')
SRCDIRS := $(shell find . -name '*.$(SRCEXT)' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))

DEBUG    = -g
INCLUDES = -I$(H55)/include
CFLAGS   = -O3 -c $(INCLUDES)
#CFLAGS   = -c $(DEBUG) $(INCLUDES)
LDFLAGS  = -lm -lz -L$(H55)/lib -lhdf5

CC       = mpicc

.PHONY: all clean distclean


all: $(BINDIR)/$(APP)

$(BINDIR)/$(APP): buildrepo $(OBJS)
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) $(OBJS) $(LDFLAGS) -o $@

$(OBJDIR)/%.o: %.$(SRCEXT)
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

clean:
	$(RM) -r $(OBJDIR)

distclean: clean
	$(RM) -r $(BINDIR)

buildrepo:
	@$(call make-repo)

define make-repo
   for dir in $(SRCDIRS); \
   do \
	mkdir -p $(OBJDIR)/$$dir; \
   done
endef
