MAKEFILE_IN = $(PWD)/Makefile.in
include $(MAKEFILE_IN)

APP      = disco

SRCEXT   = c
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)')
SRCDIRS := $(shell find . -name '*.$(SRCEXT)' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))

DEBUG    = -g -Wall
INCLUDES = -I$(H55)/include
#CFLAGS   = -O3 -c $(INCLUDES)
CFLAGS   = -c $(DEBUG) $(INCLUDES)
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
