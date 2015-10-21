#TODO: Check if Makefile.in exists
MAKEFILE_IN = $(PWD)/Makefile.in
include $(MAKEFILE_IN)

APP      = disco

SRCEXT   = c
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin
VISDIR   = python_plotting_scripts
PARDIR   = parfiles

SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)')
SRCDIRS := $(shell find . -name '*.$(SRCEXT)' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))

DEBUG    = -g -Wall
OPT      = -O3
INCLUDES = -I$(H55)/include
CFLAGS   = -c $(INCLUDES)
LDFLAGS  = -lm -lz -L$(H55)/lib -lhdf5

ifeq ($(strip $(USE_DEBUG)), 1)
	CFLAGS += $(DEBUG)
endif
ifeq ($(strip $(USE_OPT)), 1)
	CFLAGS += $(OPT)
endif

CFLAGS += $(LOCAL_C_FLAGS)
LDFLAGS += $(LOCAL_LD_FLAGS)

.PHONY: all clean distclean install

all: $(BINDIR)/$(APP)

install: $(BINDIR)/$(APP)
ifndef INSTALL_DIR
	$(error INSTALL_DIR has not been set in Makefile.in $(INSTALL_DIR))
endif
	@echo "Installing into $(INSTALL_DIR)..."
	@mkdir -p $(INSTALL_DIR)/$(BINDIR)
	@mkdir -p $(INSTALL_DIR)/$(VISDIR)
	@mkdir -p $(INSTALL_DIR)/$(PARDIR)
	@cp $(BINDIR)/$(APP) $(INSTALL_DIR)/$(BINDIR)/$(APP)
	@cp -r $(VISDIR)/* $(INSTALL_DIR)/$(VISDIR)/
	@cp -r $(PARDIR)/* $(INSTALL_DIR)/$(PARDIR)/

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
