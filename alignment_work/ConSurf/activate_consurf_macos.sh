#!/usr/bin/env bash

CONSURF_ROOT="/Users/natlinde/consurf_workspace"
BIOPERL_LIB="/opt/homebrew/opt/bioperl/libexec/lib/perl5"

export PATH="/opt/homebrew/opt/perl/bin:/opt/homebrew/bin:$PATH"
export PERL5LIB="$CONSURF_ROOT:$BIOPERL_LIB:$HOME/perl5/lib/perl5${PERL5LIB:+:$PERL5LIB}"
export CONSURFCONF="$CONSURF_ROOT/consurfrc.macos"

echo "ConSurf macOS environment active"
echo "Configuration: $CONSURFCONF"