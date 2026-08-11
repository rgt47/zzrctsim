#!/usr/bin/env bash
## tools/render.sh
##
## Render a document to a stamped PDF and stage a dated, versioned
## copy in the project's share/ directory. This is a thin wrapper
## over tools/stamp-render.R, which does the work for .Rmd, .qmd,
## and .md documents alike (it dispatches by extension internally
## and drives rmarkdown, quarto, or pandoc as required).
##
## Usage:
##   bash tools/render.sh <document.Rmd|.qmd|.md>
##
## Vendored by zzcollab; this file is not project-specific.

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <document.Rmd|.qmd|.md>" >&2
  exit 1
fi

input="$1"
if [[ ! -f "$input" ]]; then
  echo "error: file not found: $input" >&2
  exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
input_dir="$(cd "$(dirname "$input")" && pwd)"
abs_input="$input_dir/$(basename "$input")"

Rscript -e "source('$script_dir/stamp-render.R')\$value('$abs_input')"
