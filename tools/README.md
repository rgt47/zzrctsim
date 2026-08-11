# tools/

This directory holds the render-stamp helpers that zzcollab
vendors into every document-rendering compendium. They do two
things on every render. First, they stamp a discreet provenance
footer onto every page of the PDF, naming the source document and
the time of rendering. Second, they deposit a dated, git-versioned
copy of the PDF in the project's `share/` directory,
with a one-line entry added to a manifest. We have found that an
unlabeled PDF circulating among collaborators is a recurring
source of confusion; a footer that names its own source, and a
`share/` directory that accumulates an unambiguous record of what
was sent when, are the two halves of the remedy.

## What is here

- `stamp.tex` defines the footer: a small LaTeX preamble,
  included at render time, that defines three macros,
  `\stampsource`, `\stamptime`, and `\stampversion`. The default
  footer prints the first two; `\stampversion` is defined and
  available, but is left off the page because it already names the
  staged copy in `share/`. It carries `\providecommand` defaults
  so that a document still compiles when the values have not been
  supplied.
- `stamp-render.R` does the work. It renders the document,
  computes the provenance triple, injects the footer, copies the
  PDF into `share/`, and appends a row to the manifest. It
  handles `.Rmd`, `.qmd`, and `.md` input, driving rmarkdown,
  quarto, or pandoc as required.
- `render.sh` is a thin command-line wrapper over
  `stamp-render.R`.
- `README.md` is this file.

## The provenance triple

Both the footer and the staged filename derive from one triple:

- the **source document** (`report.Rmd`, `report.qmd`, ...);
- the **render time**, local;
- the **version**, from `git describe --tags --always
  --dirty=-wip`, which names the exact commit and marks an
  uncommitted tree with a `-wip` suffix.

The footer names the source document and the render time. The
staged filename names the artifact, as
`<prefix>-<date>-<time>-<version>.pdf`, where the prefix is the
document's basename, or its parent directory when the basename is
a generic stub such as `report`. Footer and filename therefore
record the same render, while each names the thing appropriate to
it.

## Where staged copies go

The staged copy is written to `analysis/report/share/` when an
`analysis/report/` tree exists, which is the zzcollab compendium
case, and to a repository-root `share/` otherwise. Each render
also appends a row to `share/MANIFEST.md`.

## Usage

For a one-off render of any supported document:

```
bash tools/render.sh analysis/report/report.Rmd
```

For R Markdown, the more convenient arrangement is to let the
document stamp itself on every render. Add the following to the
YAML header; the hook walks up to find `stamp-render.R`, so it
needs no package:

```
knit: (function(input, ...) { d <- dirname(input); while (!file.exists(file.path(d, 'tools', 'stamp-render.R')) && d != dirname(d)) d <- dirname(d); source(file.path(d, 'tools', 'stamp-render.R'))$value(input) })
```

For a multi-document Quarto project, the cleaner arrangement is a
`pre-render` script in `_quarto.yml`; the `render.sh` route
remains available for single-file Quarto renders.

## Provenance

These files are vendored from the zzcollab templates so that the
compendium builds standalone, without reaching outside the
repository. They are not project-specific; the canonical copies
live in zzcollab, and `zzc tools` reinstalls them. Project-
specific scripts may live alongside them here, but are written
per project and are not managed by zzcollab.
