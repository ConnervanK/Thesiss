// ============================================================
// Sensing the Invisible: Time-lapse Imaging of Sub-wavelength Processes using Ground Penetrating Radar
// Using Time-Lapse Data — MSc Thesis
//
// Typst port of TimeLapse_LaTeX/Thesis/27_June/main.tex
// Built on the IDEA League / RWTH Applied Geophysics template
//
// Compile with:  typst compile main.typ
// Or open in VS Code with the Tinymist extension for live preview.
// ============================================================

#import "template.typ": *

#show: thesis.with(
  title:          "Sensing the Invisible",
  subtitle:       "Time-lapse Imaging of Sub-wavelength Processes using Ground Penetrating Radar",
  subtitle-display: [Time-lapse Imaging of Sub-wavelength Processes #linebreak() using Ground Penetrating Radar],
  author:         "Conner Marcus van Kooten",
  date:           "August 3, 2026",
  supervisor-one: "Dr. Alexis Shakas",
  supervisor-two: "Dr. Johannes Aichele",
  committee-extra: "Dr. Joeri Brackenhoff",
  keywords:       "GPR, time-lapse, phase-plane, migration, sub-wavelength",
)

// ======================== FRONT MATTER ========================

#include "chapters/00_abstract.typ"
#pagebreak()

#include "frontmatter/acknowledgements.typ"
#pagebreak()

// Table of Contents
#outline(
  title: [Table of Contents],
  indent: 2em,
  depth: 3,
)
#pagebreak()

// List of Figures
#outline(
  title: [List of Figures],
  target: figure.where(kind: image),
)
#pagebreak()

// List of Tables
#outline( 
  title: [List of Tables],
  target: figure.where(kind: table),
)
#pagebreak()

#include "frontmatter/acronyms.typ"
#pagebreak() 

// ======================== MAIN MATTER =========================

// Reset page counter and chapter counter for mainmatter
#counter(page).update(1)
#counter(heading).update(0)
#set page(numbering: "1")

#include "chapters/01_introduction.typ"
#include "chapters/02_literature_review.typ"
#include "chapters/03_theory.typ"
#include "chapters/04_methodology.typ"
#include "chapters/05_hypothesis1_translation.typ"
#include "chapters/06_hypothesis2_migration_noise.typ"
#include "chapters/07_hypothesis3_complex.typ"
#include "chapters/08_discussion.typ"

// Summary (unnumbered backmatter chapter)
#include "backmatter/summary.typ"

// Bibliography
#pagebreak()
#heading(level: 1, numbering: none, outlined: true)[Bibliography]
#bibliography("references.bib", style: "chicago-author-date", title: none)

