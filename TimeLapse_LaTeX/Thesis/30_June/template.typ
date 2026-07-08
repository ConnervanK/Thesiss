// ============================================================
// IDEA League / RWTH Applied Geophysics MSc Thesis Template
// Typst port of MScThesis.cls
// Compile with: typst compile main.typ
// (Tinymist VS Code extension provides live preview)
// ============================================================

// ---- Figure-path resolver ---------------------------------
// Figures live in TimeLapse_Figures/<study>/ at the repo root.
// 30_June/ is at TimeLapse_LaTeX/Thesis/30_June/, so three
// levels up reaches the repo root.
#let _fig_base = "../../../TimeLapse_Figures/"

#let img(name, width: 100%) = {
  let dir = if name.starts-with("VTL_") { "Vertical_TimeLapse_Study/" }
    else if name.starts-with("DTL_") { "Diagonal_TimeLapse_Study/" }
    else if name.starts-with("FF_")  { "FluidFlow_Study/" }
    else if name.starts-with("TLP_") { "TimeLapse_Processing/" }
    else if name.starts-with("TL_")  { "TimeLapse_Study/" }
    else if name.starts-with("RES_") { "Resolution_Study/" }
    else if name.starts-with("TLC_") { "TimeLapse_Cleaning/" }
    else if name.starts-with("FD_")  { "FieldData_Study/" }
    else                              { "" }
  image(_fig_base + dir + name, width: width)
}

// ---- Draft-note marker ------------------------------------
#let draftnote(note) = block(
  fill: rgb("#FFF0F0"),
  stroke: (left: 3pt + rgb("#8B0000")),
  inset: (left: 0.8em, right: 0.6em, top: 0.4em, bottom: 0.4em),
  below: 0.8em,
  above: 0.8em,
  text(fill: rgb("#8B0000"), weight: "bold", size: 10pt)[
    \[DRAFT --- #note\]
  ]
)

// ---- Unnumbered chapter (still appears in TOC) ------------
#let nonumchap(title) = heading(level: 1, numbering: none)[#title]

// ---- Math shorthand macros --------------------------------
// Mirror preamble.tex custom commands
#let Dz     = $Delta z$
#let Dx     = $Delta x$
#let Dt     = $Delta t$
#let Dphi   = $Delta phi.alt$
#let DPhi   = $Delta Phi$
#let Dtheta = $Delta theta$
#let kz     = $k_z$
#let kx     = $k_x$
#let vmig   = $v_"mig"$
#let eps    = $epsilon.alt$
#let XS     = $italic("XS")$

// ---- Two-column subfigure grid ----------------------------
// Usage: #subfigs(cols: 2, img("A.png"), img("B.png"))
#let subfigs(cols: 2, ..images) = grid(
  columns: (1fr,) * cols,
  gutter: 0.8em,
  ..images.pos()
)

// ---- Three-column subfigure grid --------------------------
#let subfigs3(..images) = grid(
  columns: (1fr, 1fr, 1fr),
  gutter: 0.6em,
  ..images.pos()
)

// ---- Four-column subfigure grid --------------------------
#let subfigs4(..images) = grid(
  columns: (1fr, 1fr, 1fr, 1fr),
  gutter: 0.5em,
  ..images.pos()
)

// ---- Two-row subfigure (2x2 grid) -------------------------
#let subfigs2x2(a, b, c, d) = grid(
  columns: (1fr, 1fr),
  rows: (auto, auto),
  gutter: 0.8em,
  a, b, c, d
)

// ---- Three-row subfigure (3x2 grid) -----------------------
#let subfigs3x2(a, b, c, d, e, f) = grid(
  columns: (1fr, 1fr),
  gutter: 0.8em,
  a, b, c, d, e, f
)

// ---- Six-panel (3x2) grid ---------------------------------
#let subfigs6(..images) = grid(
  columns: (1fr, 1fr),
  gutter: 0.8em,
  ..images.pos()
)

// ---- Inline paragraph heading (replaces \paragraph{}) ----
#let para-head(title) = [*#title* ]

// ---- Main thesis template ---------------------------------
#let thesis(
  title:          "Subwavelength Imaging in Ground-Penetrating Radar",
  subtitle:       "Using Time-Lapse Data",
  author:         "Conner van Kooten",
  date:           "June 30, 2026",
  supervisor-one: "",
  supervisor-two: "",
  keywords:       "GPR, time-lapse, phase-plane, migration, sub-wavelength",
  body
) = {

  // Document metadata
  set document(
    title:    title + " — " + subtitle,
    author:   author,
    keywords: keywords.split(", "),
  )

  // Page layout (a4wide equivalent)
  set page(
    paper: "a4",
    margin: (left: 3.0cm, right: 2.5cm, top: 2.8cm, bottom: 2.8cm),
    header-ascent: 40%,
    footer-descent: 30%,
  )

  // Body typography
  set text(font: "Linux Libertine O", size: 11pt, lang: "en", hyphenate: true)
  set par(justify: true, first-line-indent: 0pt, spacing: 0.65em)

  // Equation numbering
  set math.equation(numbering: "(1)", supplement: "Equation")

  // Figure settings
  set figure(gap: 0.5em, supplement: "Figure")
  set figure.caption(separator: [. ], position: bottom)

  // Caption style: small, sans-serif, bold label
  show figure.caption: c => {
    set text(size: 9pt, font: "Linux Biolinum O")
    [*#c.supplement #context c.counter.display(c.numbering):* #c.body]
  }

  // Table stroke defaults
  set table(stroke: none, inset: (x: 0.6em, y: 0.35em))
  show table.cell.where(y: 0): set text(weight: "bold")

  // ---- Heading numbering (chapter-section with dots) ------
  // Level 1 = Chapter, Level 2 = Section, Level 3 = Subsection
  set heading(numbering: "1.1.1.1")
  show heading.where(level: 1): set heading(supplement: "Chapter")
  show heading.where(level: 2): set heading(supplement: "Section")
  show heading.where(level: 3): set heading(supplement: "Section")
  show heading.where(level: 4): set heading(supplement: "Section")

  // ---- Chapter heading (fncychap-style) -------------------
  show heading.where(level: 1): h => {
    // Force new page for each chapter
    pagebreak(weak: true)
    v(1.0em)
    line(length: 100%, stroke: 0.8pt)
    v(0.25em)
    if h.numbering != none {
      context {
        let ch-num = counter(heading).display("1")
        text(font: "Linux Biolinum O", size: 13pt, weight: "regular")[
          CHAPTER #ch-num
        ]
      }
    }
    v(0.25em)
    line(length: 100%, stroke: 0.8pt)
    v(0.4em)
    align(
      right,
      text(font: "Linux Biolinum O", size: 21pt, weight: "bold")[#h.body]
    )
    v(2.2em)
  }

  // ---- Section heading ------------------------------------
  show heading.where(level: 2): h => {
    v(1.0em, weak: false)
    text(font: "Linux Biolinum O", size: 13pt, weight: "bold")[
      #if h.numbering != none {
        context counter(heading).display("1.1") + "  "
      }
      #h.body
    ]
    v(0.45em, weak: false)
  }

  // ---- Subsection heading ---------------------------------
  show heading.where(level: 3): h => {
    v(0.75em, weak: false)
    text(font: "Linux Biolinum O", size: 12pt, weight: "bold")[
      #if h.numbering != none {
        context counter(heading).display("1.1.1") + "  "
      }
      #h.body
    ]
    v(0.35em, weak: false)
  }

  // ---- Subsubsection / paragraph (unnumbered, inline) ----
  show heading.where(level: 4): h => {
    v(0.5em, weak: false)
    text(font: "Linux Biolinum O", size: 11pt, weight: "bold")[#h.body.]
    [ ]
  }

  // ---- Quote block ----------------------------------------
  set quote(block: true)
  show quote.where(block: true): q => block(
    inset: (left: 2em, right: 2em, top: 0.4em, bottom: 0.4em),
    above: 0.8em, below: 0.8em,
    q.body
  )

  // ============================================================
  // TITLE PAGES
  // ============================================================

  // Page 1: Main title page with logo
  set page(numbering: none, header: none, footer: none)

  align(center)[
    #image("IDEA_League_Logo.png", width: 8cm)
    #v(0.5em)
    #text(size: 16pt, weight: "bold")[#smallcaps[Master of Science in Applied Geophysics]]
    #v(0.2em)
    #text(size: 14pt, weight: "regular")[#smallcaps[Research Thesis]]
  ]

  v(1fr)
  line(length: 100%, stroke: 1.5pt)
  v(0.4em)
  align(right)[
    #text(size: 26pt, weight: "bold", font: "Linux Biolinum O")[#title] \
    #text(size: 16pt, weight: "bold", font: "Linux Biolinum O")[#subtitle] \
    #v(0.3em)
    #text(size: 14pt, weight: "bold")[#author]
  ]
  v(0.4em)
  line(length: 100%, stroke: 1.5pt)
  v(0.3em)
  text(size: 11pt)[#date]
  v(1fr)

  pagebreak()

  // Page 2: Second title page
  set page(numbering: none, header: none, footer: none)

  v(2cm)
  align(center)[
    #text(size: 24pt, weight: "bold", font: "Linux Biolinum O")[#title] \
    #text(size: 14pt, weight: "bold", font: "Linux Biolinum O")[#subtitle]
    #v(0.8em)
    #text(size: 14pt)[Master of Science Thesis]
    #v(1.5em)
    #text(size: 12pt)[for the degree of Master of Science in Applied Geophysics \ by]
    #v(1.5em)
    #text(size: 14pt, weight: "bold")[#author]
    #v(0.8em)
    #text(size: 12pt)[#date]
  ]
  v(1fr)

  pagebreak()

  // Page 3: Reader signature page
  set page(numbering: none, header: none, footer: none)

  v(3em)
  align(center)[
    #text(size: 13pt)[
      IDEA LEAGUE \
      JOINT MASTER'S IN APPLIED GEOPHYSICS
    ]
    #v(2em)
    Delft University of Technology, The Netherlands \
    ETH Zürich, Switzerland \
    RWTH Aachen, Germany
  ]
  v(3em)
  align(right)[Dated: _#date _]
  v(2em)
  grid(
    columns: (auto, 1fr),
    gutter: (1em, 1em),
    [Supervisor(s):], align(right, line(length: 7cm) + linebreak() + supervisor-one),
    [],             align(right, v(0.8em) + line(length: 7cm) + linebreak() + supervisor-two),
    v(1em),         [],
    [Committee Members:], align(right, line(length: 7cm) + linebreak() + supervisor-one),
    [],             align(right, v(0.8em) + line(length: 7cm) + linebreak() + supervisor-two),
  )
  v(1fr)

  pagebreak()

  // ============================================================
  // FRONT MATTER (roman page numbers)
  // ============================================================
  counter(page).update(1)
  set page(
    numbering: "i",
    header: context {
      let h1-hits = query(heading.where(level: 1).before(here()))
      let h1-title = if h1-hits.len() > 0 { h1-hits.last().body } else { [] }
      grid(
        columns: (1fr, auto),
        align(left,  text(weight: "bold", size: 9pt)[#h1-title]),
        align(right, text(weight: "bold", size: 9pt)[#counter(page).display("i")]),
      )
      line(length: 100%, stroke: 0.4pt)
    },
    footer: context {
      line(length: 100%, stroke: 0.4pt)
      v(0.2em)
      grid(
        columns: (1fr, auto),
        align(left,  text(size: 8pt)[#datetime.today().display("[day] [month repr:long] [year]")]),
        align(right, text(size: 8pt)[]),
      )
    }
  )

  body
}

// ---- Transition to main matter (call after frontmatter) ---
#let start-mainmatter() = {
  counter(page).update(1)
  counter(heading).update(0)
  set page(numbering: "1")
}

// ---- Transition to appendix ------------------------------
#let start-appendix() = {
  counter(heading).update(0)
  set heading(numbering: (..nums) => {
    let n = nums.pos()
    let letters = ("A","B","C","D","E","F","G","H","I","J")
    if n.len() == 1 {
      letters.at(n.at(0) - 1)
    } else {
      letters.at(n.at(0) - 1) + "." + n.slice(1).map(str).join(".")
    }
  })
}
