// ============================================================
// MSc Thesis Defense Slides
// "Phase-Based Sub-Wavelength Displacement Sensing in Time-Lapse GPR"
// Built from thesis_defense_outline.md against TimeLapse_LaTeX/Thesis/30_June
//
// Compile (from the repo root, so relative figure paths resolve):
//   tinymist compile TimeLapse_LaTeX/Thesis/17_August/defense_slides.typ out.pdf
// Or open in VS Code with the Tinymist extension for live preview.
//
// Page is Typst's "presentation-16-9" preset = 29.7 x 16.7 cm.
// ============================================================

#let fig-base  = "../../../TimeLapse_Figures/"
#let logo-path = "../30_June/IDEA_League_Logo.png"

#let cDark        = rgb("#16223A")
#let cAccent       = rgb("#00748C")
#let cAccentLight = rgb("#DCEEF1")

#set page(
  paper: "presentation-16-9",
  margin: (x: 1.8cm, top: 1.3cm, bottom: 1.2cm),
  fill: white,
  numbering: "1",
)
#set text(font: "Arial", size: 16pt, fill: cDark, lang: "en")
// Spacing is controlled explicitly via #v() throughout this deck, so
// paragraph auto-spacing is zeroed out to avoid double-spacing stacked
// on top of those manual gaps.
#set par(leading: 0.5em, spacing: 0pt, justify: false)
#set list(marker: [•], spacing: 0.5em, indent: 0.3em)
#set strong(delta: 0)

// ---------------------------------------------------------
// Building blocks
// ---------------------------------------------------------

// Claim banner: the one-sentence takeaway at the top of every content slide.
#let claim(body) = block(
  width: 100%,
  fill: cAccentLight,
  inset: (x: 0.75em, y: 0.5em),
  radius: 3pt,
  below: 0.5em,
)[#text(weight: "bold", size: 15pt)[#body]]

// Dark title bar at the top of every content slide.
#let slide-title(body) = block(
  width: 100%,
  fill: cDark,
  inset: (x: 0.8em, y: 0.35em),
  below: 0.45em,
)[#text(fill: white, weight: "bold", size: 21pt)[#body]]

// One slide == one page.
#let slide(title: none, body) = {
  if title != none { slide-title(title) }
  body
  pagebreak(weak: true)
}

// ============================================================
// TITLE SLIDE
// ============================================================
#align(center + horizon)[
  #image(logo-path, height: 1.5cm)

  #v(1em)
  #text(size: 25pt, weight: "bold", fill: cDark)[
    Phase-Based Sub-Wavelength Displacement
  ] \
  #text(size: 25pt, weight: "bold", fill: cDark)[
    Sensing in Time-Lapse GPR
  ]

  #v(0.6em)
  #text(size: 15pt)[MSc Thesis Defense]

  #v(0.3em)
  #text(size: 11.5pt, style: "italic")[
    Sensing the Invisible: Time-lapse Imaging of
  ] \
  #text(size: 11.5pt, style: "italic")[
    Sub-wavelength Processes using Ground Penetrating Radar
  ]

  #v(1.1em)
  #text(size: 14pt)[Conner Marcus van Kooten]

  #v(0.7em)
  #text(size: 11.5pt)[
    Supervisors: Dr. Alexis Shakas, Dr. Johannes Aichele
  ] \
  #text(size: 11.5pt)[
    Committee: Dr. Joeri Brackenhoff
  ]

  #v(1em)
  #text(size: 10.5pt, fill: cDark.lighten(30%))[17 August 2026]
]
#pagebreak(weak: true)

// ============================================================
// SLIDE 1 --- Resolution Problems in GPR
// ============================================================
#slide(title: [Resolution Problems in GPR])[
  #claim[GPR has a hard, wavelength-scale resolution floor --- and the exact processes worth monitoring over time sit below it.]

  #grid(
    columns: (48%, 1fr),
    gutter: 1cm,
    align(center + top)[
      #image("fig_res_psf_crop.png", height: 10.6cm)
    ],
    align(top)[
      - Two reflectors closer than $tilde.op 1/2 lambda$ can't be told apart in a migrated amplitude image
      - Confirmed for all three algorithms: Kirchhoff fails at $1 lambda$, Gazdag at $1/2 lambda$, back-propagation (best) at $1/4 lambda$
      - Migration relocates energy --- it does not shrink the floor
    ],
  )
]
// Say: why the floor exists physically (Fresnel zone / pulse bandwidth), then
// let the merging PSFs make the point visually. Land hard on: the changes we
// care about -- a fracture aperture filling with fluid, mm-scale creep -- are
// structurally invisible here. This is the foundation the rest of the talk
// stands on; don't rush it.

// ============================================================
// SLIDE 2 --- Looking at Phase
// ============================================================
#slide(title: [Looking at Phase])[
  #claim[A sub-wavelength shift is invisible in amplitude, but it produces an exact, linear phase ramp --- recoverable to a fraction of a millimetre, precisely where amplitude has already failed.]

  #align(center)[
    #image(
      fig-base + "Hypothesis_1/General/H1_037_Hypothesis_1_--_Detectability_Map_Amplitude_vs_Phase.png",
      height: 7.6cm,
    )
  ]

  #v(0.3em)
  - Compare baseline vs. monitor in the 2D Fourier domain $arrow$ the shift becomes a straight phase ramp; its slope *is* the displacement
  - Recovered to $<0.35$ mm, down to $1/32 lambda$ ($tilde.op 3.5$ mm) --- orders of magnitude below the wavelength
  - Holds for a rigid point target *and* a diffuse, graded fluid front --- not a special case
]
// Say: this is the core result -- treat it as the centrepiece, not one of
// several findings. Explain the phase-ramp idea in one or two sentences (skip
// the WLS/masking machinery unless asked), then let the crossover figure
// carry the weight: phase overtakes amplitude exactly where amplitude gives
// out. This slide can reasonably take the most time of the five.

// ============================================================
// SLIDE 3 --- Adding Noise
// ============================================================
#slide(title: [Adding Noise])[
  #claim[Realistic noise doesn't break the method --- but it makes migration choice matter, and "best" depends on which failure mode you're trying to avoid.]

  #align(center)[
    #image(
      fig-base + "Hypothesis_2/Noisy/Kirchhoff/H2_003_Migrating_Pure_Noise_no_scatterers_no_signal_--_Kirchhoff_Ga.png",
      width: 100%,
    )
  ]

  #v(0.3em)
  - Kirchhoff: most accurate overall (4.5 mm mean error) --- but turns pure noise into a false, target-like focus
  - Back-propagation (sign-bit time-reversal): less accurate (11.1 mm) --- never produces that false focus
  - No unconditionally "best" algorithm --- accuracy and false-positive avoidance trade off against each other
]
// Say: open with the pure-noise image -- it's the single most convincing
// piece of evidence you have, and it needs no numbers to land. Then give the
// two headline error figures. Frame it explicitly as a trade-off decision,
// not a limitation.

// ============================================================
// SLIDE 4 --- Field Data
// ============================================================
#slide(title: [Field Data])[
  #claim[Applied blind --- unknown geometry, unknown true displacement --- to a real fluid-injection experiment, the method recovers a trajectory that tracks the known pumping schedule, confirmed by two independent processing chains.]

  #grid(
    columns: (52%, 1fr),
    gutter: 1cm,
    align(center + top)[
      #image("fig_fielddata_stages.png", height: 10.6cm)
    ],
    align(top)[
      - Real borehole GPR, fluid injected at 78.7 m depth --- no ground truth available
      - Recovered motion: upward $tilde.op$ 1.7--2 m while injecting, reversing downward once pumping stops --- matches the injection schedule
      - Gazdag and back-propagation share no downstream processing, yet agree on direction at *every* stage
    ],
  )
]
// Say: frame this as "does it survive contact with reality." The
// independent-method agreement is your strongest piece of evidence it's a
// real signal, not an artefact -- say that explicitly, don't leave it
// implicit. Leave the lateral-shift/fracture-orientation open question for
// Q&A rather than putting it on the slide.

// ============================================================
// SLIDE 5 --- Conclusion
// ============================================================
#slide(title: [Conclusion: It Works, and It's Robust])[
  #claim[Phase, not amplitude, delivers real super-resolution GPR monitoring --- demonstrated clean, stress-tested under noise, and confirmed on real field data.]

  #v(0.3em)
  - Sub-wavelength displacement is recoverable from a standard survey --- no new hardware required
  - Robust under realistic noise, with a clear, decision-ready choice between accuracy and false-positive avoidance
  - Confirmed --- procedurally and physically --- against real, uncontrolled field data

  #v(1.4em)
  #block(
    width: 100%,
    fill: cAccentLight,
    inset: (x: 0.9em, y: 0.6em),
    radius: 3pt,
  )[
    #text(size: 13pt, style: "italic")[
      Can time-lapse GPR accurately track sub-wavelength movement --- a form
      of super-resolution monitoring --- by analysing phase rather than
      amplitude?
    ]

    #v(0.35em)
    #text(size: 17pt, weight: "bold")[Yes.]
  ]
]
// Say: close by restating the research question and answering it in one
// sentence: yes. Keep future work / limitations off the slide entirely --
// hold those for questions, where they'll land as evidence of rigor rather
// than diluting the story.

// ============================================================
// CLOSING SLIDE
// ============================================================
#align(center + horizon)[
  #text(size: 27pt, weight: "bold", fill: cDark)[Thank you]

  #v(0.4em)
  #text(size: 15pt, fill: cDark.lighten(20%))[Questions?]
]
