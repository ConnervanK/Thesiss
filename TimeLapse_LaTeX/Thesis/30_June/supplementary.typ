// ============================================================
// Supplementary Material
// Subwavelength Imaging in Ground-Penetrating Radar Using Time-Lapse Data
// — MSc Thesis companion document
//
// Collects figures moved out of the main thesis chapters to keep those
// chapters concise. Each chapter here mirrors the section structure and
// numbering of its corresponding thesis chapter one-for-one; a pointer
// left behind in the main text (e.g. "§S1.1.1") locates the figure here.
//
// Compile with:  typst compile supplementary.typ
// ============================================================

#import "template.typ": *

#show: supplement.with(
  title:        "Supplementary Material",
  parent-title: "Subwavelength Imaging in Ground-Penetrating Radar — Using Time-Lapse Data",
  author:       "Conner van Kooten",
  date:         "August 3, 2026",
)

#include "supp_chapters/S1_methodology.typ"
#include "supp_chapters/S2_hypothesis1.typ"
#include "supp_chapters/S3_hypothesis2.typ"
#include "supp_chapters/S4_hypothesis3.typ"
