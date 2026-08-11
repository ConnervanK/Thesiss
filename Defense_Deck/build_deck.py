"""
build_deck.py
=============
Generate the baseline PowerPoint for the MSc defense.

    python build_deck.py
    python build_deck.py --out "Defense_2026.pptx"
    python build_deck.py --no-notes            # strip speaker notes

What you get
------------
18 content slides + 9 backup slides, 16:9, clean minimal white styling:

  * titles, the accent rule, and the sparse bullet text already typeset
  * every figure position marked by a labelled DROP ZONE giving the exact
    source path of the asset that belongs there
  * full speaker notes embedded in the notes pane of every slide
  * slides 5, 15, 16 and 18 built out properly rather than left as
    placeholders, because they are layout-heavy and fiddly to do by hand

Then you: open it, and for each drop zone, Insert > Picture (or Video) using
the path printed on the placeholder, and delete the grey box.

Design note
-----------
Nothing here uses PowerPoint's built-in layout placeholders. Everything is
positioned absolutely in inches on a blank layout, which is the only reliable
way to get identical geometry across PowerPoint versions and platforms.
"""

from __future__ import annotations

import argparse
import os

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Emu, Inches, Pt

# ==========================================================================
#  Style
# ==========================================================================

SLIDE_W, SLIDE_H = 13.333, 7.5          # inches, 16:9

DARK = RGBColor(0x16, 0x22, 0x3A)
ACCENT = RGBColor(0x00, 0x74, 0x8C)
WARM = RGBColor(0xC1, 0x44, 0x0E)
GREY = RGBColor(0x88, 0x92, 0xA6)
LIGHT = RGBColor(0xF4, 0xF6, 0xF9)
EDGE = RGBColor(0xC8, 0xD0, 0xDC)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
GOOD = RGBColor(0x1B, 0x7F, 0x4B)
BAD = RGBColor(0xB3, 0x26, 0x1E)

FONT = "Arial"

M = 0.62                                 # left/right margin
TITLE_TOP = 0.42
RULE_Y = 1.30
BODY_TOP = 1.55

# Where the animation renderer puts its output, relative to this file
ANIM_DIR = os.path.join("assets", "animations")
# Where your thesis figures live, relative to the Thesis root
FIG = os.path.join("..", "TimeLapse_Figures")
CURATED = os.path.join("..", "TimeLapse_LaTeX", "Thesis", "17_August")


# ==========================================================================
#  Primitives
# ==========================================================================

def blank(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


def textbox(slide, x, y, w, h, text, size=18, bold=False, color=DARK,
            align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, italic=False,
            space_after=6, line_spacing=1.0):
    """Multi-line textbox. `text` may be a str or a list of paragraph strings."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = 0
    tf.margin_top = tf.margin_bottom = 0

    lines = [text] if isinstance(text, str) else list(text)
    for i, line in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = line
        p.alignment = align
        p.space_after = Pt(space_after)
        p.line_spacing = line_spacing
        for r in p.runs:
            r.font.size = Pt(size)
            r.font.bold = bold
            r.font.italic = italic
            r.font.color.rgb = color
            r.font.name = FONT
    return tb


def title(slide, text, subtitle=None):
    """Slide title + thin teal rule. The deck's signature."""
    textbox(slide, M, TITLE_TOP, SLIDE_W - 2 * M, 0.75, text,
            size=32, bold=True, color=DARK)
    if subtitle:
        textbox(slide, M, TITLE_TOP + 0.62, SLIDE_W - 2 * M, 0.34, subtitle,
                size=14, color=GREY, italic=True)
    rule = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, Inches(M), Inches(RULE_Y),
        Inches(SLIDE_W - 2 * M), Pt(2.5))
    rule.fill.solid()
    rule.fill.fore_color.rgb = ACCENT
    rule.line.fill.background()
    rule.shadow.inherit = False


BULLET_STEP = 0.53          # vertical pitch between bullets, inches
BULLET_H = 0.42
SAFE_BOTTOM = 7.18          # nothing may extend past this


def bullets(slide, x, y, w, items, size=20):
    """
    Sparse bullets: a small teal square + text. Max 3, by house rule.

    The assertions are deliberate. Silent bottom-edge overflow is the single
    most common way a generated deck looks broken on a projector, and it is
    invisible in a thumbnail. Fail loudly at build time instead.
    """
    assert len(items) <= 3, "House rule: no more than three bullets per slide."
    last = y + (len(items) - 1) * BULLET_STEP + BULLET_H
    assert last <= SAFE_BOTTOM, (
        f"Bullets would overflow the slide: last bullet ends at {last:.2f}\", "
        f"limit is {SAFE_BOTTOM}\". Move them up or use fewer.")

    for i, item in enumerate(items):
        yy = y + i * BULLET_STEP
        sq = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x),
                                    Inches(yy + 0.12), Pt(9), Pt(9))
        sq.fill.solid()
        sq.fill.fore_color.rgb = ACCENT
        sq.line.fill.background()
        sq.shadow.inherit = False
        textbox(slide, x + 0.26, yy, w - 0.26, BULLET_H, item,
                size=size, color=DARK)


def dropzone(slide, x, y, w, h, label, path, kind="FIGURE"):
    """
    A labelled placeholder telling you exactly what to drop in and from where.

    kind: FIGURE | ANIMATION | DRAW
      FIGURE    -> Insert > Pictures
      ANIMATION -> Insert > Video > This Device  (set to Start: On Click)
      DRAW      -> build it by hand in PowerPoint; there is no source file
    """
    accent = {"FIGURE": ACCENT, "ANIMATION": WARM, "DRAW": GREY}[kind]

    box = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                 Inches(x), Inches(y), Inches(w), Inches(h))
    box.fill.solid()
    box.fill.fore_color.rgb = LIGHT
    box.line.color.rgb = accent
    box.line.width = Pt(1.5)
    box.line.dash_style = 4                      # dashed
    box.shadow.inherit = False
    box.text_frame.text = ""

    tag = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x + 0.14),
                                 Inches(y + 0.12), Inches(1.32), Inches(0.26))
    tag.fill.solid()
    tag.fill.fore_color.rgb = accent
    tag.line.fill.background()
    tag.shadow.inherit = False
    tf = tag.text_frame
    tf.text = kind
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    p.runs[0].font.size = Pt(10)
    p.runs[0].font.bold = True
    p.runs[0].font.color.rgb = WHITE
    p.runs[0].font.name = FONT

    body = [label]
    if path:
        body += ["", path]
    else:
        body += ["", "(no source file — build this one by hand)"]

    tb = slide.shapes.add_textbox(Inches(x + 0.2), Inches(y + h / 2 - 0.55),
                                  Inches(w - 0.4), Inches(1.1))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    for i, line in enumerate(body):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = line
        p.alignment = PP_ALIGN.CENTER
        for r in p.runs:
            r.font.name = "Consolas" if i == 2 else FONT
            r.font.size = Pt(13 if i == 0 else 9.5)
            r.font.bold = (i == 0)
            r.font.color.rgb = DARK if i == 0 else GREY
    return box


def pill(slide, x, y, w, h, text, fill, textcolor=WHITE, size=15, bold=True):
    """A rounded, filled label — used for verdicts and timeline stages."""
    s = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                               Inches(x), Inches(y), Inches(w), Inches(h))
    s.fill.solid()
    s.fill.fore_color.rgb = fill
    s.line.fill.background()
    s.shadow.inherit = False
    tf = s.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.text = text
    for p in tf.paragraphs:
        p.alignment = PP_ALIGN.CENTER
        for r in p.runs:
            r.font.size = Pt(size)
            r.font.bold = bold
            r.font.color.rgb = textcolor
            r.font.name = FONT
    return s


def panel(slide, x, y, w, h, fill=LIGHT, edge=EDGE):
    s = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                               Inches(x), Inches(y), Inches(w), Inches(h))
    s.fill.solid()
    s.fill.fore_color.rgb = fill
    s.line.color.rgb = edge
    s.line.width = Pt(1.0)
    s.shadow.inherit = False
    return s


def notes(slide, text, enabled=True):
    if enabled:
        slide.notes_slide.notes_text_frame.text = text.strip()


FOOTER_Y = 7.22


def slide_number(slide, n, total):
    textbox(slide, SLIDE_W - M - 1.2, FOOTER_Y, 1.2, 0.24,
            f"{n} / {total}", size=10, color=GREY, align=PP_ALIGN.RIGHT)


def timing(slide, secs):
    """Small timing hint, bottom-left. Delete these before presenting."""
    textbox(slide, M, FOOTER_Y, 2.0, 0.24, f"~{secs}s",
            size=10, color=GREY, italic=True)


# ==========================================================================
#  The deck
# ==========================================================================

def build(path: str, with_notes: bool = True, with_timing: bool = True):
    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W)
    prs.slide_height = Inches(SLIDE_H)

    TOTAL = 18
    n = [0]

    def new(t=None, sub=None, secs=None):
        n[0] += 1
        s = blank(prs)
        if t:
            title(s, t, sub)
        slide_number(s, n[0], TOTAL)
        if secs and with_timing:
            timing(s, secs)
        return s

    # ---------------------------------------------------------------- 1
    s = new(secs=30)
    textbox(s, M, 2.35, SLIDE_W - 2 * M, 1.2,
            "Seeing Movement Smaller Than the Wave",
            size=44, bold=True, color=DARK, align=PP_ALIGN.CENTER)
    r = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(4.6), Inches(3.62),
                           Inches(4.13), Pt(3))
    r.fill.solid(); r.fill.fore_color.rgb = ACCENT
    r.line.fill.background(); r.shadow.inherit = False
    textbox(s, M, 3.92, SLIDE_W - 2 * M, 0.5,
            "Phase-based sub-wavelength displacement sensing in time-lapse GPR",
            size=19, color=ACCENT, align=PP_ALIGN.CENTER)
    textbox(s, M, 5.25, SLIDE_W - 2 * M, 1.2,
            ["Conner  ·  MSc Thesis Defense",
             "Supervisors: [names]   ·   [institution]   ·   [date]"],
            size=14, color=GREY, align=PP_ALIGN.CENTER)
    notes(s, """
Good morning. Over the next twenty minutes I want to convince you of one thing:
a ground-penetrating radar can measure a movement that is far smaller than the
wave it uses to look. About thirty times smaller. And the reason it can is that
we have been reading the wrong half of the signal.
""", with_notes)

    # ---------------------------------------------------------------- 2
    s = new("Two Surveys, One Question", secs=55)
    dropzone(s, M, BODY_TOP + 0.15, 5.2, 4.0,
             "Cartoon: MONDAY — ground cross-section, fracture, fluid front", None, "DRAW")
    dropzone(s, M + 5.6, BODY_TOP + 0.15, 5.2, 4.0,
             "Cartoon: FRIDAY — identical, front advanced a few mm", None, "DRAW")
    pill(s, 5.95, BODY_TOP + 1.85, 0.75, 0.65, "?", WARM, size=26)
    bullets(s, M, 6.15, 11.0, ["Monday → Friday", "What moved? How far?"], size=21)
    notes(s, """
Here is the problem that motivates everything. You survey a site. You come back
later and survey it again. Something in between has moved — a fracture has
filled with water, a slope has crept, a wetting front has advanced. You want to
know what moved and by how much.

That is a TIME-LAPSE question, and I will come back to it properly in a few
minutes. But to answer it, we first have to be honest about a much older
limitation — one that applies even to a single, static survey. So let me PARK
time-lapse for now and start there.

[Coaching] Say "park it" out loud. It signals the structure and buys patience
through slides 3 and 4.
""", with_notes)

    # ---------------------------------------------------------------- 3
    s = new("Resolution Is Set by the Wavelength",
            "One survey. One snapshot. The classical limit.", secs=50)
    dropzone(s, M, BODY_TOP + 0.1, 11.0, 4.15,
             "ANIMATION 1 — two scatterers converging  (let it run to ~½λ, pause here)",
             os.path.join(ANIM_DIR, "anim1_imaging_floor.mp4"), "ANIMATION")
    bullets(s, M, 6.15, 11.0,
            ["1.5 GHz in ice → λ ≈ 113 mm",
             "Two reflectors closer than ~λ/2 merge"], size=21)
    notes(s, """
GPR builds its image from a short electromagnetic pulse. Everything about how
finely it can see is set by the length of that pulse. In my setup — 1.5 GHz in
ice — one wavelength is about eleven centimetres.

[The animation runs at λ = 1 m for readability. Say so: "I've scaled the
wavelength to one metre here so the numbers are easy to read — the physics is
scale-invariant."]

Two reflectors closer than roughly half a wavelength stop being two reflectors.
Laterally that limit comes from the Fresnel zone; vertically from the pulse
bandwidth. Same physics as the Rayleigh criterion in optics.

And crucially: MIGRATION DOES NOT FIX THIS. Migration moves recorded energy back
to where it belongs — it focuses — but you still read the result off as an
AMPLITUDE picture, and that is still bound by the same floor.
""", with_notes)

    # ---------------------------------------------------------------- 4
    s = new("Where the Image Gives Up",
            "Three migration algorithms, all hitting the same wall", secs=60)
    dropzone(s, M, BODY_TOP + 0.1, 7.6, 4.35,
             "Measured PSF collapse: 1λ → ½λ → ¼λ → ⅛λ, three methods",
             os.path.join(CURATED, "fig_res_psf_crop.png"), "FIGURE")
    panel(s, M + 8.0, BODY_TOP + 0.1, 3.0, 4.35)
    textbox(s, M + 8.25, BODY_TOP + 0.42, 2.5, 0.4,
            "Rayleigh ratio", size=16, bold=True, color=DARK)
    rows = [("1 λ", "1.77", True), ("½ λ", "0.89", False),
            ("¼ λ", "0.44", False), ("⅛ λ", "0.22", False)]
    for i, (lab, val, ok) in enumerate(rows):
        yy = BODY_TOP + 1.02 + i * 0.68
        textbox(s, M + 8.25, yy, 0.9, 0.4, lab, size=17, bold=True, color=DARK)
        pill(s, M + 9.25, yy - 0.05, 1.5, 0.44, val, GOOD if ok else BAD, size=15)
    textbox(s, M + 8.25, BODY_TOP + 3.80, 2.5, 0.5,
            "Kirchhoff. Below 1 = merged.", size=11, color=GREY, italic=True)
    bullets(s, M, 6.15, 7.4,
            ["Two static scatterers, closing in", "Below ~½λ: one lobe"], size=19)
    notes(s, """
This is that limit measured rather than asserted. Two stationary scatterers, gap
swept from one wavelength down to an eighth, migrated three different ways —
Kirchhoff, Gazdag phase-shift, and time-reversal back-propagation.

Read it top to bottom. At one wavelength, all three give two clean peaks. At a
half, Kirchhoff has already lost it — ratio 0.89, below the Rayleigh threshold
of one. By a quarter wavelength everything is a single lobe. There is genuinely
no information left in the amplitude.

Back-propagation is the best of the three — it holds on to about a quarter
wavelength — but the point is not which algorithm wins. They ALL hit a wall, and
the wall is at roughly half a wavelength. Call this the IMAGING RESOLUTION FLOOR.

[Coaching] Physically point at the ½λ Kirchhoff panel. That single panel is the
whole argument.
""", with_notes)

    # ---------------------------------------------------------------- 5
    s = new("Imaging Is Not Monitoring",
            "The distinction the rest of this talk depends on", secs=60)
    LX, RX, PW = M, M + 5.85, 5.3
    TY, PH = BODY_TOP + 0.2, 4.5
    panel(s, LX, TY, PW, PH, fill=LIGHT, edge=EDGE)
    panel(s, RX, TY, PW, PH, fill=RGBColor(0xE8, 0xF3, 0xF6), edge=ACCENT)
    div = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(M + 5.55),
                             Inches(TY), Pt(2.5), Inches(PH))
    div.fill.solid(); div.fill.fore_color.rgb = GREY
    div.line.fill.background(); div.shadow.inherit = False

    pill(s, LX + 0.3, TY + 0.28, PW - 0.6, 0.55, "IMAGING", GREY, size=19)
    pill(s, RX + 0.3, TY + 0.28, PW - 0.6, 0.55, "MONITORING", ACCENT, size=19)

    left = [("Question", "Where is it?"), ("Data", "One survey"),
            ("Axis", "Space"), ("Limit", "Rayleigh,  ~λ/2")]
    right = [("Question", "How far did it move?"), ("Data", "Two surveys"),
             ("Axis", "Time"), ("Limit", "?")]
    for i, ((k, v), (k2, v2)) in enumerate(zip(left, right)):
        yy = TY + 1.12 + i * 0.82
        textbox(s, LX + 0.32, yy, 1.5, 0.3, k, size=11, color=GREY)
        textbox(s, LX + 0.32, yy + 0.26, PW - 0.64, 0.42, v,
                size=19, bold=True, color=DARK)
        textbox(s, RX + 0.32, yy, 1.5, 0.3, k2, size=11, color=GREY)
        textbox(s, RX + 0.32, yy + 0.26, PW - 0.64, 0.42, v2,
                size=19, bold=True,
                color=WARM if v2 == "?" else DARK)
    textbox(s, RX + 0.32, TY + PH - 0.62, PW - 0.64, 0.42,
            "← the subject of this thesis", size=13, italic=True, color=WARM)
    textbox(s, M, 6.6, 11.0, 0.45,
            "Same radar. Same wavelength. Different question — so, possibly, "
            "a different limit.", size=17, color=DARK, align=PP_ALIGN.CENTER)
    notes(s, """
★ THE PIVOT SLIDE. Slow down. If they remember one slide, make it this one.

Now I want to draw a distinction that the rest of the talk depends on, so let me
be very explicit about it.

Everything so far has been IMAGING RESOLUTION. One survey, one snapshot, and the
question is: where is this thing and can I separate it from its neighbour? That
is a question about SPACE, and the answer is the half-wavelength floor we just
measured.

What I actually care about is a different question. Two surveys, and the
question is: how far did this thing move between them? That is a question about
CHANGE OVER TIME. Call it MONITORING RESOLUTION.

These are not the same quantity, and there is no law that says they must have
the same limit. That is the opening this thesis walks through.

[Build order: left panel → divider → right panel header → the "?" ]
""", with_notes)

    # ---------------------------------------------------------------- 6
    s = new("Subtracting Two Images Doesn't Help",
            "Now it IS time-lapse — and it inherits the same floor", secs=35)
    dropzone(s, M, BODY_TOP + 0.1, 11.0, 4.15,
             "ANIMATION 2 — baseline / monitor / difference, fixed colour scale",
             os.path.join(ANIM_DIR, "anim2_monitoring_floor.mp4"), "ANIMATION")
    bullets(s, M, 6.15, 11.0,
            ["Difference the two images", "Still stuck at ~½λ"], size=21)
    notes(s, """
The obvious thing to try: migrate both surveys, subtract, look at what is left.
That works — as long as the movement is big.

But the difference image inherits the same floor. Below about a quarter
wavelength the baseline and monitor lobes overlap so completely that their
difference carries no usable position information.

[Point out the FIXED colour scale — this is the honest bit. Auto-scaling a
difference image always shows SOMETHING, which is exactly how people fool
themselves. Pinned to the baseline's own peak, the difference visibly dies.]

So amplitude differencing does not buy a new limit. It re-imports the old one.
We need to stop looking at amplitude.
""", with_notes)

    # ---------------------------------------------------------------- 7
    s = new("Amplitude Reads in Steps. Phase Reads Continuously.",
            "The key insight — one scatterer, one trace", secs=60)
    dropzone(s, M, BODY_TOP + 0.1, 11.0, 4.15,
             "ANIMATION 3 — trace + phasor + staircase-vs-straight-line",
             os.path.join(ANIM_DIR, "anim3_envelope_vs_phase.mp4"), "ANIMATION")
    bullets(s, M, 6.15, 11.0,
            ["Envelope: quantised, broad", "Phase: continuous, exact"], size=21)
    notes(s, """
Here is the physical insight the whole method rests on.

A migrated point scatterer is an ENVELOPE — a broad amplitude lobe, roughly a
wavelength wide — with a wave oscillating underneath it. To read a position off
the envelope you have to locate the peak of that broad lobe on a discrete grid.
Watch panel (c): the grey curve is a STAIRCASE. Below one sample spacing it
reports zero movement.

Now the phase of the wave underneath. It advances continuously and linearly with
position. No threshold, no floor, no minimum detectable shift. Panel (b) makes
it concrete: the displacement is literally an angle.

The information was always there. Amplitude imaging quantises it away.

[Coaching] Best moment in the talk. Let the animation run once in silence before
narrating the second half.

[If challenged: "doesn't the envelope shift too?" — Yes, exactly, and that is
why I am careful to say READS IN STEPS rather than DOESN'T MOVE. The envelope
translates; what fails is separating two of them (slides 4 and 6) and reading a
sub-sample position off a discretised amplitude image.]
""", with_notes)

    # ---------------------------------------------------------------- 8
    s = new("A Shift in Space Is a Tilt in Wavenumber",
            "Fourier shift theorem: translation rotates phase, "
            "leaves amplitude untouched", secs=50)
    dropzone(s, M, BODY_TOP + 0.1, 11.0, 4.15,
             "ANIMATION 4a — image sliding | cross-spectrum phase tilting | live fit",
             os.path.join(ANIM_DIR, "anim4_kdomain_lateral.mp4"), "ANIMATION")
    bullets(s, M, 6.15, 11.0, ["Move in space", "Tilt in wavenumber"], size=21)
    notes(s, """
To turn that insight into a number, I move into the two-dimensional Fourier
domain — wavenumber space.

The Fourier shift theorem says something very convenient: translating an image
in space does not change its amplitude spectrum AT ALL. It only rotates the
phase, by an amount proportional to the wavenumber.

So on the left, the scatterer slides. On the right, watch the phase go from flat
— no shift — to a smooth ramp. The steeper that ramp, the bigger the shift. The
slope IS the displacement.

That is the conversion I need: from a displacement I cannot see, to a slope I
can measure.
""", with_notes)

    # ---------------------------------------------------------------- 9
    s = new("One Multiplication Cleans It Up",
            "The cross-spectrum: everything the two images share cancels",
            secs=95)
    # Four formula/clip pairs. Reveal the formula, click, play the clip.
    steps = [
        ("1", "m(x,z) = b(x−Δx, z−Δz)",
         "anim5_act1_model.mp4", "scatterer moves"),
        ("2", "(same relation, after migration)",
         "anim5_act2_migrated.mp4", "→ migrated"),
        ("3", "B(k) = ∫∫ b e^(−j k·r) dr        M = B · e^(−j k·Δ)",
         "anim5_act3_spectrum.mp4", "→ wavenumber   |M| = |B|"),
        ("4", "XS = B · M*  =  |B|² e^(+j k·Δ)        ∠XS = Φ = kₓΔx + k_zΔz",
         "anim5_act4_cross.mp4", "→ |XS| and ∠XS"),
    ]
    ty = BODY_TOP + 0.12
    for i, (num, formula, clip, caption) in enumerate(steps):
        yy = ty + i * 1.10
        panel(s, M, yy, 11.0, 0.95)
        pill(s, M + 0.16, yy + 0.20, 0.55, 0.55, num, ACCENT, size=17)
        textbox(s, M + 0.92, yy + 0.13, 6.3, 0.32, formula,
                size=13.5, bold=True, color=DARK)
        textbox(s, M + 0.92, yy + 0.52, 6.3, 0.30, caption,
                size=11, color=GREY, italic=True)
        box = slide.shapes if False else s.shapes
        tag = box.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(M + 7.45),
                            Inches(yy + 0.20), Inches(3.35), Inches(0.55))
        tag.fill.solid()
        tag.fill.fore_color.rgb = LIGHT
        tag.line.color.rgb = WARM
        tag.line.width = Pt(1.25)
        tag.line.dash_style = 4
        tag.shadow.inherit = False
        tf = tag.text_frame
        tf.word_wrap = True
        tf.vertical_anchor = MSO_ANCHOR.MIDDLE
        tf.text = os.path.join(ANIM_DIR, clip)
        pp = tf.paragraphs[0]
        pp.alignment = PP_ALIGN.CENTER
        pp.runs[0].font.size = Pt(9)
        pp.runs[0].font.name = "Consolas"
        pp.runs[0].font.color.rgb = DARK
    textbox(s, M, 6.30, 11.0, 0.4,
            "Build order: reveal the formula → click → the clip morphs the "
            "panels into the next domain.",
            size=12, color=GREY, italic=True, align=PP_ALIGN.CENTER)
    textbox(s, M, 6.70, 11.0, 0.42,
            "Φ(k_z , k_x)   =   k_z · Δz   +   k_x · Δx",
            size=24, bold=True, color=ACCENT, align=PP_ALIGN.CENTER)
    notes(s, """
THE DERIVATION SLIDE. Four formula/clip pairs; you control every step.

Two panels stay in the same place the whole way through — left is always the
baseline, right is always the monitor. Only what they REPRESENT changes. Say
that at the start; it is what stops the domain change from being disorienting.

[1] The complication first: the two surveys differ by a translation, and that
is all. m(x,z) = b(x−Δx, z−Δz). Click: the right-hand scatterer walks to its
monitor position.

[2] Click: both panels become their migrated versions. The sharp point becomes
a broad lobe — that lobe is the resolution floor from earlier. But the relation
still holds exactly, because migration is linear and I use the same operator on
both surveys.

[3] Now into the two-dimensional Fourier domain. The shift theorem says a
translation does not touch the amplitude spectrum at all; it only multiplies by
a phase ramp. Click: both panels become their spectra — and they are
IDENTICAL. Pause here. That is the point: a displacement is completely
invisible in the amplitude spectrum. Everything we need is hiding in the phase.

[4] So form the cross-spectrum: baseline times the complex conjugate of the
monitor. Click. The magnitude term is real and positive, so taking the phase
kills it completely, and what survives is a flat plane through the origin. Its
slope along k_z is the vertical displacement; its slope along k_x is the
lateral one.

That's it — that's the method. Everything after this is measurement.

[Coaching] Say "that's the method" out loud. It signals the conceptual work is
done and the audience can relax.

[If the projector will not play video: the four *_final.png stills in
assets/animations/ are the end state of each act. Step through those instead.]
""", with_notes)

    # ---------------------------------------------------------------- 10
    s = new("Lateral Tilts One Way. Vertical Tilts the Other.",
            "Two independent slopes of the same plane", secs=55)
    dropzone(s, M, BODY_TOP + 0.1, 11.0, 3.80,
             "ANIMATION 4b — 2×2: lateral vs vertical, image and k-domain",
             os.path.join(ANIM_DIR, "anim4_kdomain_both.mp4"), "ANIMATION")
    bullets(s, M, 5.58, 11.0,
            ["Lateral → slope in k_x", "Vertical → slope in k_z",
             "Diagonal → both, recovered separately"], size=19)
    notes(s, """
The two directions do not interfere. Purely lateral motion tilts the plane along
k_x and leaves k_z flat. Purely vertical does the reverse. A diagonal movement
tilts both, and the fit recovers each component independently — I tested that on
a 2:1 diagonal path and got the ratio back correctly.

One asymmetry worth flagging, because it surprised me. In the SPATIAL domain
these two directions behave completely differently — a migrated image separates
lateral spatial frequencies across its width, like a prism, but it does not do
the same in depth. In the FOURIER domain that asymmetry disappears: both are
just slopes of the same plane.

That is the practical argument for doing this in wavenumber space rather than
reading phase gradients off the image directly.
""", with_notes)

    # ---------------------------------------------------------------- 11
    s = new("Fit the Plane, Weight by Energy",
            "Thousands of noisy observations of one plane", secs=65)
    dropzone(s, M, BODY_TOP + 0.1, 11.0, 3.45,
             "Four-panel conceptual workflow: baseline → monitor → phase → fit",
             os.path.join(FIG, "Resolution_Study", "General",
                          "RES_018_Phase-Plane_Method__Conceptual_Workflow.png"),
             "FIGURE")
    bullets(s, M, 5.35, 7.3,
            ["Weighted least squares", "Bright bins count more",
             "Two slopes → Δx, Δz"], size=19)
    panel(s, M + 7.9, 5.35, 3.1, 1.55, fill=LIGHT)
    textbox(s, M + 8.05, 5.52, 2.8, 0.35, "equation 2 of 2", size=10,
            color=GREY, italic=True)
    textbox(s, M + 8.05, 5.88, 2.8, 0.6,
            "u = (AᵀW²A)⁻¹ AᵀW² Φ", size=16, bold=True, color=DARK)
    notes(s, """
EQUATION 2 OF 2. Do not derive it. Point at it and move on.

A real image gives me thousands of frequency bins, each an independent noisy
observation of that plane, so this is an over-determined least-squares problem.

Two things make it robust. First, a band-pass mask that keeps the fit inside the
coherent part of the wavelet spectrum, so the total phase rotation never wraps
past ±π. Second — and this matters more — I weight every bin by the
cross-spectrum energy. A bin at the antenna's peak power is far more trustworthy
than one at the edge of the band, and unweighted least squares treats them
identically. In the fit panel you can see it working: bright high-weight points
sit tight on the line, dim ones scatter.

On this clean example the recovered shift is ten millimetres against a true ten
millimetres.
""", with_notes)

    # ---------------------------------------------------------------- 12
    s = new("Phase Takes Over Exactly Where Amplitude Quits",
            "Four movement types · three migration methods · clean data", secs=75)
    dropzone(s, M, BODY_TOP + 0.05, 11.0, 4.30,
             "Detectability map — top: Rayleigh ratio | bottom: phase error %",
             os.path.join(FIG, "Hypothesis_1", "General",
                          "H1_037_Hypothesis_1_--_Detectability_Map_Amplitude_vs_Phase.png"),
             "FIGURE")
    bullets(s, M, 6.15, 6.6,
            ["Top: amplitude dies below ½λ",
             "Bottom: phase accurate to ¹⁄₃₂λ"], size=18)
    pill(s, M + 7.0, 6.18, 4.0, 0.88,
         "3–4 orders of magnitude,\nin the regime amplitude cannot reach",
         ACCENT, size=14)
    notes(s, """
★ THE MONEY SLIDE. Do not rush. Pause for two full seconds after the first line.

This is the central result, and I want to give you a moment to read it.

[STOP TALKING. Two seconds.]

Four movement types across the columns — lateral, vertical, diagonal, and a
graded fluid front. Three migration methods in colour. Read right to left:
displacement gets SMALLER as you move left.

Top row is amplitude — the Rayleigh ratio. Above the dotted line is resolvable.
Every curve crosses below it around a half to a quarter wavelength.

Bottom row is the phase method, percentage error on a log scale. Notice WHERE it
drops: at or BEFORE the point where amplitude fails. In the lateral case Gazdag
sits at a few thousandths of a percent error — a few hundredths of a millimetre
— at a displacement of four millimetres, one thirty-second of a wavelength.

The two rows are mirror images of each other, and that is the whole claim: phase
becomes reliable precisely where amplitude stops being. Not by a little. By
three to four orders of magnitude.
""", with_notes)

    # ---------------------------------------------------------------- 13
    s = new("It Works on a Fluid Front Too",
            "A harder target than a metal cylinder — deliberately", secs=50)
    dropzone(s, M, BODY_TOP + 0.1, 5.3, 3.55,
             "Graded wetting zone: 7 permittivity steps, water → air",
             os.path.join(FIG, "Hypothesis_1", "General",
                          "H1_028_FluidFlow_--_Model_Set_Up_baseline_graded_zone_true_scale.png"),
             "FIGURE")
    dropzone(s, M + 5.7, BODY_TOP + 0.1, 5.3, 3.55,
             "MAE summary — grey out all groups except FluidFlow",
             os.path.join(FIG, "Hypothesis_1", "General",
                          "H1_039_Hypothesis_1_--_MAE_Summary_Across_Movement_Types.png"),
             "FIGURE")
    pill(s, 4.2, 5.42, 4.9, 0.6, "≈ 1 mm error on a 113 mm target",
         ACCENT, size=17)
    bullets(s, M, 6.15, 11.0,
            ["Amplitude fails even earlier than for a hard reflector",
             "Phase: still ~1 mm"], size=18)
    notes(s, """
[CUT THIS SLIDE FIRST if you are running long — slide 12 already shows the
FluidFlow column.]

One fair objection to everything so far is that a perfect metal cylinder is a
very generous target. So I repeated the whole sweep with something closer to
what I actually care about: a graded wetting front — seven permittivity steps
from water down to air, standing in for fluid advancing through a fracture.

This target is HARDER. Its point-spread function is intrinsically broader, so
amplitude differencing fails even sooner — Gazdag's ratio is already below one
at two full wavelengths.

The phase fit does not care. About a millimetre of error, on a target whose own
smear is a hundred and thirteen millimetres wide. So the result is not an
artefact of the idealised point-scatterer geometry.
""", with_notes)

    # ---------------------------------------------------------------- 14
    s = new("What Migration Does to Nothing at All",
            "Laplace noise, fitted to the real field dataset", secs=65)
    dropzone(s, M, BODY_TOP + 0.1, 11.0, 3.60,
             "Pure noise in → Kirchhoff (coherent!) vs Gazdag (speckle)",
             os.path.join(FIG, "Hypothesis_2", "Noisy", "Kirchhoff",
                          "H2_003_Migrating_Pure_Noise_no_scatterers_no_signal_"
                          "--_Kirchhoff_vs.png"),
             "FIGURE")
    textbox(s, M, 5.26, 11.0, 0.32,
            "Circle one bright cluster in the Kirchhoff panel and label it: "
            "“looks like a target — isn’t one.”",
            size=11, color=WARM, italic=True)
    bullets(s, M, 5.62, 11.0,
            ["Input: pure noise, no target",
             "Kirchhoff manufactures structure", "Gazdag stays speckle"], size=18)
    notes(s, """
Real data is not clean, so I repeated all four experiments under noise — and not
arbitrary noise: a Laplace distribution fitted to the amplitude statistics of
the real field dataset, heavier-tailed than a Gaussian because field noise is
dominated by clutter spikes rather than thermal noise.

Before testing signal-plus-noise, I ran a sanity check I would recommend to
anyone doing this: migrate noise ALONE. No scatterer, no signal.

Kirchhoff turns it into that. Coherent, streaky, focused-looking structure —
manufactured from nothing. And that is structural, not bad luck: delay-and-sum
stacks many traces along travel-time curves, and stacking imposes coherence on
incoherent input BY CONSTRUCTION. Gazdag, which downward-continues in the
frequency domain, leaves it as speckle.

This matters because my pipeline crops its fitting window around the brightest
region. If the brightest region is a noise artefact, the fit is meaningless.
""", with_notes)

    # ---------------------------------------------------------------- 15
    s = new("Accuracy and Safety Pick Different Winners",
            "Mean absolute error, sub-½λ regime, all four movement types", secs=75)
    hdr_y = BODY_TOP + 0.18
    textbox(s, M + 3.55, hdr_y, 3.0, 0.4, "ACCURACY (MAE)", size=13,
            bold=True, color=GREY, align=PP_ALIGN.CENTER)
    textbox(s, M + 7.0, hdr_y, 3.6, 0.4, "FALSE-POSITIVE RISK", size=13,
            bold=True, color=GREY, align=PP_ALIGN.CENTER)
    meth = [("Kirchhoff", "4.5 mm", GOOD, "Invents coherent targets", BAD),
            ("Back-prop  (sign-bit)", "11.1 mm", GREY, "Stays incoherent", GOOD),
            ("Gazdag", "20.3 mm", BAD, "Stays incoherent", GOOD)]
    for i, (name, mae, mc, risk, rc) in enumerate(meth):
        yy = hdr_y + 0.62 + i * 1.02
        panel(s, M, yy, 11.0, 0.82)
        textbox(s, M + 0.25, yy + 0.2, 3.2, 0.45, name, size=18, bold=True,
                color=DARK)
        pill(s, M + 3.75, yy + 0.16, 2.6, 0.5, mae, mc, size=17)
        pill(s, M + 7.2, yy + 0.16, 3.2, 0.5, risk, rc, size=14)
    pill(s, 3.1, 5.42, 7.1, 0.62, "Most accurate  ≠  safest", DARK, size=20)
    bullets(s, M, 6.15, 11.0,
            ["Kirchhoff: accurate, but cross-check it",
             "Sign-bit back-prop: the conservative choice"], size=17)
    notes(s, """
So which migration should you use? I framed this originally as "which is most
noise-robust", and the honest answer is that this was the WRONG QUESTION,
because the two sensible definitions of robust give opposite answers.

On raw accuracy, Kirchhoff wins clearly — 4.5 mm mean absolute error across all
four movement types, about 2.5× better than back-propagation and 4.5× better
than Gazdag.

On false-positive risk, Kirchhoff is the ONLY method that fails. It is the one
that manufactures coherent structure from noise.

Back-propagation trades a factor of two in accuracy for that safety — but only
once you make it noise-robust, which took one specific change: injecting the
SIGN of the time-reversed wavefield rather than its amplitude. Focusing is
driven by zero-crossings, so the sign carries all the information that matters,
and it clamps every noise spike to the same ±1 as the real signal.

Gazdag loses on both counts.

My recommendation is conditional: Kirchhoff when you need accuracy and can
cross-check; sign-bit back-propagation when a false detection is the expensive
outcome.

[Caveat if pushed: this is one noise level — 10% of each B-scan's own signal
standard deviation. A sweep across levels is future work.]
""", with_notes)

    # ---------------------------------------------------------------- 16
    s = new("A Borehole, and 38 Surveys",
            "Push–pull tracer test in fractured rock  ·  6 June 2016", secs=55)
    dropzone(s, M, BODY_TOP + 0.1, 4.0, 4.0,
             "Borehole domain: water-filled channel, injection at 78.7 m",
             os.path.join(FIG, "FieldData_Study", "Compilations",
                          "FD_borehole_domain_schematic.png"), "FIGURE")
    stages = [("PUSH", "1 → 3", "injecting", "up"),
              ("CHASE", "3 → 8", "still injecting", "up"),
              ("WAIT", "8 → 20", "pumping stops", "settle"),
              ("PULL", "20 → 38", "extracting", "back")]
    sx, sw = M + 4.5, 1.52
    for i, (nm, pair, what, expect) in enumerate(stages):
        xx = sx + i * (sw + 0.22)
        pill(s, xx, BODY_TOP + 0.55, sw, 0.55, nm,
             ACCENT if i < 2 else GREY, size=15)
        textbox(s, xx, BODY_TOP + 1.22, sw, 0.32, pair, size=13, bold=True,
                color=DARK, align=PP_ALIGN.CENTER)
        textbox(s, xx, BODY_TOP + 1.58, sw, 0.5, what, size=10.5, color=GREY,
                align=PP_ALIGN.CENTER)
        pill(s, xx, BODY_TOP + 2.15, sw, 0.42, expect, LIGHT, DARK, size=12)
        if i < 3:
            textbox(s, xx + sw, BODY_TOP + 0.62, 0.22, 0.4, "→", size=17,
                    color=GREY, align=PP_ALIGN.CENTER)
    textbox(s, sx, BODY_TOP + 2.90, 6.5, 0.9,
            "No ground truth. Target geometry unknown, true displacement "
            "unknown,\nregion of influence painted by hand in napari.",
            size=13, color=WARM, italic=True)
    bullets(s, M, 6.15, 11.0,
            ["38 borehole profiles, fractured rock",
             "Injection at 78.7 m — same band as the tracked reflector"], size=17)
    notes(s, """
Finally, real data. A single-hole borehole GPR dataset from a controlled
push-pull tracer test in fractured rock — 38 profiles down the same borehole,
from the Shakas campaign.

Four stages, each with a physical expectation. Push: fluid actively injected, at
78.7 metres depth. Chase: continued injection at a lower rate. Wait: pumping
stops. Pull: fluid extracted.

The critical difference from everything before this is that here I have NO
GROUND TRUTH. I do not know the target geometry, I do not know the true
displacement, and there is no pre-labelled region to fit. I picked the region of
influence by hand, in napari, painting directly on the cross-spectrum — and I
will come back to how much that choice matters.
""", with_notes)

    # ---------------------------------------------------------------- 17
    s = new("Two Independent Pipelines, One Trajectory",
            "Cumulative displacement along the borehole", secs=80)
    dropzone(s, M, BODY_TOP + 0.1, 8.0, 3.80,
             "Stage-anchored displacement: Gazdag (solid) vs back-prop (dashed)",
             os.path.join(CURATED, "fig_fielddata_stages.png"), "FIGURE")
    panel(s, M + 8.4, BODY_TOP + 0.1, 2.6, 3.80)
    textbox(s, M + 8.6, BODY_TOP + 0.35, 2.2, 0.4, "Push stage", size=13,
            bold=True, color=GREY)
    textbox(s, M + 8.6, BODY_TOP + 0.75, 2.2, 1.1,
            ["−1.66 m  Gazdag", "−2.06 m  back-prop"], size=15, bold=True,
            color=DARK)
    textbox(s, M + 8.6, BODY_TOP + 1.92, 2.2, 0.34, "Agreement", size=13,
            bold=True, color=GREY)
    pill(s, M + 8.6, BODY_TOP + 2.28, 2.2, 0.75,
         "same sign,\nevery stage", GOOD, size=14)
    textbox(s, M + 8.6, BODY_TOP + 3.10, 2.2, 0.7,
            "No shared processing downstream of the raw B-scans.",
            size=11, color=GREY, italic=True)
    bullets(s, M, 5.58, 11.0,
            ["Up while injecting", "Back down when pumping stops",
             "Two methods, no shared processing"], size=17)
    notes(s, """
Cumulative displacement along the borehole, Gazdag solid, back-propagation
dashed.

The front moves UPWARD while fluid is pumped in, through Push and Chase. Pumping
stops, and it settles back down through Wait and Pull. Exactly the physically
expected pattern, given injection at 78.7 metres — inside the same depth band as
the reflector I am tracking.

But the part I would actually defend is the AGREEMENT between the two lines.
Gazdag and back-propagation share NO processing steps downstream of the raw
B-scans — different migration physics, different gain corrections, different
artefact suppression. They agree on the sign of the displacement in every single
stage. Two independent pipelines converging on the same trajectory is much
stronger evidence that this is a real physical signal than either estimate alone.

I want to be precise about what this does and does not show. It shows the
pipeline GENERALISES — it survives real, noisy, geometry-unknown data and returns
something stable and physically sensible. It is not a sub-wavelength validation,
because I have no independent ground truth to validate against.

═══════════════════════════════════════════════════════════════════
★ REHEARSE THIS ANSWER VERBATIM — most likely technical challenge:

Q: "Your field wavelength is ~1 m and your Push displacement is 1.66 m —
    shouldn't the phase have wrapped?"

A: The hand-picked k-space ROI is a narrow band around the spectral peak, which
   keeps the total phase rotation inside ±π across that band. The ±λ/2
   unambiguous range quoted in the synthetic study reflects the much wider
   passband used there (|k| < 1.4 k_c). The genuinely sub-wavelength stage here
   is Pull, at 0.14 m against a ~1 m wavelength.
═══════════════════════════════════════════════════════════════════
""", with_notes)

    # ---------------------------------------------------------------- 18
    s = new("What This Buys, and What It Costs", secs=100)
    cols = [
        ("✔  CLAIM", GOOD,
         ["Sub-wavelength movement,",
          "read from phase.",
          "",
          "Down to ¹⁄₃₂ λ,",
          "to a few tenths of a mm —",
          "where amplitude has",
          "already collapsed.",
          "",
          "Holds for a fluid front,",
          "not just a point scatterer."]),
        ("⚠  CAVEAT", WARM,
         ["Migration choice under noise",
          "is a trade-off, not solved.",
          "One noise level tested.",
          "",
          "ROI picked by hand is the",
          "weakest link: RANSAC shifts",
          "a k-space pick ≤ 6%,",
          "a B-scan pick up to 59%.",
          "",
          "Velocity change would mimic",
          "a vertical shift. Untested."]),
        ("→  NEXT", ACCENT,
         ["Multiple simultaneous movers",
          "— each imprints its own ramp,",
          "so a global fit blends them.",
          "Needs a localised variant.",
          "",
          "Noise-level sweep to find the",
          "accuracy / safety crossover.",
          "",
          "All 38 field profiles, not 5.",
          "Step-by-step, not 4 stages."]),
    ]
    cw = 3.45
    for i, (hdr, col, lines) in enumerate(cols):
        xx = M + i * (cw + 0.42)
        pill(s, xx, BODY_TOP + 0.15, cw, 0.55, hdr, col, size=17)
        textbox(s, xx + 0.1, BODY_TOP + 0.92, cw - 0.2, 4.05, lines,
                size=13.5, color=DARK, space_after=1, line_spacing=1.06)
    textbox(s, M, 6.72, 11.0, 0.4, "Thank you.", size=22, bold=True,
            color=DARK, align=PP_ALIGN.CENTER)
    notes(s, """
To close.

THE CLAIM. Time-lapse GPR can recover a movement about thirty times smaller than
its own wavelength, by reading the phase of the migrated image instead of its
amplitude. Down to one thirty-second of a wavelength, to a few tenths of a
millimetre, exactly in the regime where the amplitude difference image has
already collapsed — and that holds for a fluid front as well as a point
scatterer.

THE CAVEATS. Three, honestly. One: under noise the choice of migration is a real
trade-off, not a solved problem, and I tested a single noise level. Two:
choosing the region of influence by hand is the weakest link in the field
pipeline — refitting with RANSAC changes a hand-picked k-space estimate by at
most six percent, but an amplitude-gated one by up to fifty-nine. Three: a
genuine change in medium velocity between surveys would masquerade as a vertical
shift, and I have not tested that.

WHAT'S NEXT. The most interesting one is multiple simultaneously-moving
scatterers — in the wavenumber domain each imprints its own tilted ramp, and a
global fit returns a meaningless blend, so that needs a localised variant.
Beyond that: a noise-level sweep to find where the accuracy and safety criteria
cross over, and running the corrected back-propagation over all thirty-eight
field profiles instead of five.

Thank you. I'm happy to take questions.
""", with_notes)

    # ================================================================
    #  Backup slides
    # ================================================================
    sep = blank(prs)
    textbox(sep, M, 3.2, SLIDE_W - 2 * M, 1.0, "Backup",
            size=40, bold=True, color=GREY, align=PP_ALIGN.CENTER)
    textbox(sep, M, 4.2, SLIDE_W - 2 * M, 0.5,
            "Do not present. Jump here from Q&A.",
            size=15, color=GREY, italic=True, align=PP_ALIGN.CENTER)

    backups = [
        ("Phase-Plane Fit, Full Diagnostic", "“Show me an actual fit.”",
         os.path.join(FIG, "Hypothesis_1", "Gazdag",
                      "H1_038_Phase-Plane_Workflow_--_Gazdag_Lateral_132_λ.png")),
        ("Does the Weighting Actually Help?", "OLS vs WLS, hardest noisy case",
         os.path.join(FIG, "Hypothesis_2", "General",
                      "H2_008_Chapter_54_--_OLS_vs_WLS_Phase-Plane_Fitting_"
                      "on_Noisy_Data.png")),
        ("Why Laplace?", "Field noise distribution + Gaussian comparison",
         os.path.join(FIG, "Hypothesis_2", "General",
                      "H2_002_Noise_distribution_at_7_Constant_Velocity_"
                      "pre-gain_with_Lapl.png")),
        ("How Much Does the ROI Matter?", "WLS vs RANSAC, all pairs and methods",
         os.path.join(FIG, "FieldData_Study", "Compilations",
                      "FD_ransac_vs_wls_displacement_summary.png")),
        ("Does Slide 12 Survive Noise?", "Noisy detectability map",
         os.path.join(FIG, "Hypothesis_2", "General",
                      "H2_029_Hypothesis_2_--_Detectability_Map_Noisy.png")),
        ("Clean vs Noisy, Every Pair", "Dumbbell plot, log axis",
         os.path.join(FIG, "Hypothesis_2", "General",
                      "H2_030_Hypothesis_2_--_Clean_vs_Noisy_MAE_Dumbbell_Plot.png")),
        ("What Does the Raw Field Data Look Like?",
         "5 profiles × 4 representations",
         os.path.join(FIG, "FieldData_Study", "Compilations",
                      "FD_profile_migration_grid.png")),
        ("What Is Sign-Bit Time Reversal?",
         "Peak-normalised vs sign-bit excitation",
         os.path.join(FIG, "Hypothesis_2", "Back-Propagation",
                      "H2_007_Sign-Bit_Time-Reversed_Excitation_--_"
                      "B-scans_and_Spectra_Lat.png")),
        ("Vertical Motion in the k-Domain", "Animation 4, vertical mode",
         os.path.join(ANIM_DIR, "anim4_kdomain_vertical.mp4")),
    ]
    for i, (t, sub, p) in enumerate(backups):
        s = blank(prs)
        title(s, t, sub)
        kind = "ANIMATION" if p.endswith(".mp4") else "FIGURE"
        dropzone(s, M, BODY_TOP + 0.15, 11.0, 5.35, "Backup asset", p, kind)
        textbox(s, SLIDE_W - M - 1.4, FOOTER_Y, 1.4, 0.24, f"B{i+1}",
                size=10, color=GREY, align=PP_ALIGN.RIGHT)

    prs.save(path)
    return prs


# ==========================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="Defense_Deck.pptx")
    ap.add_argument("--no-notes", action="store_true")
    ap.add_argument("--no-timing", action="store_true",
                    help="omit the small ~Ns timing hints")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out = args.out if os.path.isabs(args.out) else os.path.join(here, args.out)

    prs = build(out, with_notes=not args.no_notes, with_timing=not args.no_timing)
    n = len(prs.slides._sldIdLst)
    print(f"Wrote {out}")
    print(f"  {n} slides  (18 content + 1 divider + 9 backup)")
    print(f"  16:9, {SLIDE_W}\" x {SLIDE_H}\"")
    print()
    print("Next steps:")
    print("  1. python animations/render_all.py")
    print("  2. Open the deck. For each DROP ZONE: Insert > Picture / Video,")
    print("     using the path printed on the grey box, then delete the box.")
    print("  3. Set every video to Start: On Click  (Playback tab).")
    print("  4. Slides 2, 5, 15, 16, 18 need no images — 5/15/16/18 are already")
    print("     built out; only slide 2's two cartoons are yours to draw.")


if __name__ == "__main__":
    main()
