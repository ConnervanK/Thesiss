#import "../template.typ": *

#heading(level: 1, numbering: none, outlined: true)[Acronyms]

// Scoped locally so only the acronym column (marked *bold* below) is bold --
// the document-wide `show table.cell.where(y: 0): set text(weight: "bold")`
// (template.typ) is meant for tables with a real header row; this table has
// none, so without this override its row 0 (GPR) would get *both* columns
// bolded, including the spelled-out meaning.
#[
#show table.cell.where(y: 0): set text(weight: "regular")
#table(
  columns: (auto, 1fr),
  stroke: none,
  inset: (y: 0.2em, x: 0.4em),
  [*FDTD*],   [Finite-Difference Time-Domain],
  [*FFT*],    [Fast Fourier Transform],
  [*FWHM*],   [Full Width at Half Maximum],
  [*FWI*],    [Full-Waveform Inversion],
  [*GPR*],    [Ground-Penetrating Radar],
  [*HF*],     [High Frequency],
  [*LSRTM*],  [Least-Squares Reverse-Time Migration],
  [*MAE*],    [Mean Absolute Error],
  [*MUSIC*],  [Multiple Signal Classification],
  [*OLS*],    [Ordinary Least Squares],
  [*PCA*],    [Principal Component Analysis],
  [*PEC*],    [Perfect Electric Conductor],
  [*PML*],    [Perfectly Matched Layer],
  [*PSF*],    [Point-Spread Function],
  [*RANSAC*], [Random Sample Consensus],
  [*ROI*],    [Region of Interest],
  [*RTM*],    [Reverse-Time Migration],
  [*SNR*],    [Signal-to-Noise Ratio],
  [*STFT*],   [Short-Time Fourier Transform],
  [*SVD*],    [Singular Value Decomposition],
  [*UHF*],    [Ultra-High Frequency],
  [*VRP*],    [Vertical Radar Profile],
  [*WLS*],    [Weighted Least Squares],
)
]
