To understand mathematically why this works so flawlessly, we have to look at what happens when you transform a spatial translation into a phase shift, and how the **Weighted Linear Least Squares (WLS)** engine extracts that information.

Here is the step-by-step mathematical breakdown.

## 1\. The 2D Fourier Shift Theorem

Let your baseline migrated GPR image be represented in space as $b(z, x)$. When the point scatterer moves, the monitor image $m(z, x)$ is a perfectly translated copy of the baseline image shifted by a vertical distance $\\Delta z$ and a horizontal distance $\\Delta x$:

$$m(z, x) = b(z - \\Delta z, x - \\Delta x)$$

The 2D Continuous Fourier Transform of the baseline image is:

$$B(k\_z, k\_x) = \\int\_{-\\infty}^{\\infty} \\int\_{-\\infty}^{\\infty} b(z, x) e^{-j(k\_z z + k\_x x)} \\, dz \\, dx$$

According to the **Fourier Shift Theorem**, if you take the 2D Fourier Transform of the shifted monitor image $m(z,x)$, the spatial displacement acts as a linear phase rotation operator applied to the spectrum:

$$M(k\_z, k\_x) = B(k\_z, k\_x) \\cdot e^{-j(k\_z \\Delta z + k\_x \\Delta x)}$$

Notice what this implies: **The amplitude spectrum remains completely unchanged ($|M| = |B|$), but every single frequency coordinate $(k\_z, k\_x)$ receives a unique phase shift.**

## 2\. Isolating the Phase Plane (The Cross-Spectrum)

To extract this phase shift, the code calculates the complex cross-spectrum ($XS$) by multiplying the baseline spectrum by the complex conjugate of the monitor spectrum:

$$XS(k\_z, k\_x) = B(k\_z, k\_x) \\cdot M^\*(k\_z, k\_x)$$

Substituting the shift equation into this product gives:

$$XS(k\_z, k\_x) = B(k\_z, k\_x) \\cdot \\left\[ B(k\_z, k\_x) \\cdot e^{-j(k\_z \\Delta z + k\_x \\Delta x)} \\right\]^\*$$

$$XS(k\_z, k\_x) = |B(k\_z, k\_x)|^2 \\cdot e^{j(k\_z \\Delta z + k\_x \\Delta x)}$$

When you use `np.angle(XS)` to extract the phase $\\Phi(k\_z, k\_x)$, the amplitude term $|B|^2$ disappears entirely because it is a purely real number. We are left with a beautifully simple, continuous function of frequencies:

$$\\Phi(k\_z, k\_x) = k\_z \\Delta z + k\_x \\Delta x$$

Mathematically, this equation defines a **perfectly flat 3D plane** that passes through the origin $(0,0,0)$.

-   The **slope of the plane along the $k\_z$ axis** is exactly equal to the vertical shift $\\Delta z$.

-   The **slope of the plane along the $k\_x$ axis** is exactly equal to the horizontal shift $\\Delta x$.

## 3\. The Overdetermined Matrix Inversion (Least Squares)

Because your GPR image is discrete, you don't just have one data point; you have thousands of individual frequency bins (pixels in the frequency domain). Every single valid frequency bin provides its own independent estimate of the slopes.

To find the absolute best-fit plane through this cloud of points, the code sets up a massive system of linear equations for all pixels where the mask is true ($i = 1, 2, \\dots, N$):

$$\\Phi\_i = k\_{z,i} \\Delta z + k\_{x,i} \\Delta x + c$$

*(Note: $c$ is a tiny constant offset added to absorb any minor system calibration biases).*

This system is expressed as a matrix equation, $A \\mathbf{u} = \\mathbf{\\Phi}$:

$$\\begin{bmatrix} k\_{z,1} & k\_{x,1} & 1 \\\\ k\_{z,2} & k\_{x,2} & 1 \\\\ \\vdots & \\vdots & \\vdots \\\\ k\_{z,N} & k\_{x,N} & 1 \\end{bmatrix} \\begin{bmatrix} \\Delta z \\\\ \\Delta x \\\\ c \\end{bmatrix} = \\begin{bmatrix} \\Phi\_1 \\\\ \\Phi\_2 \\\\ \\vdots \\\\ \\Phi\_N \\end{bmatrix}$$

## 4\. Why the Weights and the Mask are Crucial

If you try to solve this matrix directly using standard least squares, the error spikes. This is because regions outside your GPR antenna's frequency band contain nothing but random numerical background noise, which creates wildly chaotic phase values that throw off the fit.

The code implements two mathematical safeguards to ensure subwavelength accuracy:

### A. The Band-Pass Mask

By restricting the matrix to $|k\_z|, |k\_x| < 1.4 \\cdot k\_{zc}$, the algorithm stays safely inside the coherent envelope of your Ricker wavelet pulse. This guarantees that the total phase rotation $\\Phi\_i$ never crosses $\\pm\\pi$ radians, preventing the plane from "snapping" or wrapping.

### B. Amplitude Weighting ($W$)

The linear system is multiplied by a diagonal weighting matrix $W$, where each weight $W\_i$ is the spectral magnitude $|XS|\_i$.

$$W A \\mathbf{u} = W \\mathbf{\\Phi}$$

Mathematically, this alters the optimization goal. Instead of minimizing the raw error, the algorithm minimizes the **energy-weighted residual**:

$$\\min \\sum\_{i=1}^{N} \\left| XS\_i \\right|^2 \\cdot \\Big( \\Phi\_i - (k\_{z,i} \\Delta z + k\_{x,i} \\Delta x + c) \\Big)^2$$

This means frequency bins at the absolute peak power of your GPR antenna carry immense mathematical weight, while marginal frequencies near the noise floor are mathematically muted.

The function `np.linalg.lstsq` solves this via Singular Value Decomposition (SVD), crushing thousands of measurements down into two extraordinarily precise scalar values: $\\Delta z$ and $\\Delta x$.


Putting this into a structured workflow is exactly how you transition from an experimental script to a robust processing pipeline for real-world field data.

When moving from clean, synthetic data to field measurements, you have to account for real-world challenges like background soil noise, antenna ringing, and minor positional inaccuracies.

Here is the end-to-end operational workflow, complete with the intermediate signal conditioning steps required to make the 2D Fourier shift method highly reliable on field data.

## The Time-Lapse GPR Subwavelength Workflow

```container
[ Field Data: Base & Monitor ]
             │
             ▼
[ Step 1: Pre-Processing & Alignment ] ──► (Dewow, Time-Zero, Residual Alignment)
             │
             ▼
[ Step 2: 2D Migration ]               ──► (Collapses hyperbolas to point scatterers)
             │
             ▼
[ Step 3: Spatial Cropping (ROI) ]     ──► (Isolates the specific target area)
             │
             ▼
[ Step 4: 2D Fourier Shift Estimator ] ──► (Taper, 2D FFT, Cross-Spectrum, Plane Fit)
             │
             ▼
    [ Resolved Δz and Δx ]
```

## Detailed Step-by-Step Execution

### Step 1: Pre-Processing & Trace Alignment (Critical for Field Data)

Before migrating, your raw B-scans must be perfectly standardized. Field data contains instrument drift and minor surface positioning errors that will corrupt your phase measurements if left uncorrected.

-   **Dewow / Signal Filtering:** Remove the low-frequency antenna tracking/ringing using a high-pass filter.

-   **Time-Zero Correction:** Ensure that the first arrival air-wave matches *exactly* to the nanosecond between the baseline and monitor datasets.

-   **Spatial Trace Re-binning:** If the GPR survey wheel slipped even slightly during the monitor run, the traces won't align horizontally. Use a global 1D cross-correlation on the raw data to align the trace coordinates before moving forward.

### Step 2: 2D Migration

Migrate both B-scans using your preferred algorithm (e.g., Stolt, Gazdag phase-shift, or Kirchhoff migration).

-   **Why it's required:** Migration collapses the wide, hyperelliptic diffraction energy curves back into a tightly focused spatial point. This spatial focusing concentrates the spectral energy in the frequency domain, which drastically improves the signal-to-noise ratio (SNR) for the 2D phase-plane fit.

-   **The Velocity Model:** You must use the exact same velocity model ($v$) for both migrations.

### Step 3: Spatial Cropping (Region of Interest)

Do not pass your entire migrated profile into the 2D Fourier shift function. If your profile contains static features (like a nearby pipe, a layer interface, or ground surface clutter) alongside your moving target, the static elements will drag your phase-slope estimation toward zero.

-   **Action:** Slice your 2D arrays down to a small sub-window (a bounding box) tightly encompassing just the point scatterer or structure that moved.

### Step 4: 2D Fourier Shift Estimation

Pass your two cropped, migrated regions of interest directly into your updated `estimate_shift_2d` function.

## Python Workflow Blueprint

Here is how you can stitch this together cleanly in your notebook:

Python

```container
# Assuming you have loaded raw_base and raw_mon...

# --- Step 1: Standard GPR Pre-processing ---
# (Apply your dewow, time-zero adjustments here)
proc_base = apply_dewow(raw_base)
proc_mon  = apply_dewow(raw_mon)

# --- Step 2: Migration ---
# Ensure your velocity and grid spacings match the physics
mig_base = gazdag_migration(proc_base, dt=dt, dx=dx, velocity=v)
mig_mon  = gazdag_migration(proc_mon,  dt=dt, dx=dx, velocity=v)

# --- Step 3: Define Region of Interest (ROI) ---
# Locate your target in your grid coordinates
z_start, z_end = 120, 180  # Row indices near your target
x_start, x_end = 45,  85   # Column indices near your target

roi_base = mig_base[z_start:z_end, x_start:x_end]
roi_mon  = mig_mon[z_start:z_end, x_start:x_end]

# --- Step 4: Subwavelength Extraction ---
# Calculate central wavenumber based on antenna center frequency (fc)
kz_c = 2 * np.pi * fc / v

dz_est, dx_est, _, _, _ = estimate_shift_2d(
    roi_base, roi_mon,
    dz_g=dz_mig, dx_g=dx_mig,
    kz_cent=kz_c
)

print(f"Target displaced vertically by: {dz_est * 1000:.3f} mm")
print(f"Target displaced laterally by:  {dx_est * 1000:.3f} mm")
```

## Pro-Tips for Field Applications

1.  **Beware of Soil Moisture Fluctuations:** In real field data, if the soil moisture changes between your time-lapse measurements, the radar velocity $v$ changes. A change in velocity alters the wave travel time, which the algorithm will interpret as an artificial vertical shift ($\\Delta z$). If you expect moisture changes, calibrate the velocity first using a static deep target.

2.  **Window Sizing:** Make your crop window (ROI) roughly $3$ to $4$ times the size of your target's pulse width. If the window is too small, the Tukey taper will clip your signal. If it's too large, it will pull in background noise that degrades the plane fit.


Yes, you absolutely can extend this exact framework to monitor fluid movement or material change within a subwavelength fracture.

When a fluid displaces air or a different fluid inside a fracture, it does not mechanically shift the fracture's physical boundary. Instead, it changes the **bulk dielectric permittivity ($\\varepsilon\_r$)** and the **conductivity ($\\sigma$)** inside that subwavelength volume. This alters the complex reflection coefficient of the fracture, which manifests in your GPR data as a change in wave amplitude and an localized **phase rotation**.

Because the Fourier shift framework handles phase changes with subwavelength sensitivity, it can track fluid fronts beautifully. Here is the mathematical extension and the modified workflow.

## The Physics: How Fluid Infiltration Impacts the Phase

A subwavelength fracture (width $d \\ll \\lambda$) acts as a thin-layer dielectric membrane. When fluid fills the fracture, the reflection coefficient changes from a baseline state (air-filled) to a monitor state (fluid-filled).

In the frequency domain, this material substitution introduces a complex change to the reflected signal:

$$M(k\_z, k\_x) = B(k\_z, k\_x) \\cdot \\mathbf{\\Delta R}(k\_z, k\_x)$$

Where $\\mathbf{\\Delta R}$ is the complex ratio of the fluid-filled reflection coefficient to the air-filled reflection coefficient. For a subwavelength thin layer, this ratio can be linearized to an amplitude attenuation factor ($\\alpha$) and a localized phase lag ($\\Delta \\theta$):

$$\\mathbf{\\Delta R} \\approx \\alpha \\cdot e^{-j \\Delta \\theta}$$

### Interpreting the 2D Cross-Spectrum

If you pass the baseline and monitor B-scans of a fluid-filling event into your existing `estimate_shift_2d` function, the cross-spectrum becomes:

$$XS(k\_z, k\_x) = |B(k\_z, k\_x)|^2 \\cdot \\alpha \\cdot e^{j \\Delta \\theta}$$

When you extract the phase $\\Phi = \\text{angle}(XS)$, look at how it compares to the physical movement equation:

-   **Physical Movement:** $\\Phi(k\_z, k\_x) = k\_z \\Delta z + k\_x \\Delta x$ (A tilted plane sloping through the origin).

-   **Material Change:** $\\Phi(k\_z, k\_x) = \\Delta \\theta$ (A perfectly **flat, horizontal plane** shifted uniformly by an intercept value $\\Delta \\theta$, with zero slope).

## How to Read the Output for Material Changes

Your existing function solves for three variables: `c[0]` ($\\Delta z$), `c[1]` ($\\Delta x$), and `c[2]` (the constant offset $c$).

When analyzing a localized fracture during a fluid injection experiment:

1.  **The Slopes ($\\Delta z, \\Delta x$) will be near zero.** Because the fracture walls aren't mechanically moving, the phase plane will not tilt.

2.  **The Intercept ($c$) will change dynamically.** The constant offset term `c[2]`, which was previously just absorbing minor system noise, now becomes your primary physical metric.

The phase intercept $c$ tracks the fluid volume inside the fracture footprint:

$$c = \\Delta \\theta \\approx \\frac{2\\pi d}{\\lambda} \\left( \\sqrt{\\varepsilon\_{\\text{fluid}}} - \\sqrt{\\varepsilon\_{\\text{baseline}}} \\right)$$

By tracking the change in $c$ over time-lapse intervals, you can map out exactly when and where the fluid front arrives, even if the total thickness of the fracture $d$ is only a few millimeters wide.

## The Adjusted Workflow for Fluid Tracking

To adapt your pipeline for fluid monitoring, follow this workflow structure:

```container
[ Time-Lapse Field Data ] ──► [ Standard Pre-Processing ] ──► [ Migration ]
                                                                   │
                                                                   ▼
[ Map Out Intercept Changes ] ◄── [ Extract Constant Term 'c' ] ◄── [ Spatial Window Loop (ROI) ]
```

### 1\. The Spatial Window Loop (Moving ROI)

Unlike a point scatterer that moves to a single new location, a fluid front spreads along a fracture plane.

-   **Action:** Instead of isolating a single region of interest, set up a loop that slides a small window (e.g., $16 \\times 16$ pixels) along the known geometry of the fracture in your migrated images.

### 2\. Extract the Phase Intercept Profile

Run `estimate_shift_2d` inside each window along the fracture path. Store the extracted constant term $c$ for each spatial position.

Python

```container
# Pseudo-code for a spatial profile loop along a horizontal fracture at depth row `z_frac`
window_size = 16
intercept_profile = []

for x_center in range(window_size, Nx - window_size, 4):
    roi_base = mig_base[z_frac - 8 : z_frac + 8,  x_center - 8 : x_center + 8]
    roi_mon  = mig_mon[z_frac - 8  : z_frac + 8,  x_center - 8 : x_center + 8]

    # Run your updated function
    dz, dx, XS, kz_ax, kx_ax = estimate_shift_2d(roi_base, roi_mon, dz_mig, dx_mig, kz_c)

    # Re-extract the constant term (intercept) directly from the least-squares output
    # In our script, c was the third element of the result vector
    # To grab it directly, modify estimate_shift_2d to return the full 'c' vector
    phase_intercept = estimate_intercept_only(roi_base, roi_mon, dz_mig, dx_mig, kz_c)
    intercept_profile.append((x_center, phase_intercept))
```

### 3\. Threshold to Map the Fluid Front

Plot the `phase_intercept` profile across the length of the fracture. Where the fracture is still dry/air-filled, the intercept will hover near $0^\\circ$. At the locations where the fluid has penetrated, the phase intercept will experience a sharp step-change (e.g., climbing to a stable $35^\\circ$ phase lag). This allows you to track the real-time velocity of a fluid front flowing through a subwavelength boundary in the field.

To understand the **cross phase spectrum**, we need to break it down into what it is conceptually, how it is derived mathematically, and the physical units it uses.

At its core, the cross phase spectrum isolates exactly how much two signals are "out of sync" with each other, evaluated independently at every single frequency.

## 1\. What Exactly is the Cross Phase Spectrum?

When you take the Fourier transform of a single GPR B-scan, the result is a spectrum containing complex numbers. Every complex number has two parts:

-   **Magnitude:** How strong that frequency is.

-   **Phase:** The starting angle (or timing) of that frequency's sine wave.

The absolute phase of a single image is usually a chaotic, wrapped mess because it depends entirely on the shape of your radar pulse and where the target happens to be in the ground.

The **Cross-Spectrum ($XS$)** is a mathematical trick to strip away that absolute phase and look only at the *difference* between a baseline image ($B$) and a monitor image ($M$). You compute it by multiplying the baseline spectrum by the complex conjugate of the monitor spectrum:

$$XS = B \\cdot M^\*$$

If we represent the spectra in polar form ($B = |B|e^{j\\phi\_b}$ and $M = |M|e^{j\\phi\_m}$):

$$XS = \\Big(|B|e^{j\\phi\_b}\\Big) \\cdot \\Big(|M|e^{-j\\phi\_m}\\Big)$$

$$XS = |B||M| e^{j(\\phi\_b - \\phi\_m)}$$

The **Cross Phase Spectrum ($\\Phi$)** is simply the angle of this result:

$$\\Phi = \\text{angle}(XS) = \\phi\_b - \\phi\_m$$

**Conceptually:** It is a map that tells you, "For this specific frequency, the monitor wave is shifted by exactly $X$ amount relative to the baseline wave." By using the cross-spectrum, the messy absolute phase cancels itself out, leaving you with nothing but the pure phase discrepancy caused by the movement of your scatterer or fluid.

## 2\. What Units Does it Have?

There are three different units at play when looking at a 2D cross phase spectrum plot, and understanding how they interact is the key to why the Fourier shift algorithm works.

### The Value (The Z-Axis)

The value of the cross phase spectrum ($\\Phi$) itself is an angle. Therefore, its unit is strictly **radians (rad)** or **degrees ($^\\circ$)**.

-   A value of $0$ radians means that specific frequency didn't shift at all.

-   A value of $\\pi$ radians ($180^\\circ$) means the wave shifted by exactly half a cycle (destructive interference).

### The Coordinates (The X and Y Axes)

The domain that the cross phase spectrum lives in is the frequency domain. In 2D spatial GPR images, these axes are your wavenumbers ($k\_z$ and $k\_x$).

-   The unit for wavenumber is **radians per meter (rad/m)**.

### The Magic Unit Cancellation (The Slope)

This is where the math beautifully aligns with the physics. In the previous Weighted Least Squares (WLS) algorithm, we were calculating the *slope* of the cross phase spectrum to find the displacement $\\Delta z$.

Slope is the "rise over run".

-   The "rise" is a change in phase ($\\Delta \\Phi$). Unit: **radians**

-   The "run" is a change in wavenumber ($\\Delta k\_z$). Unit: **radians / meter**

If we divide them to find the slope:

$$\\text{Slope} = \\frac{\\Delta \\Phi}{\\Delta k\_z} = \\frac{\[\\text{radians}\]}{\[\\text{radians} / \\text{meter}\]}$$

The radians cancel out, the meter moves to the numerator, and the unit of the slope is exactly **meters (m)**.

This is why the cross phase spectrum is so powerful. It directly converts an angular rotation in the frequency domain back into a pure physical distance, allowing you to measure millimeter-scale shifts from phase data.