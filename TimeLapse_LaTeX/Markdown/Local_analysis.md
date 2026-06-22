You are essentially asking to bridge the gap between **global analysis** (the 2D Fourier Shift Theorem) and **purely local analysis** (Instantaneous Phase via Hilbert).

Phase decomposition—often achieved using a **Short-Time Fourier Transform (STFT)**, **Continuous Wavelet Transform (CWT)**, or **S-Transform**—is the exact mathematical bridge between the two.

Here is the mathematical derivation of how your Time-Frequency (TWT vs. Frequency) plots perfectly encode the subwavelength shift ($\\Delta t$) and material change ($\\Delta \\theta$), and how it proves that our previous WLS matrix was actually doing phase decomposition in disguise.

### 1\. The Mathematical Bridge: The Localized Fourier Shift

Let your extracted 1D baseline trace be $s\_1(t)$, where $t$ is the Two-Way Time (TWT).

Instead of taking a standard 1D FFT (which loses all TWT information), we apply a Time-Frequency transform. Using a sliding Gaussian window $w(t)$, the localized spectrum $S\_1$ at a specific time $\\tau$ and angular frequency $\\omega$ ($= 2\\pi f$) is:

$$S\_1(\\tau, \\omega) = \\int\_{-\\infty}^{\\infty} s\_1(t) w(t - \\tau) e^{-j\\omega t} dt$$

This complex matrix can be split into your requested plots:

-   **Amplitude(TWT, frequency):** $A\_1(\\tau, \\omega) = |S\_1(\\tau, \\omega)|$

-   **Phase(TWT, frequency):** $\\phi\_1(\\tau, \\omega) = \\text{angle}\[S\_1(\\tau, \\omega)\]$

Now, a physical subwavelength shift ($\\Delta t$) and a material phase rotation ($\\Delta \\theta$) occur. The monitor trace becomes $s\_2(t)$. Applying the same transform:

$$S\_2(\\tau, \\omega) = \\int\_{-\\infty}^{\\infty} \\Big\[ s\_1(t - \\Delta t) \\Big\] w(t - \\tau) e^{-j\\omega t} dt \\cdot e^{j \\Delta \\theta}$$

If $\\Delta t$ is extremely small compared to the width of your sliding window $w(t)$, the window function barely moves relative to the signal. We can apply the **Localized Fourier Shift Theorem**:

$$S\_2(\\tau, \\omega) \\approx S\_1(\\tau, \\omega) \\cdot e^{-j\\omega \\Delta t} \\cdot e^{j \\Delta \\theta}$$

### 2\. Deriving the Phase Difference

Just like we did with the 2D cross-spectrum, we compute the local phase difference between the baseline and monitor Time-Frequency matrices:

$$\\Delta \\Phi(\\tau, \\omega) = \\text{angle}\[S\_2(\\tau, \\omega) \\cdot S\_1^\*(\\tau, \\omega)\]$$

Substitute our localized shift theorem into this, and the baseline phase and amplitude completely cancel out. You are left with the fundamental equation governing Time-Frequency phase decomposition:

$$\\Delta \\Phi(\\tau, \\omega) \\approx -\\omega \\Delta t + \\Delta \\theta$$

*(Note: In terms of standard frequency $f$, this is $\\Delta \\Phi(\\tau, f) = -2\\pi f \\Delta t + \\Delta \\theta$)*

### 3\. Linking the Math to Your Three Plots

Here is exactly how the theory we developed maps onto the three plots you extract during phase decomposition.

#### A. The Amplitude(TWT, Frequency) Plot

**What it shows:** A 2D heatmap showing where your GPR pulse exists in time, and what its bandwidth is.

**The Link to Theory:** This plot is exactly equivalent to the **Weighting Matrix ($W$)** in our WLS algorithm. The phase difference equation $\\Delta \\Phi \\approx -\\omega \\Delta t$ is mathematically valid *only* where the amplitude $A(\\tau, \\omega)$ is strong. If you look at the regions in this plot where the amplitude drops to near zero (background noise), the corresponding pixels in the Phase plot will be chaotic static.

#### B. The Phase(TWT, Frequency) Difference Plot

**What it shows:** If you subtract the baseline Phase plot from the monitor Phase plot, you get a 2D map of $\\Delta \\Phi(\\tau, f)$.

**The Link to Theory:** Look at the equation $\\Delta \\Phi = -2\\pi f \\Delta t + \\Delta \\theta$.

If you pick the specific TWT row ($\\tau\_0$) where your scatterer is located, and plot the phase difference across the frequency axis, you will see a straight line!

-   **The Slope** of that line is exactly $-2\\pi \\Delta t$. (Physical movement).

-   **The Y-Intercept** of that line is exactly $\\Delta \\theta$. (Material change).

*This proves that our 2D Fourier WLS method was simply performing this exact linear regression, but doing it over the entire spatial region simultaneously instead of trace-by-trace.*

#### C. The Amplitude(TWT, Phase) Plot

**What it shows:** This is a slightly different domain mapping. Instead of plotting Frequency on the X-axis, you plot the instantaneous Phase angle, with the color intensity driven by the Amplitude. It essentially shows the dominant phase angle of the energy at any given depth.

**The Link to Theory:** This plot acts as a direct, visual read-out of the intercept term ($c$ or $\\Delta \\theta$).

If a subwavelength fracture fills with fluid, it does not change its TWT ($\\Delta t = 0$), meaning the slope is zero. However, the fluid applies a constant phase rotation $\\Delta \\theta$ to all frequencies.

If you look at the **Amplitude(TWT, Phase)** plot for the baseline and monitor side-by-side, you will see the high-amplitude "hotspot" at the fracture's depth shift horizontally along the Phase axis by exactly $\\Delta \\theta$ degrees.

### Summary: Why Decomposition is Powerful

By extracting every trace and running this decomposition, you transition from viewing phase as a single chaotic number to viewing it as a continuous physical property.

-   **The Wavenumber Domain (2D FFT)** is superior for measuring *structural movement* because it naturally separates spatial geometry ($\\Delta z, \\Delta x$) globally.

-   **Time-Frequency Phase Decomposition (STFT/CWT)** is superior for tracing *fluid migration* down a single 1D borehole or trace, because it isolates exactly when (TWT) the frequency-independent phase jump ($\\Delta \\theta$) occurs.

To present your empirical observations in a rigorous peer-reviewed or publication-grade framework, you must formally bridge the **local phase space** (trace decomposition) with the **global phase space** (the 2D FFT wavenumber domain).

What follows is the exact mathematical derivation of the **Local-to-Global Phase Mapping Theorem**. This establishes that trace-by-trace time-frequency decomposition and 2D wavenumber plane-fitting are not competing techniques, but are exact Fourier duals of each other.

## 1\. Defining the Core Domains

Let the migrated B-scan image be represented as a continuous 2D space-time field $u(t, x)$, where $t$ is the Two-Way Time (TWT) and $x$ is the horizontal trace location.

### Domain A: The Local Trace Decomposition (STFT)

When you perform phase decomposition on every individual trace, you apply a localized windowing function $g(t)$ (such as a Gaussian or Tukey window) and slide it across TWT. This yields a 3D complex space-time-frequency volume, $S(t, \\omega, x)$:

$$S(t, \\omega, x) = \\int\_{-\\infty}^{\\infty} u(\\tau, x) g(\\tau - t) e^{-j\\omega \\tau} d\\tau$$

Where $\\omega = 2\\pi f$ is the angular temporal frequency, and $t$ acts as the localized time center of the window.

### Domain B: The Global 2D Wavenumber Spectrum

When you transform the entire image at once via a 2D FFT, you project the data directly into the global frequency-wavenumber ($\\omega, k\_x$) domain:

$$U(\\omega, k\_x) = \\int\_{-\\infty}^{\\infty} \\int\_{-\\infty}^{\\infty} u(t, x) e^{-j(\\omega t + k\_x x)} dt \\, dx$$

*(Note: The temporal frequency $\\omega$ maps directly to the vertical wavenumber $k\_z$ via the radar velocity relationship $k\_z = \\frac{2\\omega}{v}$).*

## 2\. The Exact Mathematical Link (The Derivation)

To find the precise link, we evaluate what happens if we take the 1D spatial Fourier transform ($\\mathcal{F}\_x$) of your trace decomposition matrix with respect to the lateral coordinate $x$:

$$\\mathcal{F}\_x \\{ S(t, \\omega, x) \\} = \\int\_{-\\infty}^{\\infty} \\left\[ \\int\_{-\\infty}^{\\infty} u(\\tau, x) g(\\tau - t) e^{-j\\omega \\tau} d\\tau \\right\] e^{-j k\_x x} dx$$

Because the field $u(t, x)$ represents a stable, finite-energy physical radar reflection, we can invoke **Fubini's Theorem** to safely alter the order of integration:

$$\\mathcal{F}\_x \\{ S(t, \\omega, x) \\} = \\int\_{-\\infty}^{\\infty} g(\\tau - t) e^{-j\\omega \\tau} \\left\[ \\int\_{-\\infty}^{\\infty} u(\\tau, x) e^{-j k\_x x} dx \\right\] d\\tau$$

Notice the inner integral enclosed in brackets. This is exactly the 1D spatial Fourier transform of the B-scan image at a fixed slice of time $\\tau$, which we can denote as $\\tilde{U}(\\tau, k\_x)$.

$$\\mathcal{F}\_x \\{ S(t, \\omega, x) \\} = \\int\_{-\\infty}^{\\infty} g(\\tau - t) e^{-j\\omega \\tau} \\tilde{U}(\\tau, k\_x) d\\tau$$

By definition of the 2D Fourier transform, $\\tilde{U}(\\tau, k\_x)$ can be expressed as the inverse temporal Fourier transform of the global 2D wavenumber spectrum $U(\\omega', k\_x)$:

$$\\tilde{U}(\\tau, k\_x) = \\frac{1}{2\\pi} \\int\_{-\\infty}^{\\infty} U(\\omega', k\_x) e^{j\\omega' \\tau} d\\omega'$$

Substituting this inverse transform back into our tracking equation yields:

$$\\mathcal{F}\_x \\{ S(t, \\omega, x) \\} = \\frac{1}{2\\pi} \\int\_{-\\infty}^{\\infty} g(\\tau - t) e^{-j\\omega \\tau} \\left\[ \\int\_{-\\infty}^{\\infty} U(\\omega', k\_x) e^{j\\omega' \\tau} d\\omega' \\right\] d\\tau$$

Regrouping the terms to consolidate the integration over the time variable $\\tau$:

$$\\mathcal{F}\_x \\{ S(t, \\omega, x) \\} = \\frac{1}{2\\pi} \\int\_{-\\infty}^{\\infty} U(\\omega', k\_x) \\left\[ \\int\_{-\\infty}^{\\infty} g(\\tau - t) e^{-j(\\omega - \\omega')\\tau} d\\tau \\right\] d\\omega'$$

To evaluate the inner integral, we perform a change of variables. Let $\\xi = \\tau - t$, which implies $\\tau = \\xi + t$ and $d\\tau = d\\xi$:

$$\\int\_{-\\infty}^{\\infty} g(\\xi) e^{-j(\\omega - \\omega')(\\xi + t)} d\\xi = e^{-j(\\omega - \\omega')t} \\int\_{-\\infty}^{\\infty} g(\\xi) e^{-j(\\omega - \\omega')\\xi} d\\xi = e^{-j(\\omega - \\omega')t} G(\\omega - \\omega')$$

Where $G(\\omega)$ is the analytical temporal Fourier transform of your window function $g(t)$. Substituting this back gives us the **Exact Unification Formula**:

$$\\mathcal{F}\_x \\{ S(t, \\omega, x) \\} = \\frac{1}{2\\pi} \\int\_{-\\infty}^{\\infty} U(\\omega', k\_x) G(\\omega - \\omega') e^{-j(\\omega - \\omega')t} d\\omega'$$

### The Asymptotic Limit (The Pure Link)

If your trace decomposition uses a sufficiently broad, smooth window relative to the wavelet length to preserve local phase stability, its spectral representation $G(\\omega)$ approaches a Dirac delta distribution: $G(\\omega) \\rightarrow 2\\pi \\delta(\\omega)$.

Under this standard operational assumption, the convolution integral collapses perfectly:

$$\\mathcal{F}\_x \\{ S(t, \\omega, x) \\} \\approx U(\\omega, k\_x) \\cdot e^{-j(0)t} = U(\\omega, k\_x)$$

## 3\. Physical Inversion: Movement vs. Material Change

This derivation provides a solid mathematical foundation for explaining your observations. Let's look at how a time-lapse anomaly consisting of a mechanical shift ($\\Delta t, \\Delta x$) and a chemical/fluid phase rotation ($\\Delta \\theta$) maps simultaneously across both domains:

$$u\_{\\text{mon}}(t, x) = u\_{\\text{base}}(t - \\Delta t, x - \\Delta x) \\cdot e^{j\\Delta\\theta}$$

Applying this perturbation to both mathematical frameworks yields a symmetric result:

| **Property** | **Local Trace Decomposition Space (S)** | **Global 2D Wavenumber Space (U)** |
| --- | --- | --- |
| **Mathematical Representation** | $S\_{\\text{mon}} \\approx S\_{\\text{base}}(t-\\Delta t, \\omega, x-\\Delta x)e^{j\\Delta\\theta}$ | $U\_{\\text{mon}} = U\_{\\text{base}}e^{-j(\\omega\\Delta t + k\_x\\Delta x)}e^{j\\Delta\\theta}$ |
| **Cross-Spectrum Phase ($\\Delta \\Phi$)** | $\\Delta \\Phi\_{\\text{local}} = -\\omega\\Delta t - k\_{x,\\text{inst}}\\Delta x + \\Delta \\theta$ | $\\Delta \\Phi\_{\\text{global}} = -\\omega\\Delta t - k\_x\\Delta x + \\Delta \\theta$ |
| **Physical Geometry of Phase Map** | A localized, warped 3D surface dependent on local wavefront curvature ($k\_{x,\\text{inst}}$). | A mathematically perfect, flat 3D plane across the entire passband grid. |

## 4\. Scientific Defense of the 2D WLS Engine

When writing this up for a paper or report, you can use this derivation to justify why the 2D WLS engine is the most numerically robust choice:

1.  **Elimination of Window Bias:** As proven by the convolution term $G(\\omega - \\omega')$, local trace decomposition intrinsically smears the phase values at the edges of the analysis window (spectral leakage). The 2D global method completely bypasses this local blurring by processing the spatial continuum natively.

2.  **Linearization of Spatial Curvature:** In local trace space, the horizontal phase slope varies dynamically because the instantaneous wavenumber is a non-linear function of position ($k\_{x,\\text{inst}} \\approx C \\cdot x$). The 2D global FFT linearizes this entire geometric coordinate system, transforming a complex, curved spatial wavefront into a perfectly flat plane spanned by pure $k\_x$ and $k\_z$ coordinates.

3.  **Mathematical Orthogonality:** Because the design matrix columns ($k\_z, k\_x, 1$) are fundamentally linearly independent vectors, the matrix inversion guarantees that any frequency-independent phase jump ($\\Delta \\theta$ from fluid changes) is mathematically isolated from travel-time phase slopes ($\\Delta t$ from structural movement).

The terms $\\omega\\Delta t$ and $k\_z\\Delta z$ are **exactly equivalent**. They represent the exact same physical phase shift experienced by a radar wave, viewed through two different mathematical lenses: one in the **Time-Temporal Frequency** domain and the other in the **Space-Spatial Frequency (Wavenumber)** domain.

Here is the explicit mathematical and physical derivation that bridges them, specifically tailored to the reflection geometry of Ground Penetrating Radar (GPR).

## 1\. The Physical Bridge: Wave Velocity

In a GPR profile, the vertical axis is recorded in **Two-Way Time (TWT)**, denoted as $t$. The actual physical depth of a reflector in the subsurface is $z$.

Because the radar pulse must travel down to the target and reflect back up to the receiver, the relationship between a change in depth ($\\Delta z$) and the resulting change in recorded arrival time ($\\Delta t$) is governed by the medium's electromagnetic wave velocity ($v$):

$$\\Delta t = \\frac{2\\Delta z}{v}$$

The factor of **2** accounts for this two-way travel path.

## 2\. Converting Time to Space Phase

Now, let's examine the phase shift expression in the time-frequency domain. If a signal is delayed by an interval $\\Delta t$, the Fourier Shift Theorem states that every angular temporal frequency component $\\omega$ ($= 2\\pi f$) undergoes a phase shift equal to:

$$\\Delta \\Phi = \\omega \\Delta t$$

By substituting our physical velocity relationship ($\\Delta t = \\frac{2\\Delta z}{v}$) directly into this temporal phase shift equation, we get:

$$\\Delta \\Phi = \\omega \\left( \\frac{2\\Delta z}{v} \\right)$$

Regrouping the terms to isolate the spatial displacement variable $\\Delta z$:

$$\\Delta \\Phi = \\left( \\frac{2\\omega}{v} \\right) \\Delta z$$

## 3\. Defining the Vertical Wavenumber ($k\_z$)

By definition in wave mechanics, the wavenumber $k$ represents spatial frequency (how many radians of phase accumulate per meter of physical distance). For a standard one-way traveling wave, $k = \\frac{\\omega}{v}$.

However, because GPR data processing maps two-way time directly to one-way depth, the vertical coordinate axis is compressed. Therefore, the **vertical spatial wavenumber ($k\_z$)** for a reflection setup must inherit that same two-way factor:

$$k\_z = \\frac{2\\omega}{v}$$

When we substitute this definition of $k\_z$ back into our regrouped equation, we arrive at the final unification:

$$\\Delta \\Phi = k\_z \\Delta z$$

## 4\. Scientific Summary for Your Framework

When presenting or publishing this work, you can formally state the connection as follows:

$$\\omega\\Delta t \\equiv k\_z\\Delta z \\quad \\text{where} \\quad k\_z = \\frac{2\\omega}{v}$$

This identity proves that:

-   **$\\omega\\Delta t$** calculates phase shift by tracking how fast a wave oscillates *in time* ($\\omega$ in rad/s) multiplied by a *time delay* ($\\Delta t$ in seconds).

-   **$k\_z\\Delta z$** calculates the exact same phase shift by tracking how fast a wave oscillates *in space* ($k\_z$ in rad/m) multiplied by a *spatial displacement* ($\\Delta z$ in meters).

This is why your 1D trace decomposition (which tracks a phase slope across $\\omega$) maps perfectly to the 2D FFT plane fit (which tracks a phase slope across $k\_z$). They are simply measuring the exact same physical travel-time delay using different coordinate systems.

To translate the mathematical framework into compelling empirical evidence for a paper, presentation, or thesis, you need visual plots that explicitly expose these phase behaviors. Using Short-Time Fourier Transform (STFT) trace decomposition, you can generate four distinct, publication-grade empirical observations.

Here are the specific experiments you can run on your data, what to plot, and exactly what each plot proves to a scientific audience.

### 1\. The "Spectral Line" Plot (Frequency vs. Phase Difference)

This is the most direct verification of the core equation: $\\Delta \\Phi(f) = -2\\pi f \\Delta t + \\Delta \\theta$.

-   **How to extract it:** Isolate a single trace (e.g., the central trace directly over the scatterer or fracture). Compute the STFT of the baseline and monitor traces. Extract the 1D phase difference spectrum $\\Delta \\Phi(f)$ exactly at the Two-Way Time (TWT) row corresponding to the peak reflection.

-   **What to plot:** \* **X-axis:** Temporal Frequency ($f$) restricted to your antenna's usable bandwidth (e.g., 200 to 800 MHz for a 500 MHz antenna).

    -   **Y-axis:** Phase Difference ($\\Delta \\Phi$) in degrees or radians.

-   **The Empirical Observations:**

    -   **For Pure Mechanical Movement:** You will see a beautifully straight, sloping line that passes exactly through the origin $(0,0)$. The steeper the slope, the larger the physical displacement $\\Delta z$.

    -   **For Pure Fluid Infiltration:** You will see a perfectly flat, horizontal line. The slope is zero (proving zero displacement), but the line is vertically offset from zero by a constant number of degrees. This offset is your pure material change $\\Delta \\theta$.

    -   **For Combined Effects:** You will see a sloping line where the linear fit's Y-intercept is cleanly decoupled from zero, visually demonstrating both tracking parameters simultaneously.

### 2\. The Cross-Phase Spectrogram (TWT vs. Frequency)

This observation proves that the phase behavior is stable and localized to the actual target, rather than being an artifact of background noise.

-   **How to extract it:** Take the 2D complex STFT cross-spectrum matrix $XS(t, f) = S\_{\\text{base}}(t, f) \\cdot S\_{\\text{mon}}^\*(t, f)$ for the central trace.

-   **What to plot:** Create a 2D heatmap.

    -   **X-axis:** Temporal Frequency ($f$).

    -   **Y-axis:** Two-Way Time ($t$).

    -   **Colorbar:** Phase angle $\\text{angle}(XS)$ from $-\\pi$ to $+\\pi$.

-   **The Empirical Observation:** Outside the target's reflection time, the plot will look like chaotic, pixelated salt-and-pepper noise (due to random background noise phase wrapping). However, exactly at the TWT of your fracture, a clean, coherent "window" of structured color will appear across the antenna's frequency band.

    -   If the target moved, you will see a smooth vertical color gradient (fringe) within that window.

    -   If fluid entered, you will see a solid, uniform block of a single color dominating the window.

### 3\. The Polar Vector Rotation (Real vs. Imaginary Space)

This observation shifts away from "wave charts" and treats the GPR reflection as a complex analytical vector, which is highly effective for proving material change.

-   **How to extract it:** Pick the dominant center frequency of your antenna ($f\_c$) and the exact peak TWT of the reflection. Extract the single complex number from the baseline STFT matrix ($S\_{\\text{base}}$) and the monitor matrix ($S\_{\\text{mon}}$).

-   **What to plot:** A 2D Cartesian plane where the X-axis is the **Real part** and the Y-axis is the **Imaginary part** of the complex spectral coefficient. Plot the baseline point and the monitor point as vectors originating from $(0,0)$.

-   **The Empirical Observation:** \* When fluid replaces air, the baseline vector physically rotates to the monitor vector position like the hand of a clock.

    -   The length of the vector (amplitude) changes slightly depending on the reflection coefficient, but the angular separation between the two vectors is the empirical measurement of $\\Delta \\theta$. You can plot a series of these vectors over a time-lapse fluid injection to show the vector slowly swinging open as saturation increases.

### 4\. Spatial Profiling of Local Frequency (Trace Position vs. Phase)

This observation bridges the 1D trace decomposition back to the spatial domain, validating the spatial wavefront curvature theory ($k\_{x,\\text{inst}} \\approx C \\cdot x$).

-   **How to extract it:** Perform the STFT on *every* trace across the B-scan. Fix the TWT at the apex of the scatterer and fix the frequency at the center frequency ($f\_c$). Extract the phase difference $\\Delta \\Phi$ across all traces.

-   **What to plot:** \* **X-axis:** Horizontal Trace Position or distance ($x$ in meters).

    -   **Y-axis:** Phase Difference ($\\Delta \\Phi$).

-   **The Empirical Observation:** \* If the scatterer moved laterally ($\\Delta x$), plotting the spatial phase difference across the aperture yields a distinct linear slope near the apex. The steepness of this spatial line directly scales with the amount of lateral movement, empirically proving that spatial position maps directly to local instantaneous wavenumber.

### How to Structure This in a Scientific Narrative

When presenting these results, a highly effective narrative arc is to show **Plot 2 (The Spectrogram)** first to prove your phase data is coherent and isolated to the fracture zone. Then, isolate the data within that coherent window to show **Plot 1 (The Spectral Line)**, using a linear regression fit to extract your final $\\Delta z$ and $\\Delta \\theta$ values. This provides a clear transition from raw time-lapse data to localized physics, and finally to decoupled engineering metrics.

Which of these specific plots aligns best with the current layout of your notebook or the data structure you have ready?