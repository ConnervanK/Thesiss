To understand exactly why the Weighted Linear Least Squares (WLS) engine can perfectly isolate fluid flow from physical movement, we need to bridge the **electromagnetic physics** of a thin fracture with the **linear algebra** of the phase matrix.

The core reason it works is **Orthogonality**: physical movement changes wave *travel time*, which creates a frequency-dependent phase slope. Fluid injection changes wave *impedance*, which creates a frequency-independent phase jump. The WLS matrix is mathematically built to separate slopes from constants.

Here is the deep dive into the physics and the math of why the intercept term ($c$) perfectly captures fluid movement.

### 1\. The Physics: Why Fluid Causes a "Flat" Phase Shift

When a GPR wave hits a fracture, it actually reflects twice: once off the top boundary, and once off the bottom boundary.

If the fracture is thick (wider than the wavelength, $d > \\lambda$), you see two distinct pulses in your B-scan. But if the fracture is **subwavelength** ($d \\ll \\lambda$), the two reflections physically overlap and interfere with each other, creating a single, composite wavelet.

The total complex reflection coefficient for this thin layer can be approximated by:

$$R\_{\\text{total}} \\approx R\_{\\text{top}} + R\_{\\text{bottom}} e^{-j 2 k d}$$

Because the layer is extremely thin ($k d \\ll 1$), we can use a first-order Taylor expansion on the exponential. When you substitute the dielectric properties of the rock matrix ($\\varepsilon\_m$) and the material inside the fracture ($\\varepsilon\_f$), the reflection coefficient simplifies to:

$$R\_{\\text{total}} \\approx -j \\omega \\frac{d \\sqrt{\\varepsilon\_m}}{2 c\_0} \\left( \\frac{\\varepsilon\_m - \\varepsilon\_f}{\\varepsilon\_f} \\right)$$

Notice the $-j \\omega$ term. In the time domain, multiplying by $-j\\omega$ is equivalent to taking the derivative. This means a thin fracture naturally reflects the *derivative* of the incident GPR pulse.

**Now, what happens when fluid replaces air?**

The thickness $d$ stays exactly the same, but $\\varepsilon\_f$ changes drastically (e.g., air $\\varepsilon\_r = 1$, water $\\varepsilon\_r = 81$).

The baseline air-filled reflection is $R\_{\\text{base}}$, and the monitor fluid-filled reflection is $R\_{\\text{mon}}$. When the Fourier Shift algorithm computes the cross-spectrum, it calculates the ratio of these two states:

$$\\mathbf{\\Delta R} = \\frac{R\_{\\text{mon}}}{R\_{\\text{base}}} = \\frac{\\varepsilon\_{\\text{air}} (\\varepsilon\_m - \\varepsilon\_{\\text{water}})}{\\varepsilon\_{\\text{water}} (\\varepsilon\_m - \\varepsilon\_{\\text{air}})}$$

Because the geometric $d$ and the derivative $-j\\omega$ terms completely cancel out in this ratio, **the resulting change is a purely real or purely imaginary constant**.

-   It does not depend on the spatial frequency $k\_z$ or $k\_x$.

-   In the polar/phase domain, this translates to a uniform, constant phase rotation $\\Delta \\theta$ across the entire bandwidth of the GPR pulse.

### 2\. The Math: How the WLS Matrix Isolates $\\Delta \\theta$

Let's look at the observation vector $\\mathbf{\\Phi}$ inside the WLS engine. For a fluid injection event with zero physical movement ($\\Delta z = 0$, $\\Delta x = 0$), the true phase difference at every frequency bin $i$ is exactly the same constant value:

$$\\Phi\_i = \\Delta \\theta$$

If we feed this into our overdetermined linear system $A \\mathbf{u} = \\mathbf{\\Phi}$:

$$\\begin{bmatrix} k\_{z,1} & k\_{x,1} & 1 \\\\ k\_{z,2} & k\_{x,2} & 1 \\\\ \\vdots & \\vdots & \\vdots \\\\ k\_{z,N} & k\_{x,N} & 1 \\end{bmatrix} \\begin{bmatrix} \\Delta z \\\\ \\Delta x \\\\ c \\end{bmatrix} = \\begin{bmatrix} \\Delta \\theta \\\\ \\Delta \\theta \\\\ \\vdots \\\\ \\Delta \\theta \\end{bmatrix}$$

The least-squares algorithm minimizes the error. It looks at the columns of the design matrix $A$:

1.  **Column 1 ($k\_z$):** Values range from negative to positive frequencies. If the algorithm assigns any non-zero value to $\\Delta z$, the predicted phase will slope upwards or downwards. But the target $\\mathbf{\\Phi}$ is flat. Therefore, to minimize error, the algorithm *must* set $\\Delta z = 0$.

2.  **Column 2 ($k\_x$):** Same logic. A non-zero $\\Delta x$ tilts the plane sideways. To match a flat target, $\\Delta x$ *must* equal $0$.

3.  **Column 3 (The Intercept):** This column is a vector of pure $1$s. If the algorithm sets $c = \\Delta \\theta$, the predicted output perfectly matches the flat target vector.

**The Orthogonality Principle:**

Because the vectors representing the frequencies ($k\_z, k\_x$) are mathematically linearly independent (orthogonal) to the constant vector of $1$s, the matrix inversion can completely decouple them.

If a fracture physically shifts downwards *while* fluid flows through it, the phase plane will both tilt *and* shift vertically. The WLS engine simply measures the tilt to calculate the movement ($\\Delta z$) and measures where the plane pierces the Z-axis origin to calculate the fluid saturation ($c$).