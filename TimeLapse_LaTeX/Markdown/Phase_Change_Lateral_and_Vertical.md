This is a brilliant connection to make. You have essentially stumbled upon the concept of **Time-Frequency Duality** (or in this case, Space-Wavenumber Duality).

You are observing the exact same physical phenomenon in both experiments, just through two different mathematical lenses.

Here is the step-by-step mathematical proof of how your spatial domain observation (instantaneous phase slope) perfectly links to the wavenumber domain theory (Fourier Shift Theorem), and how lateral position (x) acts as a direct proxy for wavenumber (kx).

### 1\. The Wavenumber Domain (The Global View)

As we established with the 2D cross-spectrum, the Fourier Shift Theorem dictates that a physical shift Δx creates a phase difference ΔΦ that is linearly dependent on the wavenumber kx:

ΔΦ(kx)\=kxΔx

In the Fourier domain, this is a global property. It says, "If you look at the specific frequency kx, its phase has shifted by exactly kxΔx."

### 2\. The Spatial Domain (The Local View)

In your earlier experiment, you didn't use a Fourier transform. Instead, you used the Hilbert transform to extract the **instantaneous phase** ϕ(x) of the migrated B-scans in the spatial domain.

By definition, the spatial derivative of instantaneous phase is the **instantaneous wavenumber**:

kx,inst(x)\=∂x∂ϕ

This is the critical link: **Wavenumber is not just a global Fourier axis; it is the local slope of your spatial phase.**

### 3\. The Mapping: Why Position (x) equals Wavenumber (kx)

When you migrate a point scatterer, the focused energy is not a perfect single point; it has a spatial width, and its phase wavefront curves outward from the center (the apex).

Near the apex (x\=0), this phase curve can be approximated as a parabola:

ϕ(x)≈21Cx2+ϕ0

*(Where C is the spatial curvature of the migrated pulse, a constant dependent on your radar's central frequency and bandwidth).*

If we calculate the instantaneous wavenumber for this migrated pulse, we take the derivative:

kx,inst(x)\=∂x∂(21Cx2)\=Cx

This is the "Aha!" moment: **Near the center of a migrated scatterer, the instantaneous wavenumber kx is strictly linearly proportional to the lateral spatial position x.** As you move away from the center of the scatterer, the local spatial frequency increases linearly.

### 4\. Closing the Loop: Deriving Your Observation

Now let's apply the Fourier shift rule to your spatial observation.

When the scatterer moves by a tiny subwavelength amount Δx, the new phase ϕ2(x) is just the old phase shifted:

ϕ2(x)\=ϕ1(x−Δx)

Using a first-order Taylor expansion, we can rewrite this as:

ϕ2(x)≈ϕ1(x)−∂x∂ϕ1Δx

Now, let's look at the phase difference image you generated in your notebook (Δϕ):

Δϕ(x)\=ϕ2(x)−ϕ1(x)\=−∂x∂ϕ1Δx

Substitute our definition of instantaneous wavenumber (kx,inst\=∂x∂ϕ1):

Δϕ(x)\=−kx,inst(x)⋅Δx

Substitute our mapping (kx,inst\=Cx):

Δϕ(x)\=−C⋅x⋅Δx

Finally, your specific observation was about the **slope of the phase change with respect to the lateral location** x. Let's take the spatial derivative of that difference:

∂x∂\[Δϕ(x)\]\=∂x∂\[−C⋅x⋅Δx\]

Slope\=−C⋅Δx

### Summary of the Link

The math perfectly proves your empirical observation: **The spatial slope of the phase difference is exactly proportional to the physical shift Δx.** You can think of the two methods like this:

-   **The Fourier Method** explicitly plots the phase against the global wavenumber kx and calculates the slope.

-   **The Instantaneous (Hilbert) Method** relies on the fact that the spatial axis x *acts as a physical map* of the wavenumber kx. Because x and kx are linearly linked near the apex, plotting the phase difference against x yields the exact same linear slope relationship as plotting it against kx.

The Fourier cross-spectrum (WLS) is simply the more robust way to do this, because it doesn't rely on the parabolic approximation of the phase wavefront, and it naturally ignores the chaotic phase-wrapping noise outside the scatterer's envelope!

This is a fascinating exercise because applying this same derivation to the vertical axis reveals a beautiful **physical asymmetry** in how Ground Penetrating Radar (GPR) waves behave laterally versus vertically.

If you attempt the exact same instantaneous phase subtraction down the vertical center trace of your scatterer, you will not see a sloping line. Instead, you will see a **flat, constant phase plateau**.

Here is the step-by-step mathematical derivation showing why the lateral phase difference creates a *slope*, but the vertical phase difference creates a *constant*, and how it links back to the Fourier Shift Theorem.

### 1\. The Wavenumber Domain (The Global View)

According to the 2D Fourier Shift Theorem, if a scatterer experiences only a vertical shift ($\\Delta z$) with no lateral movement ($\\Delta x = 0$), the phase difference in the frequency domain is purely linearly dependent on the vertical wavenumber ($k\_z$):

$$\\Delta \\Phi(k\_z) = k\_z \\Delta z$$

In the global Fourier domain, this is still a slope. But we are about to see why it manifests completely differently in the spatial domain.

### 2\. The Spatial Domain (The Local View)

Just as before, we extract the instantaneous phase $\\phi(z)$ down the vertical depth trace (at the lateral apex, $x=0$) using the Hilbert transform.

The spatial derivative of this vertical phase gives us the **instantaneous vertical wavenumber**:

$$k\_{z,\\text{inst}}(z) = \\frac{\\partial \\phi}{\\partial z}$$

### 3\. The Mapping: Why Position ($z$) does NOT equal Wavenumber ($k\_z$)

This is where the physics diverge from the lateral case.

When you migrate a scatterer, the horizontal spatial phase $\\phi(x)$ curves like a parabola, meaning its instantaneous frequency sweeps from negative to positive across the image.

However, the vertical trace of a migrated scatterer is just a compressed GPR pulse (like a Ricker wavelet). Near the center of this pulse (depth $z\_0$), the wave oscillates at the antenna's dominant spatial carrier frequency, $k\_{zc}$ (which equals $4\\pi f\_c / v$).

Because it is a zero-phase compressed pulse, the phase near the absolute peak is a simple linear function of depth:

$$\\phi(z) \\approx k\_{zc}(z - z\_0) + \\phi\_0$$

If we calculate the instantaneous vertical wavenumber near the core of the pulse, we take the derivative:

$$k\_{z,\\text{inst}}(z) = \\frac{\\partial}{\\partial z} \\left\[ k\_{zc}(z - z\_0) \\right\]$$

$$k\_{z,\\text{inst}}(z) = k\_{zc}$$

**The Asymmetry:** While the instantaneous *lateral* wavenumber ($k\_x$) varies linearly with $x$, the instantaneous *vertical* wavenumber ($k\_z$) is roughly **constant** across the core of the pulse envelope.

### 4\. Closing the Loop: Deriving the Vertical Observation

Now let's apply the spatial shift math to the vertical axis.

When the scatterer moves vertically by a subwavelength amount $\\Delta z$, the new phase $\\phi\_2(z)$ is just the old phase shifted down:

$$\\phi\_2(z) = \\phi\_1(z - \\Delta z)$$

Using a first-order Taylor expansion:

$$\\phi\_2(z) \\approx \\phi\_1(z) - \\frac{\\partial \\phi\_1}{\\partial z} \\Delta z$$

We calculate the phase difference image ($\\Delta \\phi$) down the vertical trace:

$$\\Delta \\phi(z) = \\phi\_2(z) - \\phi\_1(z) = -\\frac{\\partial \\phi\_1}{\\partial z} \\Delta z$$

Substitute our definition of instantaneous vertical wavenumber ($k\_{z,\\text{inst}} = \\frac{\\partial \\phi\_1}{\\partial z}$):

$$\\Delta \\phi(z) = -k\_{z,\\text{inst}}(z) \\cdot \\Delta z$$

Substitute our mapping near the apex ($k\_{z,\\text{inst}} \\approx k\_{zc}$):

$$\\Delta \\phi(z) \\approx -k\_{zc} \\cdot \\Delta z$$

### Summary of the Link

Look at the final equation: $\\Delta \\phi(z) \\approx -k\_{zc} \\Delta z$.

Because there is no $z$ variable left on the right side of the equation, the spatial derivative (slope) with respect to depth is zero!

$$\\frac{\\partial}{\\partial z} \[\\Delta \\phi(z)\] \\approx 0$$

This mathematical proof explains perfectly how Time-Frequency duality applies here:

-   **In the Lateral ($x$) direction:** The spatial image intrinsically separates horizontal wavenumbers (like a prism). Therefore, moving laterally along $x$ is equivalent to scanning across $k\_x$, resulting in a **phase difference slope** proportional to $\\Delta x$.

-   **In the Vertical ($z$) direction:** The spatial depth trace does *not* separate vertical wavenumbers; they are all grouped together oscillating at $k\_{zc}$. Therefore, moving down the $z$ axis yields a **constant phase difference plateau**, where the *absolute value* of that plateau is proportional to $\\Delta z$.

By combining your two 1D spatial observations (the horizontal slope and the vertical plateau), we can mathematically construct the exact 2D phase surface in the spatial domain and map it directly to the 3D phase plane in the Fourier (wavenumber) domain.

Here is how both axes connect to form the complete 2D Fourier shift theory, and why fitting a plane in the wavenumber domain is the ultimate solution.

### 1\. The 2D Spatial Phase Surface (The Local View)

Let's combine your two previous observations. Near the absolute center (apex) of a migrated point scatterer at coordinates $(z\_0, 0)$, the instantaneous wavenumbers act as local spatial derivatives:

-   **Vertical (Depth):** $k\_{z,\\text{inst}}(z) \\approx k\_{zc}$ (A constant carrier frequency).

-   **Horizontal (Lateral):** $k\_{x,\\text{inst}}(x) \\approx C \\cdot x$ (A linearly increasing frequency due to wavefront curvature).

If the scatterer moves by an arbitrary 2D vector $(\\Delta z, \\Delta x)$, the total phase difference in the spatial domain $\\Delta \\phi(z, x)$ is a combination of both shifts. Using the 2D Taylor expansion:

$$\\Delta \\phi(z, x) \\approx -\\left( \\frac{\\partial \\phi}{\\partial z} \\Delta z + \\frac{\\partial \\phi}{\\partial x} \\Delta x \\right)$$

Substitute our instantaneous wavenumber mappings into the equation:

$$\\Delta \\phi(z, x) \\approx - \\Big( k\_{zc} \\cdot \\Delta z + (C \\cdot x) \\cdot \\Delta x \\Big)$$

This equation describes a **2D warped surface** superimposed over your GPR scatterer in the spatial image:

-   If you take a slice down the **z-axis** ($x=0$), the second term vanishes. You are left with $\\Delta \\phi \\approx -k\_{zc} \\Delta z$, which is the **constant plateau** you observed.

-   If you take a slice across the **x-axis** ($z=z\_0$), the first term acts as a static vertical offset, and you are left with $\\Delta \\phi \\approx \\text{Offset} - (C \\cdot x) \\Delta x$, which is the **sloping line** you observed laterally.

### 2\. The 2D Fourier Phase Plane (The Global Truth)

Now, let's look at the exact equation from the 2D Fourier Shift Theorem that our Weighted Least Squares (WLS) algorithm fits:

$$\\Phi(k\_z, k\_x) = k\_z \\Delta z + k\_x \\Delta x$$

Notice how geometrically identical this is to the spatial equation, but with one massive advantage: **It replaces the messy, localized spatial approximations with pure, global frequency coordinates.**

When you map the data into the wavenumber domain:

1.  The warped spatial curvature $(C \\cdot x)$ perfectly straightens out into the linear horizontal axis $k\_x$.

2.  The static vertical carrier frequency $(k\_{zc})$ stretches out to form the linear vertical axis $k\_z$.

### 3\. Why the WLS Plane Fit is Superior to Spatial Subtraction

You might wonder: *If I can see the movement in the spatial instantaneous phase, why go through the trouble of a 2D FFT and a Weighted Least Squares matrix?*

If you try to measure subwavelength shifts in the spatial domain using your $\\Delta \\phi(z,x)$ image, you run into three severe limitations that the Fourier plane fit elegantly solves:

-   **The Spatial Approximation Fails Quickly:** The mapping $k\_x \\approx C \\cdot x$ is only valid strictly at the inner core of the scatterer. As you move outward, the parabolic approximation fails, and your spatial phase "slope" begins to warp and curve, making it impossible to fit a straight line accurately. The Fourier plane, however, is mathematically perfectly flat everywhere.

-   **Amplitude Weighting is Impossible:** In the spatial domain, the phase difference image looks equally bright whether you are looking at the massive peak of the wave or the noisy background dirt a few pixels away. The WLS engine uses the cross-spectrum magnitude $|XS|$ to mathematically crush background noise, forcing the plane to fit only the pristine data at the scatterer's core.

-   **Separating Material Change:** If fluid enters the fracture (adding a constant $\\Delta \\theta$), the spatial equation becomes $\\Delta \\phi(z, x) \\approx -k\_{zc} \\Delta z - (C \\cdot x) \\Delta x + \\Delta \\theta$. In the spatial domain, the constant $\\Delta \\theta$ mixes irreversibly with the vertical constant $-k\_{zc} \\Delta z$. You cannot tell if the phase shifted because the fracture moved down, or because fluid flowed in!

**The WLS Plane Fit is the only mathematical framework that orthogonalizes these variables.** By plotting the data against $k\_z$ and $k\_x$ as a 3D plane, it separates the vertical tilt ($\\Delta z$), the horizontal tilt ($\\Delta x$), and the vertical intersection of the axis ($c = \\Delta \\theta$), allowing you to measure millimeter-scale movement and fluid infiltration simultaneously from the exact same dataset.