To fully appreciate why the 2D Fourier shift method is so sensitive to both mechanical movement and fluid infiltration, we have to look closely at the linear algebra governing it.

Here is the full mathematical derivation of the Weighted Linear Least Squares (WLS) engine, from the foundational equations down to the matrix inversion, and exactly how it orthogonalizes (separates) movement from material change.

## 1\. Defining the System of Equations

From the Fourier Shift Theorem and complex cross-spectrum, we established that the unwrapped phase difference $\\Phi$ at any given 2D frequency bin $(k\_z, k\_x)$ is modeled by:

$$\\Phi(k\_z, k\_x) = k\_z \\Delta z + k\_x \\Delta x + c$$

In a discrete GPR image, you do not have a continuous function; you have $N$ individual frequency bins inside your passband mask. This gives us a system of $N$ linear equations:

$$\\Phi\_1 = k\_{z,1} \\Delta z + k\_{x,1} \\Delta x + c$$

$$\\Phi\_2 = k\_{z,2} \\Delta z + k\_{x,2} \\Delta x + c$$

$$\\vdots$$

$$\\Phi\_N = k\_{z,N} \\Delta z + k\_{x,N} \\Delta x + c$$

We can express this efficiently in matrix notation as $A \\mathbf{u} = \\mathbf{\\Phi}$:

$$\\begin{bmatrix} k\_{z,1} & k\_{x,1} & 1 \\\\ k\_{z,2} & k\_{x,2} & 1 \\\\ \\vdots & \\vdots & \\vdots \\\\ k\_{z,N} & k\_{x,N} & 1 \\end{bmatrix} \\begin{bmatrix} \\Delta z \\\\ \\Delta x \\\\ c \\end{bmatrix} = \\begin{bmatrix} \\Phi\_1 \\\\ \\Phi\_2 \\\\ \\vdots \\\\ \\Phi\_N \\end{bmatrix}$$

Where:

-   **$A$** is the $N \\times 3$ "design matrix" containing our frequency coordinates.

-   **$\\mathbf{u}$** is the $3 \\times 1$ vector of the unknown parameters we want to find.

-   **$\\mathbf{\\Phi}$** is the $N \\times 1$ observation vector of measured phase differences.

Because real data contains noise, the equation is never perfectly equal. There is always an error, or residual vector $\\mathbf{e}$:

$$\\mathbf{\\Phi} = A \\mathbf{u} + \\mathbf{e}$$

## 2\. Why Ordinary Least Squares (OLS) Fails

The goal of Ordinary Least Squares (OLS) is to find the vector $\\mathbf{u}$ that minimizes the sum of the squared residuals. Mathematically, it minimizes the cost function $S(\\mathbf{u})$:

$$S(\\mathbf{u}) = \\mathbf{e}^T \\mathbf{e} = (\\mathbf{\\Phi} - A \\mathbf{u})^T (\\mathbf{\\Phi} - A \\mathbf{u})$$

Setting the derivative with respect to $\\mathbf{u}$ to zero yields the classic OLS solution:

$$\\mathbf{u} = (A^T A)^{-1} A^T \\mathbf{\\Phi}$$

**The fatal flaw in GPR data:** OLS assumes every single equation (every frequency bin) is equally reliable. In a GPR spectrum, a frequency bin at the absolute peak power of the antenna is incredibly stable. A bin at the very edge of the bandwidth is dominated by thermal background noise. OLS gives the noisy edge bin the exact same voting power as the peak signal bin, completely corrupting the plane fit.

## 3\. The Weighted Least Squares (WLS) Derivation

To fix this, we introduce a diagonal weight matrix $W$ of size $N \\times N$. The diagonal entries are the magnitudes of the cross-spectrum, $W\_{ii} = |XS\_i|$, so that high-energy bins get massive weights and low-energy bins approach zero.

$$W = \\begin{bmatrix} |XS\_1| & 0 & \\dots & 0 \\\\ 0 & |XS\_2| & \\dots & 0 \\\\ \\vdots & \\vdots & \\ddots & \\vdots \\\\ 0 & 0 & \\dots & |XS\_N| \\end{bmatrix}$$

We redefine our cost function to minimize the **weighted** sum of squared residuals:

$$S\_W(\\mathbf{u}) = \\mathbf{e}^T W^2 \\mathbf{e} = (\\mathbf{\\Phi} - A \\mathbf{u})^T W^2 (\\mathbf{\\Phi} - A \\mathbf{u})$$

To find the minimum, we expand the equation and take the gradient with respect to $\\mathbf{u}$:

$$S\_W(\\mathbf{u}) = \\mathbf{\\Phi}^T W^2 \\mathbf{\\Phi} - 2 \\mathbf{u}^T A^T W^2 \\mathbf{\\Phi} + \\mathbf{u}^T A^T W^2 A \\mathbf{u}$$

Taking the derivative $\\frac{\\partial}{\\partial \\mathbf{u}}$ and setting it to zero:

$$-2 A^T W^2 \\mathbf{\\Phi} + 2 A^T W^2 A \\mathbf{u} = 0$$

Rearranging to solve for the unknown parameters $\\mathbf{u}$:

$$A^T W^2 A \\mathbf{u} = A^T W^2 \\mathbf{\\Phi}$$

$$\\mathbf{u} = (A^T W^2 A)^{-1} A^T W^2 \\mathbf{\\Phi}$$

This is the exact, closed-form solution for the Weighted Linear Least Squares estimator.

> **Computational Note (What Python is actually doing):** > Calculating $(A^T W^2 A)^{-1}$ directly can be numerically unstable if the matrix is ill-conditioned. Instead, your Python code (`np.linalg.lstsq`) avoids squaring the condition number by multiplying both $A$ and $\\mathbf{\\Phi}$ by $W$ directly:
>
> $$(W A) \\mathbf{u} = (W \\mathbf{\\Phi})$$
>
> It then solves this pre-weighted system using Singular Value Decomposition (SVD). This achieves the exact same WLS result but with vastly superior numerical precision.

## 4\. Why This Infers Movement vs. Material Change

The beauty of formulating the problem as a 3D matrix $\\mathbf{u} = \[\\Delta z, \\Delta x, c\]^T$ is that the math forces the algorithm to completely decouple spatial movement from material/phase changes.

Here is exactly how the math interprets the physics:

### A. Inferring Movement ($\\Delta z, \\Delta x$)

By the properties of the Fourier transform, a physical translation in space strictly manifests as a phase shift that is **linearly dependent on frequency**.

-   If a scatterer moves down by $1\\text{ mm}$, the phase shift at $100\\text{ MHz}$ might be $5^\\circ$, but the phase shift at $500\\text{ MHz}$ will be exactly $25^\\circ$.

-   The WLS matrix explicitly hunts for this frequency-dependent slope ($k\_z$ and $k\_x$). It projects all the phase variance that scales with frequency strictly into the $\\Delta z$ and $\\Delta x$ variables.

### B. Inferring Material Change ($c$ or $\\Delta \\theta$)

When fluid replaces air in a subwavelength fracture, it changes the bulk dielectric permittivity. This alters the complex reflection coefficient.

-   Unlike physical movement, a change in the reflection coefficient applies a **constant phase rotation** across the entire wavelet, regardless of frequency. The phase shift at $100\\text{ MHz}$ will be $15^\\circ$, and the phase shift at $500\\text{ MHz}$ will *also* be $15^\\circ$.

-   In our equation $\\Phi = k\_z \\Delta z + k\_x \\Delta x + c$, a constant phase shift has a slope of zero. Therefore, the WLS engine mathematically forces all non-frequency-dependent phase shifts directly into the intercept term $c$.

Because the columns of our design matrix $A$ separate the linear terms ($k\_z, k\_x$) from the constant term ($1$), the WLS engine guarantees that a fluid injection event will exclusively drive the value of $c$ up or down, while leaving $\\Delta z$ and $\\Delta x$ cleanly anchored at zero.