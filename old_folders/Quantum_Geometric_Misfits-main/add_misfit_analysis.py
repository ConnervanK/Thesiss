import re

fname = '/home/marcohuy/thesis_master/latex_doc/methods/misfit_functions.tex'
with open(fname, 'r') as f:
    text = f.read()

figure3_text = r"""
\subsection{Analysis of Objective Function Landscapes}
A direct comparison of the objective function landscapes highlights the advantages of the Geometric misfit over the standard Euclidean L2 formulation. The Euclidean surface often exhibits severe local minima and sharp ridges due to the source-signature dependency and cycle-skipping. In contrast, the Variable Projection (Geometric) misfit projects out the source dependency, yielding a structurally smoother landscape that broadens the basin of attraction and enhances the robustness of gradient-based optimizers.

\begin{figure}[htbp]
    \centering
    \includegraphics[width=0.9\textwidth]{figures/misfit_analysis.png}
    \caption{Topological analysis comparing the traditional L2 Euclidean misfit landscape with the Variable Projection Geometric misfit landscape.}
    \label{fig:misfit_analysis}
\end{figure}
"""

text += '\n' + figure3_text

with open(fname, 'w') as f:
    f.write(text)

