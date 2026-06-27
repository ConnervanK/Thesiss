# Force LuaLaTeX: figure filenames under TimeLapse_Figures contain Unicode
# (lambda, Delta, Phi, fraction glyphs), which pdfTeX cannot resolve reliably.
$pdf_mode = 4;
$lualatex = 'lualatex -interaction=nonstopmode -synctex=1 -file-line-error %O %S';
$pdflatex = $lualatex;
