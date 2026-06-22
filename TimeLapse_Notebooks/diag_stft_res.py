import json, sys
sys.stdout.reconfigure(encoding='utf-8')
NB_PATH = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\TimeLapse_Processing.ipynb'
with open(NB_PATH, encoding='utf-8') as f:
    nb = json.load(f)

# Find f_c and dz_mig / v_ice definitions by scanning all code cells
for c in nb['cells']:
    if c['cell_type'] != 'code':
        continue
    src = c['source']
    if isinstance(src, list):
        src = ''.join(src)
    for kw in ['f_c =', 'f_c=', 'dz_mig =', 'dz_mig=', 'v_ice =', 'v_ice=']:
        idx = src.find(kw)
        if idx >= 0:
            end = src.find('\n', idx)
            print(f'[{c.get("id")}] {src[idx:end]}')
