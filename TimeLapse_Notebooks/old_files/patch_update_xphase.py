import json, sys

NB_PATH  = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\TimeLapse_Processing.ipynb'
SRC_PATH = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\new_xphase_code.py'
CELL_ID  = '1d78e0f4'

sys.stdout.reconfigure(encoding='utf-8')

with open(SRC_PATH, encoding='utf-8') as f:
    new_src = f.read()

with open(NB_PATH, encoding='utf-8') as f:
    nb = json.load(f)

patched = False
for c in nb['cells']:
    if c.get('id') == CELL_ID:
        c['source'] = new_src
        c['outputs'] = []
        c['execution_count'] = None
        patched = True
        break

if not patched:
    print('ERROR: cell not found', file=sys.stderr)
    sys.exit(1)

with open(NB_PATH, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f'OK — cell {CELL_ID} updated with longer STFT window ({len(new_src)} chars)')
