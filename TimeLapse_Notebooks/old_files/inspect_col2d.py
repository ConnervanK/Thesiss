import json, sys

NB_PATH = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\TimeLapse_Processing.ipynb'
CELL_ID = '2790c368'

sys.stdout.reconfigure(encoding='utf-8')

with open(NB_PATH, encoding='utf-8') as f:
    nb = json.load(f)

for c in nb['cells']:
    if c.get('id') == CELL_ID:
        src = c['source']
        if isinstance(src, list):
            src = ''.join(src)

        # Get the full detrend helpers block
        keyword = 'Carrier-phase detrending helpers'
        idx = src.find(keyword)
        start = src.rfind('\n', 0, idx) + 1
        print(repr(src[start:start+1200]))
        break
