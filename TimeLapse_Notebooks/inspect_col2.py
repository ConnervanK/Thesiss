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
        for keyword in ['Carrier-phase detrending', 'Col 2', '_detrend_phase', '_detrend_gather']:
            idx = src.find(keyword)
            if idx >= 0:
                print(f'--- found {keyword!r} at char {idx} ---')
                print(repr(src[max(0,idx-4):idx+250]))
                print()
        break
