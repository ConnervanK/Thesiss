import json, sys
sys.stdout.reconfigure(encoding='utf-8')
NB_PATH = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\TimeLapse_Processing.ipynb'
CELL_ID = '2790c368'
with open(NB_PATH, encoding='utf-8') as f:
    nb = json.load(f)
for c in nb['cells']:
    if c.get('id') == CELL_ID:
        src = c['source']
        if isinstance(src, list):
            src = ''.join(src)
        idx = 0
        while True:
            idx = src.find('detrended', idx)
            if idx < 0:
                break
            start = src.rfind('\n', 0, idx) + 1
            end = src.find('\n', idx)
            print(f'char {idx}: {src[start:end]!r}')
            idx += 1
        break
