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
        # Print larger windows around each block
        for keyword, window in [
            ('Carrier-phase detrending helpers', 1000),
            ('Col 2: CLSSA phase spectrum', 600),
        ]:
            idx = src.find(keyword)
            if idx >= 0:
                # find start of line
                start = src.rfind('\n', 0, idx) + 1
                print(f'=== {keyword!r} ===')
                print(repr(src[start:start+window]))
                print()
        break
