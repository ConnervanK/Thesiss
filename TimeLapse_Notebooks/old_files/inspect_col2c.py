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

        # Find the Col 2 plotting block inside the loop
        keyword = 'Col 2: CLSSA phase spectrum'
        # Find the second occurrence (inside the loop)
        idx1 = src.find(keyword)
        idx2 = src.find(keyword, idx1+1)
        print(f'First occurrence at {idx1}, second at {idx2}')
        if idx2 >= 0:
            start = src.rfind('\n', 0, idx2) + 1
            print(repr(src[start:start+700]))
        break
