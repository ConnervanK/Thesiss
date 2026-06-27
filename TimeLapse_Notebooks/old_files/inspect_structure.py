import json, sys
sys.stdout.reconfigure(encoding='utf-8')
NB_PATH = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\TimeLapse_Processing.ipynb'
with open(NB_PATH, encoding='utf-8') as f:
    nb = json.load(f)

for idx, c in enumerate(nb['cells']):
    src = c['source']
    if isinstance(src, list):
        src = ''.join(src)
    first_line = src.split('\n', 1)[0][:90]
    print(f'[{idx:3}] id={c.get("id")!r:14} type={c["cell_type"]:8} | {first_line}')
