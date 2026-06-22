import json, sys
sys.stdout.reconfigure(encoding='utf-8')
NB_PATH = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\TimeLapse_Processing.ipynb'
with open(NB_PATH, encoding='utf-8') as f:
    nb = json.load(f)
for c in nb['cells']:
    if c.get('id') == 'c6b39323':
        src = c['source']
        if isinstance(src, list):
            src = ''.join(src)
        print(src)
        break
