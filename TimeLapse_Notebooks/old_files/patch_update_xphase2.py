import json, sys

NB_PATH    = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\TimeLapse_Processing.ipynb'
MD_PATH    = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\new_xphase_md.md'
CODE_PATH  = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\new_xphase_code.py'
MD_ID      = 'dc220149'
CODE_ID    = '1d78e0f4'

sys.stdout.reconfigure(encoding='utf-8')

with open(MD_PATH, encoding='utf-8') as f:
    md_src = f.read()
with open(CODE_PATH, encoding='utf-8') as f:
    code_src = f.read()

with open(NB_PATH, encoding='utf-8') as f:
    nb = json.load(f)

patched = {'md': False, 'code': False}
for c in nb['cells']:
    if c.get('id') == MD_ID:
        c['source'] = md_src
        patched['md'] = True
    elif c.get('id') == CODE_ID:
        c['source'] = code_src
        c['outputs'] = []
        c['execution_count'] = None
        patched['code'] = True

if not all(patched.values()):
    print(f'ERROR: missing cells: {patched}', file=sys.stderr)
    sys.exit(1)

with open(NB_PATH, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f'OK — markdown {MD_ID} and code {CODE_ID} updated to CLSSA + zoom/clip')
