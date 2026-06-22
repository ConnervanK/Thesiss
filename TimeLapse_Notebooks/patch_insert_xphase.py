import json, sys, secrets

NB_PATH    = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\TimeLapse_Processing.ipynb'
MD_PATH    = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\new_xphase_md.md'
CODE_PATH  = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\new_xphase_code.py'
AFTER_ID   = '2790c368'   # Section 7 code cell — insert new cells right after this
OLD_SEC8_MD_ID   = '27f4b168'
OLD_SEC8_CODE_ID = 'c6b39323'

sys.stdout.reconfigure(encoding='utf-8')

with open(MD_PATH, encoding='utf-8') as f:
    md_src = f.read()
with open(CODE_PATH, encoding='utf-8') as f:
    code_src = f.read()

with open(NB_PATH, encoding='utf-8') as f:
    nb = json.load(f)

new_md_cell = {
    'cell_type': 'markdown',
    'id': secrets.token_hex(4),
    'metadata': {},
    'source': md_src,
}
new_code_cell = {
    'cell_type': 'code',
    'id': secrets.token_hex(4),
    'metadata': {},
    'execution_count': None,
    'outputs': [],
    'source': code_src,
}

cells = nb['cells']
insert_at = None
for idx, c in enumerate(cells):
    if c.get('id') == AFTER_ID:
        insert_at = idx + 1
        break
if insert_at is None:
    print('ERROR: AFTER_ID not found', file=sys.stderr)
    sys.exit(1)

cells[insert_at:insert_at] = [new_md_cell, new_code_cell]

# Renumber the old "Section 8" -> "Section 9"
renamed = 0
for c in cells:
    if c.get('id') == OLD_SEC8_MD_ID:
        src = c['source']
        if isinstance(src, list):
            src = ''.join(src)
        src = src.replace('## Section 8 —', '## Section 9 —', 1)
        c['source'] = src
        renamed += 1
    elif c.get('id') == OLD_SEC8_CODE_ID:
        src = c['source']
        if isinstance(src, list):
            src = ''.join(src)
        src = src.replace('Section 8: Localized Fourier Shift', 'Section 9: Localized Fourier Shift', 1)
        c['source'] = src
        c['outputs'] = []
        c['execution_count'] = None
        renamed += 1

with open(NB_PATH, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f'OK — inserted new markdown id={new_md_cell["id"]}, code id={new_code_cell["id"]} '
      f'at index {insert_at}; renamed {renamed} old Section-8 cells to Section 9')
