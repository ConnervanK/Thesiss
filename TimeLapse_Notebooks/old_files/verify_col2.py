import json, sys
sys.stdout.reconfigure(encoding='utf-8')
NB_PATH = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\TimeLapse_Processing.ipynb'
CELL_ID = '2790c368'
with open(NB_PATH, encoding='utf-8') as f:
    nb = json.load(f)
for c in nb['cells']:
    if c.get('id') == CELL_ID:
        src = c['source'] if isinstance(c['source'], str) else ''.join(c['source'])
        checks = {
            'plt.subplots(4, 5':       True,   # 4 rows
            'figsize=(28, 20)':        True,   # taller figure
            '_ROW_SPECS':              True,   # row spec list
            'col4_type == \'2d\'':     True,   # 2-D branch
            'dphi_2d_0':               True,   # 2-D heatmap for x_sc0
            'dphi_2d_m':               True,   # 2-D heatmap for x_sc_mon
            'dphi_1d_0':               True,   # 1-D slice for x_sc0
            'dphi_1d_m':               True,   # 1-D slice for x_sc_mon
            '_dphi_2d':                True,   # 2-D helper
            '_dphi_1d':                True,   # 1-D helper
            'BEFORE  (baseline survey)': True, # row labels
            'AFTER   (monitor survey)':  True,
            '_detrend_gather':         True,   # gather detrend kept
            '_detrend_phase':          False,  # removed
            '_carrier_2d':             False,  # removed
        }
        all_ok = True
        for term, expected in checks.items():
            found = term in src
            ok = found == expected
            all_ok = all_ok and ok
            print(f'{"OK" if ok else "FAIL"} | present={found} expected={expected} | {term!r}')
        print()
        print('ALL OK' if all_ok else 'SOME CHECKS FAILED')
        break
