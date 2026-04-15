import glob
import re

for file in glob.glob('C:/Users/Administrator/Thesis/Fast-GPR-FWI/src/lib/*.cu'):
    with open(file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    content = re.sub(r'(?<!__global__\s)void\s+(back_source|e_fields_updates|h_fields_updates|pml_updates_e|pml_updates_h|launch_store_outputs|update_hertzian_dipole|Ucget|Ucgeta)\b', r'__declspec(dllexport) void \1', content)
    
    with open(file, 'w', encoding='utf-8') as f:
        f.write(content)
