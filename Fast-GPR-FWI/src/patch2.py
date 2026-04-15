import glob
import re

for file in glob.glob('C:/Users/Administrator/Thesis/Fast-GPR-FWI/src/lib/*.cu'):
    with open(file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Remove existing ones
    content = content.replace('__declspec(dllexport)', '')
    content = content.replace('__declspec(dllexport) ', '')
    
    # Add dllexport to all void functions not preceded by __global__ or __device__
    content = re.sub(r'(?<!__global__\s)(?<!__device__\s)(?<!__global__)(?<!__device__)void\s+([A-Za-z0-9_]+)\s*\(', r'__declspec(dllexport) void \1(', content)
    
    with open(file, 'w', encoding='utf-8') as f:
        f.write(content)
