import os
import sys
import ctypes

import devito.data.allocators as allocators

orig_alloc = allocators.PosixAllocator._alloc_C_libcall

def win32_alloc_C_libcall(self, size, ctype):
    msvcrt = ctypes.cdll.msvcrt
    msvcrt._aligned_malloc.restype = ctypes.c_void_p
    msvcrt._aligned_malloc.argtypes = [ctypes.c_size_t, ctypes.c_size_t]
    
    alignment = 64
    pointer = msvcrt._aligned_malloc(size, alignment)
    if not pointer:
        raise MemoryError("Couldn't allocate memory using `_aligned_malloc`")
        
    print(f"Allocated pointer {pointer} of size {size}")
    
    c_pointer = ctypes.cast(pointer, ctypes.POINTER(ctype))
    def memfree():
        try:
            msvcrt._aligned_free.restype = None
            msvcrt._aligned_free.argtypes = [ctypes.c_void_p]
            msvcrt._aligned_free(pointer)
        except Exception:
            pass
        
    return c_pointer, {"memfree": memfree}
    
allocators.PosixAllocator._alloc_C_libcall = win32_alloc_C_libcall

from devito import configuration
configuration['log-level'] = 'DEBUG'

from devito import Grid, Function
grid = Grid(shape=(10, 10))
m = Function(name='loc', grid=grid)
try:
    print(m.data)
except Exception as e:
    import traceback
    traceback.print_exc()



