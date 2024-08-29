
import numpy as np
import struct

def load_lut_from_file(path: str, dim: int):
    with open(path, 'rb') as lut:
        
        
        size = [int.from_bytes(lut.read(4), 'little') for x in range(dim)]
        minmax = [(struct.unpack('f', lut.read(4))[0], struct.unpack('f', lut.read(4))[0]) for x in range(dim)]
        
        buffer = lut.read()


        print(size)
        print(minmax)

        flatarray = np.frombuffer(buffer, dtype=np.float32)
        flatarray = np.reshape(flatarray, size)

        return flatarray
