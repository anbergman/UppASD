import ovflib
import ctypes
import numpy as np

### Load Library
_ovf = ovflib.LoadOVFLibrary()

### Return codes
OK      = -1
ERROR   = -2
INVALID = -3

### File formats
FILEFORMAT_BIN  = -53
FILEFORMAT_BIN4 = -54
FILEFORMAT_BIN8 = -55
FILEFORMAT_TEXT = -56
FILEFORMAT_CSV  = -57

class ovf_segment(ctypes.Structure):
    ### Some properties
    _fields_ = [
        ("title",            ctypes.c_char_p),
        ("comment",          ctypes.c_char_p),
        ("valuedim",         ctypes.c_int),
        ("valueunits",       ctypes.c_char_p),
        ("valuelabels",      ctypes.c_char_p),
        ("meshtype",         ctypes.c_char_p),
        ("meshunits",        ctypes.c_char_p),
        ("pointcount",       ctypes.c_int),
        ("n_cells",          ctypes.c_int*3),
        ("N",                ctypes.c_int),
        ("bounds_min",       ctypes.c_float*3),
        ("bounds_max",       ctypes.c_float*3),
        ("lattice_constant", ctypes.c_float),
        ("bravais_vectors",  ctypes.c_float*3*3)
    ]

    def __init__(self, title="", comment="", valuedim=1, valueunits="", valuelabels="", meshtype="", meshunits="",
                bounds_min=[0.0, 0.0, 0.0], bounds_max=[0.0, 0.0, 0.0], lattice_constant=0.0,
                bravais_vectors=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                pointcount=0, n_cells=[1,1,1]):

        self.title       = title
        self.comment     = comment
        self.valuedim    = valuedim
        self.valueunits  = valueunits
        self.valuelabels = valuelabels
        self.meshtype    = meshtype
        self.meshunits   = meshunits
        self.pointcount  = pointcount
        for i in range(3):
            self.n_cells[i] = n_cells[i]
        self.N = n_cells[0]*n_cells[1]*n_cells[2]
        for i in range(3):
            self.bounds_min[i] = bounds_min[i]
            self.bounds_max[i] = bounds_max[i]
            for j in range(3):
                self.bravais_vectors[i][j]  = bravais_vectors[i][j]
        self.lattice_constant = lattice_constant

### --------------------------------------------------------------

### Read a segment header
_ovf_read_segment_header = _ovf.ovf_read_segment_header
_ovf_read_segment_header.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ovf_segment)]
_ovf_read_segment_header.restype  = ctypes.c_int

### Read a segment with float precision
_ovf_read_segment_data_4 = _ovf.ovf_read_segment_data_4
_ovf_read_segment_data_4.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ovf_segment), ctypes.POINTER(ctypes.c_float)]
_ovf_read_segment_data_4.restype  = ctypes.c_int

### Read a segment with double precision
_ovf_read_segment_data_8 = _ovf.ovf_read_segment_data_8
_ovf_read_segment_data_8.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ovf_segment), ctypes.POINTER(ctypes.c_double)]
_ovf_read_segment_data_8.restype  = ctypes.c_int

### Write a segment with float precision (overwrite file)
_ovf_write_segment_4 = _ovf.ovf_write_segment_4
_ovf_write_segment_4.argtypes = [ctypes.c_void_p, ctypes.POINTER(ovf_segment), ctypes.POINTER(ctypes.c_float), ctypes.c_int]
_ovf_write_segment_4.restype  = ctypes.c_int

### Write a segment with double precision (overwrite file)
_ovf_write_segment_8 = _ovf.ovf_write_segment_8
_ovf_write_segment_8.argtypes = [ctypes.c_void_p, ctypes.POINTER(ovf_segment), ctypes.POINTER(ctypes.c_double), ctypes.c_int]
_ovf_write_segment_8.restype  = ctypes.c_int

### Append a segment with float precision
_ovf_append_segment_4 = _ovf.ovf_append_segment_4
_ovf_append_segment_4.argtypes = [ctypes.c_void_p, ctypes.POINTER(ovf_segment), ctypes.POINTER(ctypes.c_float), ctypes.c_int]
_ovf_append_segment_4.restype  = ctypes.c_int

### Append a segment with double precision
_ovf_append_segment_8 = _ovf.ovf_append_segment_8
_ovf_append_segment_8.argtypes = [ctypes.c_void_p, ctypes.POINTER(ovf_segment), ctypes.POINTER(ctypes.c_double), ctypes.c_int]
_ovf_append_segment_8.restype  = ctypes.c_int

### Fetch the latest message
_ovf_latest_message = _ovf.ovf_latest_message
_ovf_latest_message.argtypes = [ctypes.c_void_p]
_ovf_latest_message.restype  = ctypes.c_char_p

### --------------------------------------------------------------


class _ovf_file(ctypes.Structure):
    ### Some properties
    _fields_ = [
        ("filename",     ctypes.c_char_p),
        ("found",        ctypes.c_bool),
        ("is_ovf",       ctypes.c_bool),
        ("n_segments",   ctypes.c_int),
        ("_file_handle", ctypes.c_void_p)
    ]

    def read_segment_header(self, index, segment):
        return int(_ovf_read_segment_header(ctypes.addressof(self), ctypes.c_int(index), ctypes.pointer(segment)))

    def read_segment_data(self, index, segment, data):
        if data.dtype == np.dtype('f'):
            datap = data.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            return int(_ovf_read_segment_data_4(ctypes.addressof(self), ctypes.c_int(index), ctypes.pointer(segment), datap))
        elif data.dtype == np.dtype('d'):
            datap = data.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            return int(_ovf_read_segment_data_8(ctypes.addressof(self), ctypes.c_int(index), ctypes.pointer(segment), datap))
        else:
            print("ovf.py read_segment_data: not able to use data type ", data.dtype)

    def write_segment(self, segment, data, fileformat=-56):
        if data.dtype == np.dtype('f'):
            datap = data.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            return int(_ovf_write_segment_4(ctypes.addressof(self), ctypes.pointer(segment), datap, fileformat))
        elif data.dtype == np.dtype('d'):
            datap = data.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            return int(_ovf_write_segment_8(ctypes.addressof(self), ctypes.pointer(segment), datap, fileformat))
        else:
            print("ovf.py read_segment_data: not able to use data type ", data.dtype)

    def append_segment(self, segment, data, fileformat=-56):
        if data.dtype == np.dtype('f'):
            datap = data.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            return int(_ovf_append_segment_4(ctypes.addressof(self), ctypes.pointer(segment), datap, fileformat))
        elif data.dtype == np.dtype('d'):
            datap = data.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            return int(_ovf_append_segment_8(ctypes.addressof(self), ctypes.pointer(segment), datap, fileformat))
        else:
            print("ovf.py read_segment_data: not able to use data type ", data.dtype)

    def get_latest_message(self):
        return str(_ovf_latest_message(ctypes.addressof(self)))


### Setup State
_ovf_open = _ovf.ovf_open
_ovf_open.argtypes = [ctypes.c_char_p]
_ovf_open.restype = ctypes.POINTER(_ovf_file)
def open(filename):
    return _ovf_open(ctypes.c_char_p(filename.encode('utf-8')))


### Delete State
_ovf_close = _ovf.ovf_close
_ovf_close.argtypes = [ctypes.POINTER(_ovf_file)]
_ovf_close.restype = ctypes.c_int
def close(p_file):
    return int(_ovf_close(p_file))


### ovf file wrapper class to be used in 'with' statement
class ovf_file():
    """Wrapper Class for an OVF file"""

    ### Functions to make 'with' statement work
    def __init__(self, filename):
        self._p_file = open(filename)

    def __enter__(self):
        return self._p_file.contents

    def __exit__(self, exc_type, exc_value, traceback):
        close(self._p_file)