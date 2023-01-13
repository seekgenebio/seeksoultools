import os

__srcdir = os.path.dirname(os.path.abspath(__file__))

ADAPTERS = ["AAAAAAAAAAAA", ]
MINLEN = 60

CHEMISTRY = {
    '__SO01V3':{
        'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8L8B8L10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
        'linker': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'Linker1.txt'),
                   os.path.join(__srcdir, 'barcode', 'SO01V3', 'Linker2.txt'),),
    },
    '__nolinker':{
        'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
    },
    "__P3CBGB":{
        'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
    },
    "DDV1":{
       'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
    },
    "DDV2":{
       'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
    },
    "DD5V1":{
       'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
        'sc5p': True,
    },
    "MM":{
       'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
    },
    "MM-D":{
       'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
    },
    "DD-Q":{
       'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
    }
}

