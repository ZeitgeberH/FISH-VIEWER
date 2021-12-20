import numpy as np

from pyqtgraph.parametertree.parameterTypes import QtEnumParameter as enum
from pyqtgraph.Qt import QtWidgets

dlg = QtWidgets.QFileDialog

cfg = {
    'list': {
        'limits': {
            'type': 'checklist',
            'limits': ['All', 'Coronal', 'Sagittal']
        }
    },
    'file': {
        'acceptMode': {
            'type': 'list',
            'limits': list(enum(dlg.AcceptMode, dlg).enumMap)
        },
        'fileMode': {
            'type': 'list',
            'limits': list(enum(dlg.FileMode, dlg).enumMap)
        },
        'viewMode': {
            'type': 'list',
            'limits': list(enum(dlg.ViewMode, dlg).enumMap)
        },
        'dialogLabel': {
            'type': 'list',
            'limits': list(enum(dlg.DialogLabel, dlg).enumMap)
        },
        'relativeTo': {
            'type': 'str',
            'value': None
        },
        'directory': {
            'type': 'str',
            'value': None
        },
        'windowTitle': {
            'type': 'str',
            'value': None
        },
        'nameFilter': {
            'type': 'str',
            'value': None
        }
    },
    'float': {
        'Float Information': {
            'type': 'str',
            'readonly': True,
            'value': 'Note that all options except "finite" also apply to "int" parameters',
        },
        'step': {
            'type': 'float',
            'limits': [0, None],
            'value': 1,
        },
        'limits': {
            'type': 'list',
            'limits': {'[0, None]': [0, None], '[1, 5]': [1, 5]},
        },
        'suffix': {
            'type': 'list',
            'limits': ['Hz', 's', 'm'],
        },
        'siPrefix': {
            'type': 'bool',
            'value': True
        },
        'finite': {
            'type': 'bool',
            'value': True,
        },
        'dec': {
            'type': 'bool',
            'value': False,
        },
        'minStep': {
            'type': 'float',
            'value': 1.0e-12,
        },
    },

    'checklist': {
        'limits': {
            'type': 'checklist',
            'limits': ['one', 'two', 'three', 'four'],
        },
        'exclusive': {
            'type': 'bool',
            'value': False,
        }
    },


    'slider': {
        'step': {
            'type': 'float',
            'limits': [0, None],
            'value': 1, },
        'format': {
            'type': 'str',
            'value': '{0:>3}',
        },
        'precision': {
            'type': 'int',
            'value': 2,
            'limits': [1, None],
        },
        'span': {
            'type': 'list',
            'limits': {'linspace(-pi, pi)': np.linspace(-np.pi, np.pi), 'arange(10)**2': np.arange(10) ** 2},
        },

        'How to Set': {
            'type': 'list',
            'limits': ['Use span', 'Use step + limits'],
        }
    },

}
