# -*- coding: utf-8 -*-

from gr_allenFISH import GeneRenderAPI

## global setting in the file tab
CET_linearCmap = ['CET-L'+str(j) for j in range(1, 20)]
geapi_0 = GeneRenderAPI()
abAtlas_acronym= sorted(list(geapi_0.atlas_df['acronym']))
 
download_options = [
    {'name': 'Gene & Brain region & species', 'type':'group', 'children': [
        {'name': 'Gene', 'type': 'str', 'value': ''},
        {'name': 'Region', 'type': 'list', 'values':abAtlas_acronym,'value': 'ECT1'},
        {'name': 'Species','type': 'list', 'values': ['Mouse', 'Human'], 'value': 0},
        {'name': 'Synchronizing to atlas', 'type': 'bool', 'value':True},
        {'name': 'Sync', 'type': 'action'},
    ]
    },
    {'name': 'View options', 'type':'group', 'children': [
        {'name': 'Plane','type': 'list', 'values': ['All', 'Sagittal','Coronal'], 'value': 0},
        {'name': 'Expression', 'type': 'list', 'values':['color mask','approximated intensity'],'value': 'color mask'},
        {'name': 'Maximal resolution', 'type': 'bool', 'value': True},
        {'name': 'Image quality at lower resolution', 'type':'group', 'children': [
            {'name': 'Downsample', 'type': 'int', 'value': 4,\
            'limits': (0, 4),'default': 4, 'step': 1,'siPrefix': False},
            {'name': 'quality percentage', 'type': 'int', 'value': 50,\
            'limits': (0, 100),'default': 50, 'step': 10,'siPrefix': False},
                ]
        }
        ],
    },
    {'name': 'image-to-Atlas', 'type':'group', 'children': [
        {'name': 'Syn to atlas after double clicking', 'type': 'bool', 'value': True},
        {'name': 'atlas section offset', 'type': 'int', 'value': 0,\
             'limits': (-10, 10),'default': 0, 'step': 1,'siPrefix': False},       
        {'name': 'x offset', 'type': 'int', 'value': 0,\
             'limits': (-1000, 1000),'default': 0, 'step': 10,'siPrefix': False},
        {'name': 'y offset', 'type': 'int', 'value': 0,\
             'limits': (-1000, 1000),'default': 0, 'step': 10,'siPrefix': False},
        {'name': 'scale factor', 'type': 'float', 'value': 1,\
             'limits': (0.5, 2),'default': 1, 'step': 0.1,'siPrefix': False},
        # {'name': 'alpha ISH', 'type': 'float', 'value': 0.25,\
        #      'limits': (0, 1),'default': 0.2, 'step': 0.1,'siPrefix': False},
        {'name': 'Expression mask opacity', 'type': 'float', 'value': 0.65,\
             'limits': (0, 1),'default': 0.65, 'step': 0.05,'siPrefix': False},
    ]
    },
        
]