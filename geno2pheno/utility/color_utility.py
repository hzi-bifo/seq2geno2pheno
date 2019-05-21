__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu ehsaneddin.asgari@helmholtz-hzi.de"

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random


class ColorUtility(object):
    '''
        Developed by Ehsaneddin Asgari
    '''
    color_schemes = {'short':[ 'blue', 'red','green', 'gold', 'cyan'],
                 'large':['#ff0505', '#f2a041', '#cdff05', '#04d9cb', '#45a8ff', '#8503a6', '#590202', '#734d02', '#4ab304',
                  '#025359', '#0454cc', '#ff45da', '#993829', '#ffda45', '#1c661c', '#05cdff', '#1c2f66', '#731f57',
                  '#b24a04', '#778003', '#0e3322', '#024566', '#0404d9', '#e5057d', '#66391c', '#31330e', '#3ee697',
                  '#2d7da6', '#20024d', '#33011c'] + list(({'aliceblue': '#F0F8FF', 'antiquewhite': '#FAEBD7',
                                                            'aqua': '#00FFFF', 'aquamarine': '#7FFFD4',
                                                            'azure': '#F0FFFF', 'beige': '#F5F5DC', 'bisque': '#FFE4C4',
                                                            'black': '#000000', 'blanchedalmond': '#FFEBCD',
                                                            'blue': '#0000FF', 'blueviolet': '#8A2BE2',
                                                            'brown': '#A52A2A', 'burlywood': '#DEB887',
                                                            'cadetblue': '#5F9EA0', 'chartreuse': '#7FFF00',
                                                            'chocolate': '#D2691E', 'coral': '#FF7F50',
                                                            'cornflowerblue': '#6495ED', 'cornsilk': '#FFF8DC',
                                                            'crimson': '#DC143C', 'cyan': '#00FFFF',
                                                            'darkblue': '#00008B', 'darkcyan': '#008B8B',
                                                            'darkgoldenrod': '#B8860B', 'darkgray': '#A9A9A9',
                                                            'darkgreen': '#006400', 'darkkhaki': '#BDB76B',
                                                            'darkmagenta': '#8B008B', 'darkolivegreen': '#556B2F',
                                                            'darkorange': '#FF8C00', 'darkorchid': '#9932CC',
                                                            'darkred': '#8B0000', 'darksalmon': '#E9967A',
                                                            'darkseagreen': '#8FBC8F', 'darkslateblue': '#483D8B',
                                                            'darkslategray': '#2F4F4F', 'darkturquoise': '#00CED1',
                                                            'darkviolet': '#9400D3', 'deeppink': '#FF1493',
                                                            'deepskyblue': '#00BFFF', 'dimgray': '#696969',
                                                            'dodgerblue': '#1E90FF', 'firebrick': '#B22222',
                                                            'floralwhite': '#FFFAF0', 'forestgreen': '#228B22',
                                                            'fuchsia': '#FF00FF', 'gainsboro': '#DCDCDC',
                                                            'ghostwhite': '#F8F8FF', 'gold': '#FFD700',
                                                            'goldenrod': '#DAA520', 'gray': '#808080',
                                                            'green': '#008000', 'greenyellow': '#ADFF2F',
                                                            'honeydew': '#F0FFF0', 'hotpink': '#FF69B4',
                                                            'indianred': '#CD5C5C', 'indigo': '#4B0082',
                                                            'ivory': '#FFFFF0', 'khaki': '#F0E68C',
                                                            'lavender': '#E6E6FA', 'lavenderblush': '#FFF0F5',
                                                            'lawngreen': '#7CFC00', 'lemonchiffon': '#FFFACD',
                                                            'lightblue': '#ADD8E6', 'lightcoral': '#F08080',
                                                            'lightcyan': '#E0FFFF', 'lightgoldenrodyellow': '#FAFAD2',
                                                            'lightgreen': '#90EE90', 'lightgray': '#D3D3D3',
                                                            'lightpink': '#FFB6C1', 'lightsalmon': '#FFA07A',
                                                            'lightseagreen': '#20B2AA', 'lightskyblue': '#87CEFA',
                                                            'lightslategray': '#778899', 'lightsteelblue': '#B0C4DE',
                                                            'lightyellow': '#FFFFE0', 'lime': '#00FF00',
                                                            'limegreen': '#32CD32', 'linen': '#FAF0E6',
                                                            'magenta': '#FF00FF', 'maroon': '#800000',
                                                            'mediumaquamarine': '#66CDAA', 'mediumblue': '#0000CD',
                                                            'mediumorchid': '#BA55D3', 'mediumpurple': '#9370DB',
                                                            'mediumseagreen': '#3CB371', 'mediumslateblue': '#7B68EE',
                                                            'mediumspringgreen': '#00FA9A',
                                                            'mediumturquoise': '#48D1CC', 'mediumvioletred': '#C71585',
                                                            'midnightblue': '#191970', 'mintcream': '#F5FFFA',
                                                            'mistyrose': '#FFE4E1', 'moccasin': '#FFE4B5',
                                                            'navajowhite': '#FFDEAD', 'navy': '#000080',
                                                            'oldlace': '#FDF5E6', 'olive': '#808000',
                                                            'olivedrab': '#6B8E23', 'orange': '#FFA500',
                                                            'orangered': '#FF4500', 'orchid': '#DA70D6',
                                                            'palegoldenrod': '#EEE8AA', 'palegreen': '#98FB98',
                                                            'paleturquoise': '#AFEEEE', 'palevioletred': '#DB7093',
                                                            'papayawhip': '#FFEFD5', 'peachpuff': '#FFDAB9',
                                                            'peru': '#CD853F', 'pink': '#FFC0CB', 'plum': '#DDA0DD',
                                                            'powderblue': '#B0E0E6', 'purple': '#800080',
                                                            'red': '#FF0000', 'rosybrown': '#BC8F8F',
                                                            'royalblue': '#4169E1', 'saddlebrown': '#8B4513',
                                                            'salmon': '#FA8072', 'sandybrown': '#FAA460',
                                                            'seagreen': '#2E8B57', 'seashell': '#FFF5EE',
                                                            'sienna': '#A0522D', 'silver': '#C0C0C0',
                                                            'skyblue': '#87CEEB', 'slateblue': '#6A5ACD',
                                                            'slategray': '#708090', 'snow': '#FFFAFA',
                                                            'springgreen': '#00FF7F', 'steelblue': '#4682B4',
                                                            'tan': '#D2B48C', 'teal': '#008080', 'thistle': '#D8BFD8',
                                                            'tomato': '#FF6347', 'turquoise': '#40E0D0',
                                                            'violet': '#EE82EE', 'wheat': '#F5DEB3', 'white': '#FFFFFF',
                                                            'whitesmoke': '#F5F5F5', 'yellow': '#FFFF00',
                                                            'yellowgreen': '#9ACD32'}).keys()),
                 'medium':['#ff0505', '#f2a041', '#cdff05', '#04d9cb', '#45a8ff', '#8503a6', '#590202', '#734d02', '#4ab304',
                  '#025359', '#0454cc', '#ff45da', '#993829', '#ffda45', '#1c661c', '#05cdff', '#1c2f66', '#731f57',
                  '#b24a04', '#778003', '#0e3322', '#024566', '#0404d9', '#e5057d', '#66391c', '#31330e', '#3ee697',
                  '#2d7da6', '#20024d', '#33011c']}
    fg_color_map=['#e5bd89' , '#f8e076' , '#da4f38' , '#ffde3b' ,'#990f02' , '#fe7d68' ,'#66042d' , '#301131' ,'#281e5d' , '#042c36' , '#52b2c0' , '#74b62e' , '#32612d' ,'#522915' , '#2d1606' , '#7e7d9c' , '#f59bbf' , '#ffc30b' , '#fffada' , '#01381b' , '#e96210' , '#6c2d7e' , '#050100' , '#5d2c04' , '#6a6880']
    bg_color_map=['#d5ba9c' , '#b79266' , '#9b7b55' , '#c97f80' , '#b75556' , '#ae4041' , '#fe664e' , '#cc4a34' , '#cf1d01' , '#ff9934' , '#e57100' , '#ff7f00' , '#d8d828' , '#ffff01' , '#ffff01' , '#355f3b' , '#4b6e50' , '#738f76' , '#aeadcd' , '#a4a3cd' , '#7b7a9a' , '#b28e98' , '#e5b6c6' , '#ffcada']

    def __init__(self):
        print ('Color Utility v 1.0.0')

    @staticmethod
    def gen_hex_colour_code():
        return '#'+''.join([random.choice('0123456789ABCDEF') for x in range(6)])





