__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu ehsaneddin.asgari@helmholtz-hzi.de"

import random
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace, RectFace, random_color, ProfileFace
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches



class VisualizeCircularTree(object):
    '''
        Developed by Ehsaneddin Asgari
    '''
    def __init__(self, nwk_format):
        self.nwk=nwk_format


    def create_circle(self, filename, title, name2color=None, name2class_dic=None, class2color_dic=None, vector=None, ignore_branch_length=True):

        plt.clf()
        axis_font = {'size':'3'}
        plt.rc('xtick', labelsize=0.1)
        plt.rc('ytick', labelsize=0.1)
        plt.rc({'font.size':0.1})

        # legend creation
        if name2class_dic and class2color_dic:
            leg=[]
            for cls,color in class2color_dic.items():
                leg.append(mpatches.Patch(color=color, label=cls))

        t = Tree(self.nwk)
        # iterate over tree leaves only
        for l in t.iter_leaves():
            ns = NodeStyle()
            if name2color:
                ns["bgcolor"]=name2color[l.name] if l.name in name2color else 'white'
            elif name2class_dic and class2color_dic:
                ns["bgcolor"] = class2color_dic[name2class_dic[l.name]] if l.name in name2class_dic else 'white'
            # Gray dashed branch lines
            #ns["hz_line_type"] = 1
            #ns["hz_line_color"] = "#cccccc"
            #
            l.img_style = ns
            F=TextFace(l.name)
            F.ftype='Times'
            #if vector:
            #    if l.name in vector:
            l.add_features(profile = [random.random() for x in range(10)])#vector[l.name])
            l.add_features(deviation = [0 for x in range(10)])#len(vector[l.name]))])
            l.add_face(ProfileFace(max_v=1, min_v=0.0, center_v=0.5, width=200, height=40, style='heatmap', colorscheme=5), column=0, position='aligned')
        # Create an empty TreeStyle
        ts = TreeStyle()

        # Set our custom layout function
        ts.layout_fn = VisualizeCircularTree.layout

        # Draw a tree
        ts.mode = "c"

        # We will add node names manually
        ts.show_leaf_name = False
        # Show branch data
        ts.show_branch_length = True
        ts.show_branch_support = True
        ts.force_topology=ignore_branch_length
        ts.title.add_face(TextFace(title, fsize=20, ftype='Times'), column=15)

        # legend creation
        if name2class_dic and class2color_dic:
            for k , (cls, col) in enumerate(class2color_dic.items()):
                x=RectFace(8,8, 'black', col)
                #x.opacity=0.5
                ts.legend.add_face(x, column=8)
                ts.legend.add_face(TextFace(' '+cls+'   ', fsize=9,ftype='Times'), column=9)

        t.render(filename+'.pdf',tree_style=ts,dpi=5000)


    @staticmethod
    def layout(node):
        if node.is_leaf():
            # Add node name to laef nodes
            N = AttrFace("name", fsize=14, fgcolor="black")
            faces.add_face_to_node(N, node, 0)
        if "weight" in node.features:
            # Creates a sphere face whose size is proportional to node's
            # feature "weight"
            C = CircleFace(radius=node.weight, color="RoyalBlue", style="sphere")
            # Let's make the sphere transparent
            C.opacity = 0.3
            # And place as a float face over the tree
            faces.add_face_to_node(C, node, 0, position="float")

