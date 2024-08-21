import numpy as np

import node as n
import structs as s


def print_nodes(lst):
    for node in lst:
        node.Print_node()

class Element:

    def __init__(self, b_l, node2, t_l, node4, node5, node6, b_r, node8, t_r, index=None):
        self.corner_nodes = np.array((b_l, t_l, b_r, t_r))
        self.corner_nodes.sort()  # sorted to be b_l, t_l, b_r, t_r
        self.all_nodes = np.array((b_l, node2, t_l, node4, node5, node6, b_r, node8, t_r))
        self.all_nodes.sort()  # sorted bottom to top, left to right
        self.markers = list()
        self.index = index
        # node1 is bottom left node
        #print('set bottom left', end='')
        #node1.Print_node()
        b_l.Set_element(self)
        self.defining_node = self.corner_nodes[0]
        self.area = (b_r.x - b_l.x) * (t_l.y - b_l.y)
        self.mass = 0

    def Print_element(self):
        print("Nodes: ", end='')
        for node in self.corner_nodes:
            print(node.x, ",", node.y, sep='', end='; ')
        print()
        if (len(self.markers) != 0):
            print("Markers: ", end='')
            for marker in self.markers:
                print(marker.x, ",", marker.y, sep='', end='; ')
            print()

    def Add_marker(self, marker):
        self.markers.append(marker)

    def Remove_marker(self, marker):
        self.markers.remove(marker)

    def Set_mass(self, mass):
        self.mass = mass

    def Set_ip(self, ip):
        self.ip = ip

    # def sort_nodes(self, node_lst):
    #     self.nodes = list()
    #     self.nodes.append(n.Node(-1, -1))
    #     for node in range(4):
    #         old_len = len(self.nodes)
    #         k = old_len - 1
    #         while (len(self.nodes) == old_len):
    #             if (node_lst[0] > self.nodes[k]):
    #                 self.nodes.insert(k + 1, node_lst[0])
    #                 node_lst.pop(0)
    #             else:
    #                 k -= 1
    #     self.nodes.pop(0)


    # def __init__(self, b_l, t_l, b_r, t_r):
    #     self.corner_nodes = ([b_l, t_l, b_r, t_r])
    #     self.corner_nodes.sort()       # sorted to be b_l, t_l, b_r, t_r
