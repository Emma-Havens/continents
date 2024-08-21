import math
import random
import time
import warnings
warnings.filterwarnings("error")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy

import element as e
import node as n
import marker as m
import structs as s


class Space:

    def __init__(self, lenRow, lenCol, width, height):
        self.width = width
        self.height = height
        self.lenRow = lenRow
        self.lenCol = lenCol
        self.num_x_nodes = lenRow * 2 + 1  # num nodes in x direction
        self.num_y_nodes = lenCol * 2 + 1  # num nodes in y direction
        self.delta_x = width / lenRow
        self.delta_y = height / lenCol

        s.IDNum.reset()

        self.nodes = np.full(self.num_x_nodes * self.num_y_nodes, None)
        i = 0
        x_scalar = self.width / (self.num_x_nodes - 1)
        y_scalar = self.height / (self.num_y_nodes - 1)
        x_start = self.width / -2
        y_start = self.height / -2
        for x in range(self.num_x_nodes):
            for y in range(self.num_y_nodes):
                # self.nodes[i] = n.Node(x, y, i)
                self.nodes[i] = n.Node(x_start + x * x_scalar, y_start + y * y_scalar, i)
                i += 1

        self.x_min = self.nodes[0].x
        self.x_max = self.nodes[self.num_x_nodes * self.num_y_nodes - 1].x
        self.y_min = self.nodes[0].y
        self.y_max = self.nodes[self.num_y_nodes - 1].y

        self.elements = np.full(lenRow * lenCol, None)
        self.defining_nodes = np.full(lenRow * lenCol, None)
        i = 0
        for y in range(lenCol):
            for x in range(lenRow):
                # nodes for a grid
                b_l = self.nodes[2 * x * self.num_y_nodes + 2 * y]
                node2 = self.nodes[(2 * x * self.num_y_nodes + 2 * y) + 1]
                t_l = self.nodes[(2 * x * self.num_y_nodes + 2 * y) + 2]
                node4 = self.nodes[(2 * (x + 1) - 1) * self.num_y_nodes + 2 * y]
                node5 = self.nodes[((2 * (x + 1) - 1) * self.num_y_nodes + 2 * y) + 1]
                node6 = self.nodes[((2 * (x + 1) - 1) * self.num_y_nodes + 2 * y) + 2]
                b_r = self.nodes[2 * (x + 1) * self.num_y_nodes + 2 * y]
                node8 = self.nodes[(2 * (x + 1) * self.num_y_nodes + 2 * y) + 1]
                t_r = self.nodes[(2 * (x + 1) * self.num_y_nodes + 2 * y) + 2]
                self.elements[i] = e.Element(b_l, node2, t_l, node4, node5, node6, b_r, node8, t_r, i)
                self.defining_nodes[i] = b_l
                i += 1
        self.defining_nodes.sort()
        # self.print_node_list(self.defining_nodes)
        # self.print_nodes_and_elements()

    #region MARKER CREATION AND LOCATION

    def Create_markers(self, markerMultiplier=1):
        numMarkers = int(self.lenRow * self.lenCol * markerMultiplier)
        self.markers = np.full(numMarkers, None)
        for i in range(len(self.markers)):
            x = random.uniform(self.x_min, self.x_max)
            y = random.uniform(self.y_min, self.y_max)
            self.markers[i] = m.Marker(round(x, 2), round(y, 2),
                                       self.x_min, self.x_max, self.y_min, self.y_max)
        # self.print_markers()

        self.Locate_markers()

    def Add_marker(self):
        self.markers = np.full(1, None)
        x = self.x_max / 2
        y = 0
        self.markers[0] = m.Marker(round(x, 2), round(y, 2),
                                   self.x_min, self.x_max, self.y_min, self.y_max)
        self.Locate_markers()

    def Node_index_of_marker(self, marker):
        # index = self.binary_search(marker) - 1
        # print('locating marker ', end='')
        # marker.Print_marker()
        # print('at index', index, ': ', end='')
        # self.nodes[index].Print_node()
        # return index
        return self.binary_search(marker) - 1

    def Node_index_given_index(self, marker, index):
        pass

    def binary_search(self, marker):
        only_x = lambda pos : pos.x
        only_y = lambda pos : pos.y
        end_x = self.my_bisect_right(self.defining_nodes, marker, only_x)
        start_x = end_x - self.lenCol
        return self.my_bisect_right(self.defining_nodes, marker, only_y, start_x, end_x)

    def my_bisect_right(self, node_lst, marker, key, lo=0, hi=None):
        if hi is None:
            hi = len(node_lst)

        while lo < hi:
            mid = (lo + hi) // 2
            if key(marker) < key(node_lst[mid]):
                hi = mid
            else:
                lo = mid + 1
        return lo

    def fix_node_index(self, index, marker):
        new_index = index
        if (marker.y == self.nodes[index].y):
            new_index = new_index - 1
        if (not self.nodes[new_index].defining and marker.x == self.nodes[index].x):
            new_index = new_index - self.num_y_nodes
        if (not self.nodes[new_index].defining):
            new_index = new_index + 1
            if (not self.nodes[new_index].defining):
                raise ValueError("Cannot find closest defining node to (",
                                 marker.x, marker.y, ") came up with ", new_index)
        return new_index

    def Locate_markers(self):
        for marker in self.markers:
            index = self.Node_index_of_marker(marker)
            node = self.defining_nodes[index]
            try:
                element = node.bottom_left_of
            except AttributeError:
                print('trying to fix Node index')
                # print('old index:', index)
                index = self.fix_node_index(index, marker)
                # print('new index:', index)
                node = self.defining_nodes[index]
                element = node.bottom_left_of
            marker.Place_in_element(element)
            element.Add_marker(marker)
        # self.print_markers()

    def Relocate_marker(self, marker, moved=False):
        # prev_index = marker.container.defining_node.index
        index = self.Node_index_of_marker(marker)
        node = self.defining_nodes[index]
        try:
            element = node.bottom_left_of   # if this Errors, self.Move_markers handles it
        except AttributeError:
            if moved:
                print('trying to fix Node index')
                # print('old index:', index)
                index = self.fix_node_index(index, marker)
                # print('new index:', index)
                node = self.defining_nodes[index]
                element = node.bottom_left_of
            else:
                raise AttributeError
        marker.container.Remove_marker(marker)
        marker.Place_in_element(element)
        element.Add_marker(marker)

    #endregion

    #region VELOCITY

    def Set_velocity_eq(self, vel):     # lambda equation, must take in 2 var, must return an s.Vel
        self.velEq = vel
        for node in self.nodes:
            node.Set_Vel(self.velEq(node.x, node.y))

    def Set_advection(self, magnitude, n):
        L_x = self.x_max
        L_y = self.y_max
        # self.print_nodes()
        # print(L_x, ',', L_y)
        self.velEq = lambda x, y : s.Vel(magnitude * np.cos((n * np.pi * x) / (2 * L_x)) * np.sin((n * np.pi * y) / (2 * L_y)),
                                  -magnitude * np.sin((n * np.pi * x) / (2 * L_x)) * np.cos((n * np.pi * y) / (2 * L_y)))
        for node in self.nodes:
            node.Set_Vel(self.velEq(node.x, node.y))

    def Interpolate(self, left_v, left, right_v, right, p):
        return left_v * (right - p) / (right - left) + \
                right_v * (p - left) / (right - left)

    def Bilinear_interpolation_vel(self, b_l, b_r, t_l, t_r, p):
        x_vel_below = self.Interpolate(b_l.vel.v_x, b_l.x, b_r.vel.v_x, b_r.x, p.x)
        x_vel_above = self.Interpolate(t_l.vel.v_x, t_l.x, t_r.vel.v_x, t_r.x, p.x)
        x_vel = self.Interpolate(x_vel_below, b_l.y, x_vel_above, t_l.y, p.y)
        y_vel_left = self.Interpolate(b_l.vel.v_y, b_l.y, t_l.vel.v_y, t_l.y, p.y)
        y_vel_right = self.Interpolate(b_r.vel.v_y, b_r.y, t_r.vel.v_y, t_r.y, p.y)
        y_vel = self.Interpolate(y_vel_left, b_l.x, y_vel_right, b_r.x, p.x)
        return s.Vel(x_vel, y_vel)

    def Set_marker_velocities(self):
        for marker in self.markers:
            n = marker.container.corner_nodes      # sorted to be b_l, t_l, b_r, t_r
            marker.Set_Vel(self.Bilinear_interpolation_vel(n[0], n[2], n[1], n[3], marker))
        # self.print_markers()

    def Move_markers(self, step):
        for marker in self.markers:
            try:
                marker.Move_marker(step)
                self.Relocate_marker(marker)    # raises AttributeError if node is not defining
            except AttributeError:
                # print('keeping in bounds')
                # print('old marker: ', end='')
                # marker.Print_marker()
                marker.Keep_in_bounds()
                # print('new marker: ', end='')
                # marker.Print_marker()
                self.Relocate_marker(marker, True)

    #endregion

    #region INTEGRATION FUNCTIONS
    def Set_density_at_nodes(self, vel):     # lambda equation, must take in 2 var, must return a value
        self.densityEq = vel
        for node in self.nodes:
            node.Set_Density(self.densityEq(node.x, node.y))

    def Set_density(self, density):         # lambda equation, must take in 2 var, must return a value
        self.densityEq = density
        try:
            for marker in self.markers:
                marker.Set_Density(self.densityEq(marker.x, marker.y))
        except AttributeError:
            raise AttributeError("Markers were not created before density was assigned")

    def Update_Node_density(self):
        num_markers_per_node = np.zeros((len(self.nodes), 1))
        sum_per_node = np.zeros((len(self.nodes), 1))
        for element in self.elements:
            n = element.corner_nodes   # sorted to be b_l, t_l, b_r, t_r
            o_x = n[0].x + (n[2].x - n[0].x) / 2
            o_y = n[0].y + (n[1].y - n[0].y) / 2
            # print(n[2].x)
            # print(o_x, o_y)
            for marker in element.markers:
                if marker.x < o_x:
                    if marker.y < o_y:
                        node = n[0]     # p_11
                    else:
                        node = n[1]     # p_12
                else:
                    if marker.y < o_y:
                        node = n[2]     # p_21
                    else:
                        node = n[3]     # p_22
                # if n[0].x <= marker.x < o_x:
                #     if n[0].y <= marker.y < o_y:
                #         node = n[0]     # p_11
                #     elif o_y <= marker.y < n[1].y:
                #         node = n[1]     # p_12
                # elif o_x <= marker.x <= n[2].x:
                #     if n[0].y <= marker.y < o_y:
                #         node = n[2]     # p_21
                #     elif o_y <= marker.y <= n[1].y:
                #         node = n[3]     # p_22
                # else:
                #     print('UPDATE NODE DENSITY')
                #     marker.Print_marker()
                #     print('is not contained in')
                #     element.Print_element()
                #     print('at time step', time)

                num_markers_per_node[node.index] += 1
                sum_per_node[node.index] += marker.density

        for i in range(len(self.nodes)):
            try:
                self.nodes[i].density = np.divide(sum_per_node[i], num_markers_per_node[i])
                # print(self.nodes[i].density)
            except:
                #self.nodes[i].density = -1
                #print('special', self.nodes[i].density)
                raise ValueError('No Markers are within distance to Node', i)

    def Bilinear_interpolation_density(self, b_l, b_r, t_l, t_r, p):
        below = self.Interpolate(b_l.density, b_l.x, b_r.density, b_r.x, p.x)
        above = self.Interpolate(t_l.density, t_l.x, t_r.density, t_r.x, p.x)
        density = self.Interpolate(below, b_l.y, above, t_l.y, p.y)
        return density

    def gauss_quad_solver_no_interp(self, f, a, b, w):
        areas = np.zeros((len(self.elements), 1))
        loc_len_x = 2
        loc_len_y = 2
        for e in range(len(self.elements)):
            n = self.elements[e].corner_nodes   # sorted to be b_l, t_l, b_r, t_r
            self.print_node_list(n)
            x_11, x_12, x_21, x_22 = n[0].x, n[1].x, n[2].x, n[3].x
            y_11, y_12, y_21, y_22 = n[0].y, n[1].y, n[2].y, n[3].y
            ip_sum = 0
            len_ratio_x = (x_21 - x_11) / loc_len_x
            len_ratio_y = (y_12 - y_11) / loc_len_y
            e_or_x = x_11 + (x_21 - x_11) / 2
            e_or_y = y_11 + (y_12 - y_11) / 2
            for i in range(len(a)):
                x_ip = e_or_x + a[i] * len_ratio_x
                y_ip = e_or_y + b[i] * len_ratio_y
                print(x_ip, ',', y_ip)
                ip_sum += f(x_ip, y_ip) * w[i]  * len_ratio_x * len_ratio_y
                print(e,':', ip_sum)
            areas[e] = ip_sum
        return areas
    
    def gauss_quad_solver(self, a, b, w):
        areas = np.zeros((len(self.elements), 1))
        loc_len_x = 2
        loc_len_y = 2
        for e in range(len(self.elements)):
            n = self.elements[e].corner_nodes   # sorted to be b_l, t_l, b_r, t_r
            self.print_node_list(n)
            x_11, x_12, x_21, x_22 = n[0].x, n[1].x, n[2].x, n[3].x
            y_11, y_12, y_21, y_22 = n[0].y, n[1].y, n[2].y, n[3].y
            ip_sum = 0
            len_ratio_x = (x_21 - x_11) / loc_len_x
            len_ratio_y = (y_12 - y_11) / loc_len_y
            e_or_x = x_11 + (x_21 - x_11) / 2
            e_or_y = y_11 + (y_12 - y_11) / 2
            for i in range(len(a)):
                x_ip = e_or_x + a[i] * len_ratio_x
                y_ip = e_or_y + b[i] * len_ratio_y
                p = m.Marker(x_ip, y_ip)
                print(x_ip, ',', y_ip)
                ip_sum += self.Bilinear_interpolation_density(n[0], n[2], n[1], n[3], p) * w[i] * len_ratio_x * len_ratio_y
                print(e,':', ip_sum)
            areas[e] = ip_sum
        return areas

    def gauss_quad_no_interp(self, f):
        a = [-1/np.sqrt(3), -1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]
        b = [-1 / np.sqrt(3), 1 / np.sqrt(3), -1 / np.sqrt(3), 1 / np.sqrt(3)]
        w = [1, 1, 1, 1]
        self.gauss_area = self.gauss_quad_solver_no_interp(f, a, b, w)
        
    def Gauss_quad(self):
        a = [-1/np.sqrt(3), -1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]
        b = [-1 / np.sqrt(3), 1 / np.sqrt(3), -1 / np.sqrt(3), 1 / np.sqrt(3)]
        w = [1, 1, 1, 1]
        self.gauss_area = self.gauss_quad_solver(a, b, w)

    def Int_area(self, f):
        self.int_area = np.zeros((len(self.elements), 1))
        for e in range(len(self.elements)):
            n = self.elements[e].corner_nodes  # sorted to be b_l, t_l, b_r, t_r
            x_1, x_2 = n[0].x, n[2].x
            y_1, y_2 = n[0].y, n[1].y
            density = f(x_1, x_2, y_1, y_2)
            self.int_area[e] = density

    def Integrate_Element_Density(self, element, a, b, w):
        loc_len_x = 2
        loc_len_y = 2

        n = element.corner_nodes   # sorted to be b_l, t_l, b_r, t_r
        x_11, x_12, x_21, x_22 = n[0].x, n[1].x, n[2].x, n[3].x
        y_11, y_12, y_21, y_22 = n[0].y, n[1].y, n[2].y, n[3].y
        ip_sum = 0
        len_ratio_x = (x_21 - x_11) / loc_len_x
        len_ratio_y = (y_12 - y_11) / loc_len_y
        e_or_x = x_11 + (x_21 - x_11) / 2
        e_or_y = y_11 + (y_12 - y_11) / 2
        for i in range(len(a)):
            x_ip = e_or_x + a[i] * len_ratio_x
            y_ip = e_or_y + b[i] * len_ratio_y
            p = m.Marker(x_ip, y_ip)
            ip_sum += self.Bilinear_interpolation_density(n[0], n[2], n[1], n[3], p) * w[i] * len_ratio_x * len_ratio_y
        return ip_sum

    def Update_Element_mass(self):
        a = [-1 / np.sqrt(3), -1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3)]
        b = [-1 / np.sqrt(3), 1 / np.sqrt(3), -1 / np.sqrt(3), 1 / np.sqrt(3)]
        w = [1, 1, 1, 1]
        total_mass = 0
        for element in self.elements:
            mass = self.Integrate_Element_Density(element, a, b, w)
            element.Set_mass(mass)
            total_mass += mass
        # print('Total mass:', total_mass)
        return total_mass

    #endregion

    #region STOKES SOLVER
    def Local_ip_arr(self, ip_x, ip_y, shp_func):
        ip_val = np.zeros((len(ip_x), len(shp_func)), dtype=np.float32)
        for i in range(len(ip_x)):
            for j in range(len(shp_func)):
                ip_val[i][j] = shp_func[j](ip_x[i], ip_y[i])
        return ip_val

    def Set_local_ip(self):
        ip_x_dir = np.array([-1 * np.sqrt((3 + 2 * np.sqrt(6 / 5)) / 7), -1 * np.sqrt((3 - 2 * np.sqrt(6 / 5)) / 7),
                             np.sqrt((3 - 2 * np.sqrt(6 / 5)) / 7), np.sqrt((3 + 2 * np.sqrt(6 / 5)) / 7)])
        ip_y_dir = ip_x_dir
        ip_y, ip_x = np.meshgrid(ip_x_dir, ip_y_dir)
        ip_x = ip_x.reshape((1, (len(ip_x_dir) * len(ip_y_dir))))[0]
        ip_y = ip_y.reshape((1, (len(ip_x_dir) * len(ip_y_dir))))[0]

        # N1X = 0.5 * eta1 * (eta1 - 1)     = 0 iff eta1 = -1
        # N2X = -(eta1 + 1) * (eta1 - 1)    = 0 iff eta1 = 0
        # N3X = 0.5 * eta1 * (eta1 + 1)     = 0 iff eta1 = 1
        # N1Y = 0.5 * eta2 * (eta2 - 1)     = 0 iff eta2 = -1
        # N2Y = -(eta2 + 1) * (eta2 - 1)    = 0 iff eta2 = 0
        # N3Y = 0.5 * eta2 * (eta2 + 1)     = 0 iff eta2 = 1

        shp = [lambda x, y: 0.5 * x * (x - 1) * 0.5 * y * (y - 1),  # -1, -1
               lambda x, y: 0.5 * x * (x - 1) * -1 * (y + 1) * (y - 1),  # -1, 0
               lambda x, y: 0.5 * x * (x - 1) * 0.5 * y * (y + 1),  # -1, 1
               lambda x, y: -1 * (x + 1) * (x - 1) * 0.5 * y * (y - 1),  # 0, -1
               lambda x, y: -1 * (x + 1) * (x - 1) * -1 * (y + 1) * (y - 1),  # 0, 0
               lambda x, y: -1 * (x + 1) * (x - 1) * 0.5 * y * (y + 1),  # 0, 1
               lambda x, y: 0.5 * x * (x + 1) * 0.5 * y * (y - 1),  # 1, -1
               lambda x, y: 0.5 * x * (x + 1) * -1 * (y + 1) * (y - 1),  # 1, 0
               lambda x, y: 0.5 * x * (x + 1) * 0.5 * y * (y + 1)]  # 1, 1

        shp_deriv_wrt_x = [lambda x, y: (x - 0.5) * 0.5 * y * (y - 1),  # -1, -1
               lambda x, y: (x - 0.5) * -1 * (y + 1) * (y - 1),  # -1, 0
               lambda x, y: (x - 0.5) * 0.5 * y * (y + 1),  # -1, 1
               lambda x, y: -2 * x * 0.5 * y * (y - 1),  # 0, -1
               lambda x, y: -2 * x * -1 * (y + 1) * (y - 1),  # 0, 0
               lambda x, y: -2 * x * 0.5 * y * (y + 1),  # 0, 1
               lambda x, y: (x + 0.5) * 0.5 * y * (y - 1),  # 1, -1
               lambda x, y: (x + 0.5) * -1 * (y + 1) * (y - 1),  # 1, 0
               lambda x, y: (x + 0.5) * 0.5 * y * (y + 1)]  # 1, 1

        shp_deriv_wrt_y = [lambda x, y: 0.5 * x * (x - 1) * (y - 0.5),  # -1, -1
               lambda x, y: 0.5 * x * (x - 1) * -2 * y,  # -1, 0
               lambda x, y: 0.5 * x * (x - 1) * (y + 0.5),  # -1, 1
               lambda x, y: -1 * (x + 1) * (x - 1) * (y - 0.5),  # 0, -1
               lambda x, y: -1 * (x + 1) * (x - 1) * -2 * y,  # 0, 0
               lambda x, y: -1 * (x + 1) * (x - 1) * (y + 0.5),  # 0, 1
               lambda x, y: 0.5 * x * (x + 1) * (y - 0.5),  # 1, -1
               lambda x, y: 0.5 * x * (x + 1) * -2 * y,  # 1, 0
               lambda x, y: 0.5 * x * (x + 1) * (y + 0.5)]  # 1, 1

        self.shp_at_ip = self.Local_ip_arr(ip_x, ip_y, shp)
        self.shp_deriv = np.array((self.Local_ip_arr(ip_x, ip_y, shp_deriv_wrt_x),
                                  self.Local_ip_arr(ip_x, ip_y, shp_deriv_wrt_y)))

        ip_weights = np.array([(18 - np.sqrt(30)) / 36, (18 + np.sqrt(30)) / 36,
                               (18 + np.sqrt(30)) / 36, (18 - np.sqrt(30)) / 36])
        self.ip_weights = np.ndarray.flatten(np.outer(np.transpose(ip_weights), ip_weights))
        # print('trans', np.transpose(np.atleast_2d(ip_weights)))
        # print('not trans:', ip_weights)
        # print('ip_weights:', self.ip_weights)


    def Set_element_ip(self):
        self.Set_local_ip()
        # print(self.shp_at_ip)
        # print('len:', len(self.shp_at_ip))

        for element in self.elements:
            element_ip = np.full(len(self.shp_at_ip), None)
            # print('element_ip:', element_ip)
            # print('element_ip', len(element_ip))
            # print('all nodes')
            # self.print_node_list(element.all_nodes)

            # self.local_ip must traverse nodes in same way as element.all_nodes
            for i in range(len(self.shp_at_ip)):
                x_coor = 0
                y_coor = 0
                for j in range(len(self.shp_at_ip[0])):
                    cur_node = element.all_nodes[j]
                    x_coor += cur_node.x * self.shp_at_ip[i][j]
                    y_coor += cur_node.y * self.shp_at_ip[i][j]
                element_ip[i] = m.Marker(x_coor, y_coor)

            element.Set_ip(element_ip)

    def Set_Phases(self, density=None):
        if not density:
            density = np.array([20, 10])
        self.density = density
        # one_row = np.ones(self.lenCol, dtype=int)
        # bl = math.floor(.7 * self.lenCol)
        # one_row[bl::] = 0
        # self.phases = np.tile(one_row, self.lenRow)

        self.phases = np.ones(self.elements.shape, dtype=int)
        index = math.floor(0.8 * self.lenCol) * self.lenRow
        self.phases[index::] = 0

    def Stokes_solver(self, start_time):
        Ra = 1
        D = np.array([[4/3, -2/3, 0],
                     [-2/3, 4/3, 0],
                     [0, 0, 1]])
        G = np.array([[0], [-1]])
        Mu = 1
        sdof = 2 * len(self.nodes)
        nel = len(self.elements)
        A_all = np.zeros((sdof, sdof))
        Q_all = np.zeros((sdof, 3*nel))
        invM_all = np.zeros((3*nel, 3*nel))
        Rhs_all = np.zeros((sdof, 1))
        print('initialize globals', time.time() - start_time)
        cur_time = start_time

        for i in range(len(self.elements)):
            element = self.elements[i]
            rho = self.density[self.phases[i]]
            n_coor = np.array((np.fromiter((n.x for n in element.all_nodes), dtype=float),
                              np.fromiter((n.y for n in element.all_nodes), dtype=float)))
            n_coor = np.transpose(n_coor)
            A_elem = np.zeros((18, 18))
            Q_elem = np.zeros((18, 3))
            M_elem = np.zeros((3, 3))
            Rhs_elem = np.zeros((2, 9))
            P = np.array([[1, 1, 1],
                         [n_coor[0, 0], n_coor[6, 0], n_coor[8, 0]],
                         [n_coor[0, 1], n_coor[6, 1], n_coor[8, 1]]])
            P_inv = np.linalg.inv(P)
            print('initialize element', element.index, time.time() - cur_time)
            cur_time = time.time()

            for i in range(len(self.shp_at_ip)):
                ip_shp_deriv = np.array((self.shp_deriv[0][i], self.shp_deriv[1][i]))
                Pb = np.matmul(self.shp_at_ip[i], n_coor)
                Pb = np.concatenate((np.array([1]), Pb), axis=None)
                Pb = np.transpose(np.atleast_2d(Pb))
                Pi = np.matmul(P_inv, Pb)
                Pi_trans = np.transpose(np.atleast_2d(Pi))
                J = np.matmul(ip_shp_deriv, n_coor)
                invJ = np.linalg.inv(J)
                detJ = np.linalg.det(J)
                dNdx = np.matmul(np.transpose(invJ), ip_shp_deriv)
                B = np.zeros((3, len(ip_shp_deriv[0]) * 2))
                B_vol = np.zeros((len(ip_shp_deriv[0]) * 2, 1))
                weight = self.ip_weights[i] * detJ
                # print('ip', i, time.time() - cur_time)
                # cur_time = time.time()

                B[0,::2] = dNdx[0]
                B[1, 1::2] = dNdx[1]
                B[2, ::2] = dNdx[1]
                B[2, 1::2] = dNdx[0]
                B_vol[::2, 0] = dNdx[0]
                B_vol[1::2, 0] = dNdx[1]
                # print('B', time.time() - cur_time)
                # cur_time = time.time()

                A_elem = A_elem + weight * Mu * np.matmul(np.transpose(B), np.matmul(D, B))
                Q_elem = Q_elem - weight * np.matmul(B_vol, Pi_trans)
                M_elem = M_elem + weight * np.matmul(Pi, Pi_trans)

                Rhs_elem[0] = Rhs_elem[0] + weight * rho * G[0] * self.shp_at_ip[i]   # element density matters here
                Rhs_elem[1] = Rhs_elem[1] + weight * rho * G[1] * self.shp_at_ip[i]
                # print('AQMRhs', time.time() - cur_time)
                # cur_time = time.time()

            invM_elem = np.linalg.inv(M_elem)
            A_elem = A_elem + 1000 * np.matmul(np.matmul(Q_elem, invM_elem), np.transpose(Q_elem))

            adjusted_index = lambda n: n.index * 2
            global_nodes = np.array([adjusted_index(n) for n in element.all_nodes], dtype=int)
            A_index = np.zeros(len(global_nodes) * 2, dtype=int)
            A_index[::2] = global_nodes[:]
            A_index[1::2] = A_index[::2] + 1
            e_index = element.index
            Q_index = np.array([e_index, e_index + 1, e_index + 2], dtype=int)

            A_all[np.ix_(A_index, A_index)] += A_elem
            Q_all[np.ix_(A_index, Q_index)] += Q_elem
            invM_all[np.ix_(Q_index, Q_index)] += invM_elem
            Rhs_all[A_index] += np.transpose(np.atleast_2d(Rhs_elem.flatten('F')))
            # print('AQMRhs global', time.time() - cur_time)
            # cur_time = time.time()
            i = 1

        Vel = np.zeros(2*len(self.nodes))
        Pressure = np.zeros(3*len(self.elements))
        print('init Vel Pres', time.time() - cur_time)
        cur_time = time.time()

        #boundary conditions
        offset = (self.num_x_nodes - 1) * self.num_y_nodes * 2
        num_vel = self.num_x_nodes * self.num_y_nodes * 2
        bc_ind = np.concatenate((np.arange(0, self.num_y_nodes * 2, 2, dtype=int),                            # x vel left nodes
                           np.arange(offset, offset + self.num_y_nodes * 2, 2, dtype=int),          # x vel right nodes
                           np.arange(1, num_vel, self.num_y_nodes * 2, dtype=int),                  # y vel bottom nodes
                           np.arange(self.num_y_nodes * 2 - 1, num_vel, self.num_y_nodes * 2, dtype=int)))     # y vel top nodes
        bc_val = np.zeros(bc_ind.shape, dtype=np.float32)
        Vel[bc_ind] = bc_val
        print('boundary conditions', time.time() - cur_time)
        cur_time = time.time()

        free = np.arange(0, num_vel, dtype=int)
        free = np.delete(free, bc_ind)
        Rhs = Ra * Rhs_all - np.atleast_2d(np.matmul(A_all, Vel)).T
        A = A_all[np.ix_(free, free)]

        print('pre cholesky', time.time() - cur_time)
        cur_time = time.time()
        L = np.linalg.cholesky(A)
        L_T = L.T
        print('post cholesky', time.time() - cur_time)
        cur_time = time.time()
        PF = 1000

        div_max_goal = 1e-10
        div_max = 1
        iter_max = 5
        iter = 0
        while div_max > div_max_goal and iter < iter_max:
            intermed = np.linalg.solve(L, Rhs[free])
            Vel[free] = np.linalg.solve(L_T, intermed).T
            intermed = np.matmul(Q_all.T, np.atleast_2d(Vel).T)
            div = np.matmul(invM_all, intermed)
            Rhs -= PF * np.matmul(Q_all, div)
            Pressure += PF * div.T.ravel()
            div_max = max(abs(div))
            iter += 1

        print('loop exit', time.time() - cur_time)
        cur_time = time.time()
        self.Assign_to_space(Vel, Pressure)

    def Assign_to_space(self, Vel, Pres):
        for i in range(len(self.nodes)):
            vel = s.Vel(Vel[i * 2], Vel[i * 2 + 1])
            self.nodes[i].vel = vel

        for i in range(len(self.elements)):
            pressure = (Pres[i * 3], Pres[i * 3 + 1], Pres[i * 3 + 2])
            self.elements[i].pressure = pressure


    #endregion

    #region CONTINENT FORMATION



    #endregion

    #region PLOTTING FUNCTIONS
    def Get_Element_Ids(self):
        idArr = np.zeros((self.lenCol, self.lenRow))
        i = 0
        # print("Col:", self.lenCol, "Row:", self.lenRow)
        for y in range(self.lenCol):
            for x in range(self.lenRow):
                idArr[y][x] = self.elements[i].index
                i += 1
        return idArr

    def Plot_elements(self):
        plt.figure(figsize=(20, 20))
        plt.rc('font', size=25)
        cmap = mpl.cm.viridis
        norm = mpl.colors.Normalize(vmin=0, vmax=len(self.elements))
        scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        # plt.colorbar(scalarMap)
        plt.xlabel("Width")
        plt.ylabel("Depth")

        # plot elements
        idArr = self.Get_Element_Ids()
        extent = self.x_min, self.x_max, self.y_min, self.y_max
        plt.imshow(idArr, cmap=cmap, norm=norm, origin='lower', extent=extent)

        # plot nodes
        for node in self.nodes:
            plt.plot(node.x, node.y, marker='D', linestyle=':', markersize=12, color='white')

        plt.show()

    def Plot_elements_and_ips(self):
        plt.figure(figsize=(20, 20))
        plt.rc('font', size=25)
        cmap = mpl.cm.viridis
        norm = mpl.colors.Normalize(vmin=0, vmax=len(self.elements))
        scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        # plt.colorbar(scalarMap)
        plt.xlabel("Width")
        plt.ylabel("Depth")

        # plot elements
        idArr = self.Get_Element_Ids()
        extent = self.x_min, self.x_max, self.y_min, self.y_max
        plt.imshow(idArr, cmap=cmap, norm=norm, origin='lower', extent=extent)

        # plot nodes
        for node in self.nodes:
            plt.plot(node.x, node.y, marker='D', linestyle=':', markersize=12, color='white')

        for element in self.elements:
            for ip in element.ip:
                plt.plot(ip.x, ip.y, marker='.', markersize=10, color='white')

        plt.show()

    def Plot_elements_and_markers(self):
        plt.figure(figsize=(20, 20))
        plt.rc('font', size=25)
        cmap = mpl.cm.viridis
        norm = mpl.colors.Normalize(vmin=0, vmax=len(self.elements))
        scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        # plt.colorbar(scalarMap)
        plt.xlabel("Width")
        plt.ylabel("Depth")

        # plot elements
        idArr = self.Get_Element_Ids()
        extent = self.x_min, self.x_max, self.y_min, self.y_max
        plt.imshow(idArr, cmap=cmap, norm=norm, origin='lower', extent=extent)

        # plot nodes
        for node in self.nodes:
            plt.plot(node.x, node.y, marker='D', linestyle=':', markersize=12, color='white')

        # plot markers
        for marker in self.markers:
            elementColor = scalarMap.to_rgba(marker.container.index)
            plt.plot(marker.x, marker.y, marker='X',
                     markersize=30, color=elementColor, markeredgecolor='white', markeredgewidth=4)

        plt.show()

    def Plot_elements_and_marker_velocities(self):
        plt.figure(figsize=(20, 20))
        plt.rc('font', size=25)
        cmap = mpl.cm.viridis
        norm = mpl.colors.Normalize(vmin=0, vmax=len(self.elements))
        scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        plt.colorbar(scalarMap)
        plt.xlabel("Width")
        plt.ylabel("Depth")

        # plot elements
        idArr = self.Get_Element_Ids()
        extent = self.x_min, self.x_max, self.y_min, self.y_max
        plt.imshow(idArr, cmap=cmap, norm=norm, origin='lower', extent=extent)

        # plot markers
        for marker in self.markers:
            elementColor = scalarMap.to_rgba(marker.container.index)
            history = marker.pos_hist.keys()
            # for past in history:
            #     tuple_pos = marker.pos_hist[past]
            #     plt.scatter(tuple_pos[0], tuple_pos[1], color='r')
            plt.scatter(marker.x, marker.y, color='k')
            plt.quiver(marker.x, marker.y, marker.vel.v_x, marker.vel.v_y, color='w', angles='xy',
          scale_units='xy', scale=1)

        plt.show()

    def Plot_node_velocities(self):
        plt.figure(figsize=(20, 20))
        plt.rc('font', size=25)
        # cmap = mpl.cm.viridis
        # norm = mpl.colors.Normalize(vmin=0, vmax=len(self.elements))
        # scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        # plt.colorbar(scalarMap)
        plt.xlabel("Width")
        plt.ylabel("Depth")

        for i in range(len(self.elements)):
            element = self.elements[i]
            for j in range(len(element.all_nodes)):
                node = element.all_nodes[j]
                if self.phases[i] == 1:
                    col = 'r'
                else:
                    col = 'b'
                # plt.quiver(node.x, node.y, node.vel.v_x, node.vel.v_y, color='k', angles='xy', scale_units='xy', scale=40000)
                plt.scatter(node.x, node.y, color=col)

        plt.show()

    def vel_mag(self, x, y):
        vel = self.velEq(x, y)
        return np.sqrt(np.power(vel.v_x, 2) + np.power(vel.v_y, 2))

    def Plot_interpolated_markers(self):
        plt.figure(figsize=(20, 20))
        plt.rc('font', size=25)

        #plot global velocity
        lin_x = np.linspace(self.x_min, self.x_max, self.num_x_nodes * 2)
        lin_y = np.linspace(self.y_min, self.y_max, self.num_y_nodes * 2)
        X, Y = np.meshgrid(lin_x, lin_y)
        # print('X values:')
        # print(X)
        Z = np.zeros((self.num_y_nodes * 2, self.num_x_nodes * 2))
        # print('Z values:')
        # print(Z)
        v_max = 0
        for i in range(self.num_x_nodes * 2):
            for j in range(self.num_y_nodes * 2):
                x = X[j, i]
                y = Y[j, i]
                z = self.vel_mag(x, y)
                if z > v_max:
                    v_max = z
                Z[j, i] = z
                # print('x, y, z:', x, y, z)

        cmap = mpl.cm.viridis
        norm = mpl.colors.Normalize(vmin=0, vmax=v_max)
        scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        plt.colorbar(scalarMap)
        plt.xlabel("Width")
        plt.ylabel("Depth")
        # print('Z values:')
        # print(Z)
        plt.contourf(X, Y, Z, 20)

        vel_to_float = lambda z: np.sqrt(np.power(z.v_x, 2) + np.power(z.v_y, 2))

        #plot markers
        # print("Marker values:")
        for marker in self.markers:
            vel_mag = vel_to_float(marker.vel)
            # print(vel_mag)
            elementColor = scalarMap.to_rgba(vel_mag)
            # marker.Print_marker()
            plt.plot(marker.x, marker.y, marker='X', linestyle=':', markersize=13, color=elementColor)

        for node in self.nodes:
            vel_mag = vel_to_float(node.vel)
            elementColor = scalarMap.to_rgba(vel_mag)
            # node.Print_node()
            plt.plot(node.x, node.y, marker='D', linestyle=':', markersize=12, color=elementColor)

        plt.show()

    def Plot_density(self):
        plt.figure(figsize=(20, 20))
        plt.rc('font', size=25)

        # get min and max
        den_min = self.markers[0].density
        den_max = self.markers[0].density
        for marker in self.markers:
            if marker.density < den_min:
                den_min = marker.density
            elif marker.density > den_max:
                den_max = marker.density
        for node in self.nodes:
            if node.density < den_min:
                den_min = node.density
            elif node.density > den_max:
                den_max = node.density

        cmap = mpl.cm.viridis
        # print(den_min, den_max)
        norm = mpl.colors.Normalize(vmin=den_min, vmax=den_max)
        scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        plt.colorbar(scalarMap)
        plt.xlabel("Width")
        plt.ylabel("Depth")
        # print('Z values:')
        # print(Z)

        #plot markers
        # print("Marker values:")
        for marker in self.markers:
            # print(vel_mag)
            elementColor = scalarMap.to_rgba(marker.density)
            # marker.Print_marker()
            plt.plot(marker.x, marker.y, marker='X', linestyle=':', markersize=13, color=elementColor)

        for node in self.nodes:
            elementColor = scalarMap.to_rgba(node.density)
            # node.Print_node()
            plt.plot(node.x, node.y, marker='D', linestyle=':', markersize=12, color=elementColor)

        plt.show()

    def Get_shaped_arr(self, arr):
        shapedArr = np.zeros((self.lenCol, self.lenRow))
        i = 0
        # print("Col:", self.lenCol, "Row:", self.lenRow)
        for y in range(self.lenCol):
            for x in range(self.lenRow):
                shapedArr[y][x] = arr[i]
                i += 1
        return shapedArr

    def Plot_gauss_area(self):
        plt.figure(figsize=(20, 20))
        plt.rc('font', size=25)

        # determine range of colormap
        min = self.gauss_area[0]
        max = self.gauss_area[0]
        for area in self.gauss_area:
            if area < min:
                min = area
            elif area > max:
                max = area

        cmap = mpl.cm.viridis
        norm = mpl.colors.Normalize(vmin=min, vmax=max)
        scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        plt.colorbar(scalarMap)
        plt.xlabel("Width")
        plt.ylabel("Depth")

        # plot elements
        areaArr = self.Get_shaped_arr(self.gauss_area)
        extent = self.x_min, self.x_max, self.y_min, self.y_max
        plt.imshow(areaArr, cmap=cmap, norm=norm, origin='lower', extent=extent)

        # plot nodes
        for node in self.nodes:
            plt.plot(node.x, node.y, marker='D', linestyle=':', markersize=12, color='white')

        plt.show()

    def Plot_element_area(self):
        plt.figure(figsize=(20, 20))
        plt.rc('font', size=25)

        # 1D areas array and determine range of colormap
        areas = np.zeros((len(self.elements), 1))
        min = self.elements[0].area
        max = self.elements[0].area
        for i in range(len(self.elements)):
            area = self.elements[i].area
            areas[i] = area
            if area < min:
                min = area
            elif area > max:
                max = area

        cmap = mpl.cm.viridis
        norm = mpl.colors.Normalize(vmin=min, vmax=max)
        scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        plt.colorbar(scalarMap)
        plt.xlabel("Width")
        plt.ylabel("Depth")

        # plot elements
        areaArr = self.Get_shaped_arr(areas)
        extent = self.x_min, self.x_max, self.y_min, self.y_max
        plt.imshow(areaArr, cmap=cmap, norm=norm, origin='lower', extent=extent)

        # plot nodes
        for node in self.nodes:
            plt.plot(node.x, node.y, marker='D', linestyle=':', markersize=12, color='white')

        plt.show()

    def Plot_element_mass(self):
        plt.figure(figsize=(20, 20))
        plt.rc('font', size=25)

        # 1D areas array and determine range of colormap
        masses = np.zeros((len(self.elements), 1))
        min = self.elements[0].mass
        max = self.elements[0].mass
        for i in range(len(self.elements)):
            mass = self.elements[i].mass
            masses[i] = mass
            if mass < min:
                min = mass
            elif mass > max:
                max = mass

        cmap = mpl.cm.viridis
        norm = mpl.colors.Normalize(vmin=min, vmax=max)
        scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        plt.colorbar(scalarMap)
        plt.xlabel("Width")
        plt.ylabel("Depth")

        # plot elements
        areaArr = self.Get_shaped_arr(masses)
        extent = self.x_min, self.x_max, self.y_min, self.y_max
        plt.imshow(areaArr, cmap=cmap, norm=norm, origin='lower', extent=extent)

        # plot nodes
        for node in self.nodes:
            plt.plot(node.x, node.y, marker='D', linestyle=':', markersize=12, color='white')

        plt.show()

    def Plot_interpolation_error(self):
        plt.figure(figsize=(20, 20))

        marker_error = np.zeros((len(self.markers), 1))
        marker_to_node_dist = np.zeros((len(self.markers), 1))
        index = 0
        for marker in self.markers:
            global_vel = self.velEq(marker.x, marker.y)
            # marker.Print_marker()
            # global_vel.Print_Vel()
            # print(marker.vel.v_y)
            # print(global_vel.v_y)
            marker_error[index] = np.sqrt(np.power(marker.vel.v_x - global_vel.v_x, 2)
                                          + np.power(marker.vel.v_y - global_vel.v_y, 2))
            nearest_list = marker.container.corner_nodes
            nearest_dist = 1
            for node in nearest_list:
                dist = np.sqrt(np.power(node.x - marker.x, 2)
                                          + np.power(node.y - marker.y, 2))
                if dist < nearest_dist:
                    nearest_dist = dist
            marker_to_node_dist[index] = nearest_dist

            index += 1

        for i in range(len(marker_error)):
            plt.plot(marker_to_node_dist[i], marker_error[i], marker='o', markersize=12)

        plt.show()

    def Plot_gauss_error(self, max=0, min=0):
        plt.figure(figsize=(20, 20))
        plt.rc('font', size=25)

        # make difference array
        gaussArr = self.Get_shaped_arr(self.gauss_area)
        intArr = self.Get_shaped_arr(self.int_area)
        diffArr = gaussArr - intArr

        # determine range of colormap
        if (max == 0 and min == 0):
            min = diffArr[0][0]
            max = diffArr[0][0]
            for row in diffArr:
                for area in row:
                    print('area:', area)
                    if area < min:
                        min = area
                    elif area > max:
                        max = area

        cmap = mpl.cm.viridis
        norm = mpl.colors.Normalize(vmin=min, vmax=max)
        scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        plt.colorbar(scalarMap)
        plt.xlabel("Width")
        plt.ylabel("Depth")

        # plot elements
        extent = self.x_min, self.x_max, self.y_min, self.y_max
        plt.imshow(diffArr, cmap=cmap, norm=norm, origin='lower', extent=extent)

        plt.show()

        return (max, min)

    #endregion

    #region PRINT FUNCTIONS
    def print_nodes_and_elements(self):
        self.print_nodes()
        self.print_elements()

    def print_nodes(self):
        print("Nodes:")
        for node in self.nodes:
            node.Print_node()

    def print_node_list(self, nodes):
        print("Nodes:")
        for node in nodes:
            node.Print_node()

    def print_elements(self):
        print("Elements:")
        for element in self.elements:
            element.Print_element()

    def print_markers(self):
        print("Markers:")
        for marker in self.markers:
            marker.Print_marker()

    #endregion

    #region TESTING FUNCTIONS
    def Shift_nodes_by(self, x, y):
        for node in self.nodes:
            node.Move_node_by(x, y)
        #self.print_nodes_and_elements()

    def Nodes_know_elements(self):
        num_anchors = 0
        for node in self.nodes:
            # node.Print_node()
            if (node.defining):
                num_anchors += 1
                # print('Defining node:', end='')
                try:
                    # node.Print_node()
                    if (node == node.bottom_left_of.corner_nodes[0]):
                        #print('comparison')
                        #node.Print_node()
                        #node.bottom_left_of.nodes[0].Print_node()
                        pass
                    else:
                        return False
                except (e):
                    node.Print_node()
                    raise AttributeError("has no attribute")
        correct_num_anchors = len(self.defining_nodes)
        return num_anchors == correct_num_anchors

    def Confirm_marker_placement(self):
        for marker in self.markers:
            if (not hasattr(marker, 'container')):
                return False
        return True

    def Contains_marker(self, element, marker):
        elementNodes = element.corner_nodes
        min_x = elementNodes[0].x
        max_x = elementNodes[2].x
        min_y = elementNodes[0].y
        max_y = elementNodes[1].y
        return min_x <= marker.x <= max_x and min_y <= marker.y <= max_y

    def Validate_marker(self, marker):
        for i in range(len(self.elements)):
            if self.Contains_marker(self.elements[i], marker):
                if (not marker.container.index == self.elements[i].index and not
                        self.elements[i].markers.count(marker) == 0):
                    return False
        return True

    def Validate_markers(self):
        for marker in self.markers:
            if (not self.Validate_marker(marker)):
                return False
        return True

    def Markers_in_Elements(self, time):
        for marker in self.markers:
            element = marker.container
            if not self.Contains_marker(element, marker):
                print('MARKERS IN ELEMENTS')
                marker.Print_marker()
                print('is not contained in')
                element.Print_element()
                print('at time step', time)

    def Place_markers_at_corners(self):
        self.markers = np.array([m.Marker(0, 0),
                                 m.Marker(0, self.lenCol),
                                 m.Marker(self.lenRow, 0),
                                 m.Marker(self.lenRow, self.lenCol)])

    def Interpolation_consistency(self):
        self.Create_markers(1)
        constant_vel = lambda x, y: s.Vel(.5, .5)
        self.Set_velocity_eq(constant_vel)
        self.Set_marker_velocities()
        self.Plot_interpolated_markers()

    #endregion

#region EXECUTION FUNCTIONS
def Test_given_space(test):
    # tests that elements have a pointer to a node and not a copy
    init_node_y = test.nodes[0].y
    test.Shift_nodes_by(0, 1)
    assert test.nodes[0].y == init_node_y + 1
    assert test.elements[0].corner_nodes[0].y == test.nodes[0].y
    test.Shift_nodes_by(0, -1)

    # test that all nodes know what elements they are in
    assert test.Nodes_know_elements() == True

    # test that locate markers locates all markers
    assert test.Confirm_marker_placement() == True

    # test that locate marker does node location same as old locator
    assert test.Validate_markers() == True

    # test inconveniently placed markers
    test.Place_markers_at_corners()
    test.Locate_markers()
    assert test.Validate_markers() == True

def Test_Space():
    for i in range(2, 51):
        for j in range(2, 51):
            print('NEW SPACE:', i, j)
            test = Space(i, j, i, j)
            test.Create_markers(3)
            test.Locate_markers()
            Test_given_space(test)
    print("TESTS PASS")

def Advection_Sim(i, j, n, time_end, num_steps):
    space = Space(i, j, 30, 30)
    space.Create_markers(35)
    # space.Add_marker()
    space.Set_advection(1, n)

    time_step = time_end / num_steps

    time = [time_step for t in range(int(num_steps))]

    # print('Tf:', time_end)
    # print('Delta T:', time_step)

    for step in time:
        space.Set_marker_velocities()
        # store or plot values
        # space.Plot_elements_and_marker_velocities()

        space.Move_markers(step)

    print('done')
    space.Plot_elements_and_marker_velocities()

def Advection_error_generator(i):
    space = Space(i, i, 50, 50)
    space.Add_marker()
    marker_hist = space.markers[0].pos_hist
    space.Set_advection(1, 1)
    completed_loop = False
    cur_time = 0.0

    # first step to gain history
    space.Set_marker_velocities()
    space.Move_markers(1.0)
    cur_time += 1.0

    while not completed_loop:
        if (marker_hist.get(cur_time - 1.0)[1] > 0 and
            marker_hist.get(cur_time)[1] < 0):
            completed_loop = True
            start_x = marker_hist.get(0.0)[0]
            # interpolate for x value the moment y = 0
            final_x = space.Interpolate(marker_hist.get(cur_time)[0], marker_hist.get(cur_time)[1],
                                        marker_hist.get(cur_time - 1.0)[0], marker_hist.get(cur_time - 1.0)[1], 0)
            error = final_x - start_x

        space.Set_marker_velocities()
        space.Move_markers(1.0)
        cur_time += 1.0

    # print(i, 'error:', error)
    # print_list = [10, 15, 20, 30, 40, 60, 80]
    # if print_list.count(i) != 0:
    #     space.Plot_elements_and_marker_velocities()
    return (error, space.delta_x)

def Plot_advection_error():
    errors = np.zeros((90, 1))
    element_size = np.zeros((90, 1))
    index = 0
    for i in range(10, 100):
        difference, delta_x = Advection_error_generator(i)
        # absolute error
        # errors[index] = difference
        # error proportional to cell size
        errors[index] = np.power(difference, 2) / np.power(delta_x, 2)
        element_size[index] = delta_x
        index += 1

    plt.figure(figsize=(20, 20))
    plt.rc('font', size=25)
    plt.xlabel("Delta x")
    plt.ylabel("Error")

    for i in range(len(errors)):
        plt.plot(element_size[i], errors[i], marker='o', markersize=12)

    plt.show()

def Test_Interpolation(i, j, n):
    space = Space(i, j, i, j)
    space.Create_markers(1)
    space.Set_advection(1, n)
    # constant_vel = lambda x, y: s.Vel(.5, .5)
    # space.Set_velocity_eq(constant_vel)
    space.Set_marker_velocities()

    space.Plot_interpolated_markers()
    # space.Plot_interpolation_error()

def Check_gauss_area():
    space = Space(5, 5, 5, 5)
    f_0 = lambda x, y: 1
    f_1 = lambda x, y: x + y
    for i in range(space.num_y_nodes):
        space.nodes[i].Move_node_by(.5, 0)
    space.gauss_quad_no_interp(f_0)
    print(space.gauss_area)
    print(space.print_nodes())
    space.Plot_gauss_area()

def Gauss_error():
    space = Space(20, 20, 20, 20)
    density_2 = lambda x, y: np.sin(x) + np.cos(y)
    space.Set_density_at_nodes(density_2)
    space.Gauss_quad()
    int_2 = lambda a, b, c, d: (np.cos(a) - np.cos(b)) * (d - c) + \
                               (b - a) * (np.sin(d) - np.sin(c))
    space.Int_area(int_2)
    space.Plot_gauss_area()
    max, min = space.Plot_gauss_error()

    density_1 = lambda x, y: x + y
    space.Set_density_at_nodes(density_1)
    space.Gauss_quad()
    int_1 = lambda a, b, c, d: (np.power(b, 2) - np.power(a, 2)) / 2 * (d - c) + \
                               (b - a) / 2 * (np.power(d, 2) - np.power(c, 2))
    space.Int_area(int_1)
    space.Plot_gauss_area()
    space.Plot_gauss_error(max, min)

def Gauss_error_no_interp():
    space = Space(20, 20, 20, 20)
    density_2 = lambda x, y: np.sin(x) + np.cos(y)
    space.gauss_quad_no_interp(density_2)
    int_2 = lambda a, b, c, d: (np.cos(a) - np.cos(b)) * (d - c) + \
                               (b - a) * (np.sin(d) - np.sin(c))
    space.Int_area(int_2)
    space.Plot_gauss_area()
    max, min = space.Plot_gauss_error()

    density_1 = lambda x, y: x + y
    space.gauss_quad_no_interp(density_1)
    int_1 = lambda a, b, c, d: (np.power(b, 2) - np.power(a, 2)) / 2 * (d - c) + \
                               (b - a) / 2 * (np.power(d, 2) - np.power(c, 2))
    space.Int_area(int_1)
    space.Plot_gauss_area()
    space.Plot_gauss_error(max, min)

def Gauss_error_interp_comp():
    space = Space(20, 20, 20, 20)
    density_2 = lambda x, y: np.sin(x) + np.cos(y)
    space.Set_density_at_nodes(density_2)
    space.Gauss_quad()
    int_2 = lambda a, b, c, d: (np.cos(a) - np.cos(b)) * (d - c) + \
                               (b - a) * (np.sin(d) - np.sin(c))
    space.Int_area(int_2)
    space.Plot_gauss_area()
    max, min = space.Plot_gauss_error()

    density_1 = lambda x, y: np.sin(x) + np.cos(y)
    space.gauss_quad_no_interp(density_1)
    int_1 = lambda a, b, c, d: (np.cos(a) - np.cos(b)) * (d - c) + \
                               (b - a) * (np.sin(d) - np.sin(c))
    space.Int_area(int_1)
    space.Plot_gauss_area()
    max, min = space.Plot_gauss_error(max, min)

def Density():
    space = Space(10, 10, 10, 10)
    space.Create_markers(10)
    d_1 = lambda x, y: np.sin(x) + np.cos(y) + 5
    space.Set_density(d_1)
    space.Update_Node_density()
    space.Plot_density()

def Calculate_mass(lenRow, main, spot):
    space = Space(lenRow, lenRow, 100, 100)
    space.Create_markers(35)
    heavy_x = random.randint(-50, 40)
    heavy_y = random.randint(-50, 40)
    # print('heavy range:', heavy_x, heavy_y)
    for marker in space.markers:
        if heavy_x <= marker.x and marker.x <= heavy_x + 10 \
                and heavy_y <= marker.y and marker.y <= heavy_y + 10:
            marker.density = spot
        else:
            marker.density = main
    space.Update_Node_density()
    mass = space.Update_Element_mass()

    # space.Plot_density()
    # space.Plot_element_area()
    # space.Plot_element_mass()
    return mass

# def generate_main(lenRow, totalBaseMass):
#     main = list()
#     for entry in lenRow:
#         main.append(totalBaseMass / np.power(entry, 2))
#     return main
#
# def generate_spot(main, totalBaseMass):
#     spot = list()
#     five_percent = totalBaseMass * 0.05
#     for baseVal in main:
#         spot.append(five_percent + baseVal)
#     return spot

def Mass_distribution_comparison():
    totalBaseMass = 36000
    lenRow = [5, 10, 15, 20, 25, 30]
    main = 3.6
    spot = 21.6

    numIt = 20
    mass = np.zeros((len(lenRow), numIt))

    for i in range(len(lenRow)):
        print(i)
        for j in range(numIt):
            mass[i][j] = Calculate_mass(lenRow[i], main, spot)

    plt.figure(figsize=(20, 20))
    plt.rc('font', size=25)
    plt.xlabel("Width of Space")
    plt.ylabel("Total Mass")

    for i in range(len(mass)):
        plt.plot(lenRow[i], totalBaseMass + totalBaseMass * 0.05, marker='s', markersize=17, color='black')
        for j in range(len(mass[0])):
            plt.plot(lenRow[i], mass[i][j], marker='o', markersize=10, color='magenta')

    plt.show()

def Mass_Conservation(i, j, n, time_end, num_steps):
    space = Space(i, j, 100, 100)
    space.Create_markers(35)
    space.Set_advection(1, n)

    time_step = time_end / num_steps
    time = [time_step for t in range(int(num_steps))]

    # set marker densities
    total_mass = 22000
    heaviest_x = random.randint(-50, 40)
    heaviest_y = random.randint(-50, 40)
    heavy_x = random.randint(-50, 0)
    heavy_y = random.randint(-50, 0)
    light_x = random.randint(-50, 20)
    light_y = random.randint(-50, 20)
    for marker in space.markers:
        if heaviest_x <= marker.x <= heaviest_x + 10 \
                and heaviest_y <= marker.y <= heaviest_y + 10:
            marker.density += 4
        if heavy_x <= marker.x <= heavy_x + 50 \
                and heavy_y <= marker.y <= heavy_y + 50:
            marker.density += 1
        if light_x <= marker.x <= light_x + 30 \
                and light_y <= marker.y <= light_y + 30:
            marker.density -= 1

        marker.density += 2

    masses = list()

    # simulation
    accrued_time = 0
    for step in time:
        # space.Markers_in_Elements(accrued_time)

        space.Set_marker_velocities()

        space.Update_Node_density()
        mass = space.Update_Element_mass()
        masses.append(mass)
        # space.Plot_density()

        space.Move_markers(step)
        accrued_time += step

    space.Plot_density()

    plt.figure(figsize=(20, 20))
    plt.rc('font', size=25)
    plt.xlabel("Time")
    plt.ylabel("Total Mass")
    five_percent = total_mass * 0.05
    plt.xlim(0, time_end)
    plt.ylim(total_mass - five_percent, total_mass + five_percent)

    plt.plot([0, time_end], [total_mass, total_mass], color='black')
    t = 0
    for m in masses:
        plt.plot(t, m, marker='o', markersize=10, color='magenta')
        t += time_step

    plt.show()

def Test_Stokes_Time():
    smallest, largest = 10, 11
    time_arr = np.empty(largest - smallest)
    bins = np.arange(smallest, largest + 1)

    time_index = 0
    for i in range(smallest, largest):
        print(i)
        space = Space(i, i, 100, 100)
        space.Set_local_ip()
        space.Set_Phases()
        # space.Set_element_ip()    # do I need this?
        start = time.time()
        print('0')
        space.Stokes_solver(start)
        stop = time.time()
        print('time', stop - start)
        time_arr[time_index] = stop - start
        time_index += 1
        space.Plot_node_velocities()
        space.print_nodes()

    plt.stairs(time_arr, bins, fill=True)
    plt.rc('font', size=25)
    plt.xlabel("Number of Elements in x")
    plt.xticks(np.arange(bins[0], bins[-1] + 1))
    plt.ylabel("Time to run Stokes solver")
    # plt.show()


#endregion

# Density()
# Calculate_mass(10, 1, 2)
# Mass_distribution_comparison()
# Mass_Conservation(10, 10, 1, 80, 160)

# Plot_advection_error()

# Gauss_error_interp_comp()

# Check_gauss_area()

Test_Stokes_Time()

# space = Space(2, 2, 1, 1)
# space.Create_markers(3)
# space.Plot_elements_and_markers()

# space.Set_element_ip()
# space.Stokes_solver()
#space.Plot_elements_and_ips()

# Test_Space()
# Test_Interpolation(10, 10, 1)
# test = Space(5, 5)
# test.Interpolation_consistency()
# Advection_Sim(10, 10, 1, 250, 70)
# Advection_Sim(40, 40, 1, 250, 70)


# 1 million elements takes 10 seconds
#test = Space(5, 5)
#test.Create_markers()
#test.Plot_elements_and_markers()



# test.print_nodes_and_elements()
# test.Place_markers_at_corners()
# test.print_markers()
# test.Create_markers()
# test.Locate_markers()
# test.print_markers()
# test.Plot_elements_and_markers()

# v = lambda x, y : (2 * x, round(np.square(y), 2))
# v1 = lambda x, y : (round(np.sin(2 * np.pi * y / 3), 2), 0)
# test.Set_velocity_eq(v)
# test.Set_marker_velocities()
#
# test.Plot_elements_and_marker_velocities()




### 4 NODED ELEMENTS ###
# def __init__(self, lenRow, lenCol, width, height):
#     self.width = width
#     self.height = height
#     self.lenRow = lenRow
#     self.lenCol = lenCol
#     self.num_x_nodes = lenRow + 1  # num nodes in x direction
#     self.num_y_nodes = lenCol + 1  # num nodes in y direction
#     self.delta_x = width / lenRow
#     self.delta_y = height / lenCol
#
#     s.IDNum.reset()
#
#     self.nodes = np.full(self.num_x_nodes * self.num_y_nodes, None)
#     i = 0
#     x_scalar = self.width / (self.num_x_nodes - 1)
#     y_scalar = self.height / (self.num_y_nodes - 1)
#     x_start = self.width / -2
#     y_start = self.height / -2
#     for x in range(self.num_x_nodes):
#         for y in range(self.num_y_nodes):
#             # self.nodes[i] = n.Node(x, y, i)
#             self.nodes[i] = n.Node(x_start + x * x_scalar, y_start + y * y_scalar, i)
#             i += 1
#
#     self.x_min = self.nodes[0].x
#     self.x_max = self.nodes[self.num_x_nodes * self.num_y_nodes - 1].x
#     self.y_min = self.nodes[0].y
#     self.y_max = self.nodes[self.num_y_nodes - 1].y
#
#     self.elements = np.full(lenRow * lenCol, None)
#     i = 0
#     for y in range(lenCol):
#         for x in range(lenRow):
#             # nodes for a grid
#             b_l = self.nodes[x * self.num_y_nodes + y]
#             t_l = self.nodes[(x * self.num_y_nodes + y) + 1]
#             b_r = self.nodes[(x + 1) * self.num_y_nodes + y]
#             t_r = self.nodes[((x + 1) * self.num_y_nodes + y) + 1]
#             self.elements[i] = e.Element(b_l, t_l, b_r, t_r)
#             i += 1
#     # self.print_nodes_and_elements()



### FLEXIBLE INIT ###
# def __init__(self, numRow, numCol):
#     self.numRow = numRow
#     self.numCol = numCol
#     self.num_x_nodes = numRow + 1
#     self.num_y_nodes = numCol + 1
#
#     self.nodes = np.full(self.num_x_nodes * self.num_y_nodes, None)
#     i = 0
#     for x in range(self.num_x_nodes):
#         for y in range(self.num_y_nodes):
#             self.nodes[i] = node.Node(x, y)
#             i += 1
#
#     self.elements = np.full(numRow * numCol, None)
#     i = 0
#     for y in range(numCol):
#         for x in range(numRow):
#             # nodes for a grid
#             node1 = self.nodes[x * self.num_y_nodes + y]
#             node2 = self.nodes[(x * self.num_y_nodes + y) + 1]
#             node3 = self.nodes[(x + 1) * self.num_y_nodes + y]
#             node4 = self.nodes[((x + 1) * self.num_y_nodes + y) + 1]
#             self.elements[i] = element.Element(node1, node2, node3, node4)
#             i += 1
#     self.print_nodes_and_elements()


### OLD CODE ###
#     def __init__(self, rowpts, colpts):
#         self.rowpts = rowpts
#         self.colpts = colpts
#
#         self.points = np.zeros((rowpts, colpts))
#
#     def Plot_space(self):
#         # x = np.linspace(0, self.rowpts, self.rowpts - 1)
#         # y = np.linspace(0, self.colpts, self.colpts - 1)
#         # X, Y = np.meshgrid(x, y)
#
#         #self.points = linear(X, Y)
#         #self.Set_points()
#         self.Index_space()
#
#         plt.figure(figsize=(20, 20))
#         #plt.contourf(X, Y, self.points)
#         #plt.plot(X, Y, marker=".", color='k', linestyle='none')
#         extent = np.min(0), np.max(self.rowpts), np.min(0), np.max(self.colpts)
#         plt.imshow(self.points, extent=extent)
#         plt.colorbar()
#         plt.xlabel("Width")
#         #plt.xticks(x)
#         plt.ylabel("Depth")
#         #plt.yticks(y)
#         plt.show()
#
#     def Plot_points(self):
#         self.points = np.array([[1, 1, 1, 1, 1],
#                                [0, 0, 1, 1, 0],
#                                [0, 0, 0, 0, 0],
#                                [0, 0, 0, 0, 0],
#                                [0, 0, 0, 0, 0]])
#
#     def Index_space(self):
#         num = (self.rowpts) * (self.colpts)
#         self.index = np.arange(1, num)
#         i = 0
#         for x in range(self.rowpts):
#             for y in range(self.colpts):
#                 self.points[x][y] = i
#                 i += 1
#         print(self.points)
#
# def linear(x, y):
#     return y > np.square(x)
#
# def sin(x, y):
#     return np.sin(x) + np.sin(y)
#
#
# test = Space(5, 5)
# #test.Index_space()
# test.Plot_space()



