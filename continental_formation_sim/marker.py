import structs as s

class Marker:

    def __init__(self, x, y, bound_min_x=None, bound_max_x=None, bound_min_y=None, bound_max_y=None):
        self.x = x
        self.y = y
        self.vel = s.Vel(0, 0)
        self.density = 0
        self.container = None

        self.bound_min_x = bound_min_x
        self.bound_max_x = bound_max_x
        self.bound_min_y = bound_min_y
        self.bound_max_y = bound_max_y

        self.pos_hist = { 0.0: (x, y) }
        self.most_recent_time = 0.0

    def __lt__(self, other):
        if (self.x < other.x):
            return True
        elif (self.x == other.x and self.y < other.y):
            return True
        return False

    def __gt__(self, other):
        if (self.x > other.x):
            return True
        elif (self.x == other.x and self.y > other.y):
            return True
        return False

    def __le__(self, other):
        if (self.x < other.x):
            return True
        elif (self.x == other.x and self.y <= other.y):
            return True
        return False

    def __ge__(self, other):
        if (self.x > other.x):
            return True
        elif (self.x == other.x and self.y >= other.y):
            return True
        return False

    def Place_in_element(self, element):
        self.container = element

    def Set_Vel(self, v_tuple):
        self.vel = v_tuple

    def Set_Density(self, den):
        self.density = den

    def Move_marker(self, time_step):
        self.x += self.vel.v_x * time_step
        self.y += self.vel.v_y * time_step
        self.most_recent_time += time_step
        self.pos_hist[self.most_recent_time] = (self.x, self.y)

    def Keep_in_bounds(self):
        if (self.x < self.bound_min_x):
            self.x = self.bound_min_x
        elif (self.x > self.bound_max_x):
            self.x = self.bound_max_x

        if (self.y < self.bound_min_y):
            self.y = self.bound_min_y
        elif (self.y > self.bound_max_y):
            self.y = self.bound_max_y

    def Print_marker(self):
        try:
            print('M:(', self.x, ', ', self.y, ') E: ', self.container.index, ' V: <',
                  round(self.vel.v_x, 2), ', ', round(self.vel.v_y, 2), '>', sep='')
        except AttributeError:
            try:
                print('M:(', self.x, ', ', self.y, ') E: ', self.container.index, sep='')
            except AttributeError:
                print('M/E:(', self.x, ', ', self.y, ')', sep='')
        #self.container.Print_element()