import structs as s

class Node:

    def __init__(self, x, y, index=None):
        self.x = x
        self.y = y
        self.vel = s.Vel(0, 0)
        self.density = 0
        self.defining = False
        self.index = index
        # should elements have unique ids the nodes can also have?

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

    def Set_element(self, element):
        self.bottom_left_of = element
        self.defining = True
        #self.Print_node()

    def Set_Vel(self, v_tuple):
        self.vel = v_tuple

    def Set_Density(self, den):
        self.density = den

    def Print_node(self):
        try:
            print("N:", self.x, ",", self.y, ' Def: ', self.defining, ' V: <',
                  round(self.vel.v_x, 2), ', ', round(self.vel.v_y, 2), '>', sep="")
        except AttributeError:
            print("N:", self.x, ",", self.y, ' Def: ', self.defining, sep="")

    def Move_node_by(self, x, y):
        self.x += x
        self.y += y