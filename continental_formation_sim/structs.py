class Vel:

    def __init__(self, v_x, v_y):
        self.v_x = v_x
        self.v_y = v_y

    def Print_Vel(self):
        print('<', self.v_x, ", ", self.v_y, ">", sep='')


class IDNum:

    counter = -1

    @staticmethod
    def getID():
        IDNum.counter += 1
        return IDNum.counter

    @staticmethod
    def reset():
        IDNum.counter = -1