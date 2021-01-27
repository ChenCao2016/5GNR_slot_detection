import numpy as np
import array
import scipy.signal as signal

class IQsampleIn:

    class input:
        filename = None
        printInfo = None

    class output:
        data = None
        length = None

    def __init__(self):
        pass
    
    def run(self):
        f = open(self.input.filename, 'rb')
        packData = f.read()
        f.close()
        IQsample = np.array(array.array('f',packData),dtype=np.complex)
        self.output.data = np.add(IQsample[0::2], 1j*IQsample[1::2])
        self.output.length = np.size(self.output.data)
        self.info()
        return

    def info(self):
        if self.input.printInfo:
            print("\n#" + self.__class__.__name__ + "#")
            print(self.input.filename)
            print("Number of output samples: " + str(self.output.length))
            print("---------------------------------------------------")

class reSampler:
    class input:
        data = None
        targetRatio = None
        printInfo = None

    class output:
        data = None
        length = None

    def __init__(self):
        pass
    
    def run(self):
        self.output.data = signal.resample(self.input.data, int(np.size(self.input.data)*self.input.targetRatio))
        self.output.length = np.size(self.output.data)
        self.info()
        return

    def info(self):
        if self.input.printInfo:
            print("\n#" + self.__class__.__name__ + "#")
            print("Number of output samples: " + str(self.output.length))
            print("---------------------------------------------------")