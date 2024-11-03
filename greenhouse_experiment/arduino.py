from threading import Thread
import serial
import struct
import queue

SERIAL_PORT = '/dev/cu.usbmodem14101'
SERIAL_BAUDRATE = 115200
DATA_FMT = '<ccccf'
DATA_START = '<'.encode()
DATA_END = '>>'.encode()


class Arduino:
    def __init__(self, name='Arduino', port=SERIAL_PORT, baudrate=SERIAL_BAUDRATE, data_fmt=DATA_FMT):
        self.name = name
        self.port = port
        self.baudrate = baudrate
        self.dataFMT = data_fmt
        self.dataQueue = queue.Queue()

        self.serial = None
        self.thread = None
        self.isConnected = False

    def connect(self):
        print('Connecting to {} ({}) at baudrate {}...'.format(self.name, str(self.port), str(self.baudrate)))
        try:
            self.serial = serial.Serial(self.port, self.baudrate, timeout=10)
            self.isConnected = True
            print('Serial connection established! (RPi=>{})'.format(self.name))
        except Exception as e:
            self.isConnected = False
            print(e)
            return

        #if self.thread is None:
        self.thread = Thread(target=self.__thread)
        self.thread.start()

    def __thread(self):
        self.serial.reset_input_buffer()
        line = ""

        while True:
            try:
                x = self.serial.read(1)
            except Exception as e:
                self.isConnected = False
                print("FATAL ERROR reading from {}: {}".format(self.name, e))
                return

            if x == DATA_START:
                try:
                    rawData = self.serial.read(struct.calcsize(self.dataFMT))
                except Exception as e:
                    self.isConnected = False
                    print("FATAL ERROR reading from {}: {}".format(self.name, e))
                    return

                if self.serial.read(DATA_END.__len__()) == DATA_END:
                    try:
                        data_frame = struct.unpack(self.dataFMT, rawData)
                        data = (''.join(data.decode() for data in data_frame[0:4]), data_frame[4])

                        self.dataQueue.put(data)
                        print("Data received: {}".format(data))

                    except Exception as e:
                        print("ERROR receiving data frame: {}".format(rawData))

                else:
                    print("ERROR receiving data frame: {}".format(rawData))

            else:
                try:
                    line += x.decode()
                except:
                    print("ERROR unable to decode serial character: {}".format(x))

                if (line.__len__() > 0 and line[-1] == '\n') or line.__len__() > 100:
                    log = line.strip('\r\n')
                    line = ""

                    if log.__len__() > 0:
                        print("{}: {}".format(self.name, log))
                    else:
                        print(log)

    def measure(self):
        self.serial.write('GO!'.encode())
