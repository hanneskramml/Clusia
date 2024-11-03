import serial
from  copy import copy
from struct import *


ser = serial.Serial(
  port='/dev/cu.usbmodem14101',
    baudrate=9600
)

def get_next_data_block(next_f):
    if not hasattr(get_next_data_block, "data_block"):
        get_next_data_block.data_block = []
    while (1):
        try:
            current_item = next_f()
            if current_item == '^':
                next_item = next_f()
                if next_item == '^':
                    get_next_data_block.data_block.append(next_item)
                else:
                    out = copy(get_next_data_block.data_block)
                    get_next_data_block.data_block = []
                    get_next_data_block.data_block.append(next_item)
                    return out
            else:
                get_next_data_block.data_block.append(current_item)
        except :
            break


for i in range(1000): # just so that the program ends - could be in a while loop
    data_ =  get_next_data_block(ser.read)
    try:
        print(unpack('=ffffffffffff', ''.join(data_)))
    except:
        continue
