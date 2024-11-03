import serial, struct

HEADER_BEGIN = bytes([0xAA, 0x55])
HEADER_END = bytes([0x55, 0xAA])
DATAFMT = '<ff'

with serial.Serial('/dev/cu.usbmodem14101', baudrate=115200) as ser:
    ser.reset_input_buffer()

    while True:
        line = bytearray()
        while True:
            c = ser.read(1)
            if c:
                line += c
                #print(line)
                if line[-len(HEADER_BEGIN):] == HEADER_BEGIN:
                    # height_m, boom_pos_m = struct.unpack(
                    #     DATAFMT, ser.read(struct.calcsize(DATAFMT)),
                    # )
                    # print(height_m, boom_pos_m)
                    data = ser.read_until(HEADER_END)
                    print("Data received: {}".format(data[:-2]))
                    break
                #if line[-len(HEADER2):] == HEADER2:
                #    print(ser.readline())
                #    break
            else:
                break

        print(bytes(line[:-2]))

        # ser.read_until(HEADER)
        # ser.readline()
        # mode = ser.read(1)
        # print(mode)
        # print(HEADER_DATA)
        # #print(HEADER)
        # print(mode == HEADER_DATA)
        #
        #
        # if mode == HEADER_DATA:
        #     height_m, boom_pos_m = struct.unpack(DATAFMT, ser.read(struct.calcsize(DATAFMT)),)
        #     print(height_m, boom_pos_m)
        # else:
        #     print(ser.read_until(HEADER))

        # mode = ser.read(1)
        # print(mode)
        #
        # if mode == HEADER2:
        #     height_m, boom_pos_m = struct.unpack(
        #         DATAFMT, ser.read(struct.calcsize(DATAFMT)),
        #     )
        #     print(height_m, boom_pos_m)
        # else:
        #     print(ser.readline())


