import time
import datetime
import schedule
from arduino import Arduino
from excel import Excel
from owncload import Cloud
from config import RPi, ArduinoDARK, ArduinoHL
from data import Measurement


class ClusiaRPi:
    def __init__(self):
        self.arduinoDARK = Arduino(name=ArduinoDARK.NAME, port=ArduinoDARK.PORT)
        self.arduinoHL = Arduino(name=ArduinoHL.NAME, port=ArduinoHL.PORT)
        self.cloud = Cloud()

    def test(self):
        self.arduinoDARK.connect()
        self.arduinoHL.connect()

        while True:
            self.arduinoDARK.measure()
            self.arduinoHL.measure()
            time.sleep(10)

    def run(self):
        self.arduinoDARK.connect()
        self.arduinoHL.connect()

        schedule.every(RPi.MEASUREMENT_INTERVAL).minutes.at(":00").do(self.measurement, device=self.arduinoDARK, config=ArduinoDARK)
        schedule.every(RPi.MEASUREMENT_INTERVAL).minutes.at(":30").do(self.measurement, device=self.arduinoHL, config=ArduinoHL)

        while True:
            schedule.run_pending()
            time.sleep(1)

    # TODO: Timeout
    def measurement(self, device, config):
        date = datetime.date.today()
        timeofday = datetime.datetime.now().time().strftime('%H:%M')
        print("Requesting measurement from {} ({} {})...".format(config.NAME, date, timeofday))

        if device.isConnected:
            device.measure()
            time.sleep(1.5)
        else:
            print("ERROR! Device not connected: {} ({})".format(device.name, device.port))
            return

        #while self.arduinoHL.dataQueue.qsize() <= config.SENSORS.__len__():

        print("Validating data from {}...".format(device.name))
        data = Measurement(config, date, timeofday)

        while not device.dataQueue.empty():
            sensor_code, sensor_value = device.dataQueue.get()

            if config.SENSORS.get(sensor_code) is not None:
                data.set(code=sensor_code, value=sensor_value)
            else:
                print("ERROR! Sensor not configured: {}".format(sensor_code))
                continue

        print("Writing data...")
        data.save()

        # TODO: Move to separate thread
        print("Uploading files...")
        self.cloud.upload(RPi.ENV_FILE)
        self.cloud.upload(RPi.LEAF_FILE)
        self.cloud.upload(RPi.SOIL_FILE)


if __name__ == '__main__':
    clusiaRPi = ClusiaRPi()
    clusiaRPi.run()

