import os
import csv
from config import RPi


class Measurement:
    def __init__(self, config, date, time):
        self.config = config

        self.env_header = ['date', 'timeofday', 'condition', 'CO2', 'Light', 'Temp', 'Humidity']
        self.leaf_header = ['date', 'timeofday', 'condition', 'rep', 'LeafTemp', 'AmbientTemp', 'Delta']
        self.soil_header = ['date', 'timeofday', 'condition', 'rep', 'SoilTemp', 'SoilMoisture']

        self.env_format = {'date': date, 'timeofday': time, 'condition': config.CONDITION, 'CO2': "NA", 'Light': "NA", 'Temp': "NA", 'Humidity': "NA"}
        self.leaf_format = {'date': date, 'timeofday': time, 'condition': config.CONDITION, 'rep': "NA", 'LeafTemp': "NA", 'AmbientTemp': "NA", 'Delta': "NA"}
        self.soil_format = {'date': date, 'timeofday': time, 'condition': config.CONDITION, 'rep': "NA", 'SoilTemp': "NA", 'SoilMoisture': "NA"}

        self.env_data = {}
        self.leaf_data = {}
        self.soil_data = {}

        if not os.path.exists(RPi.ENV_FILE):
            env_file = csv.DictWriter(open(RPi.ENV_FILE, 'wt'), fieldnames=self.env_header)
            env_file.writeheader()

        if not os.path.exists(RPi.LEAF_FILE):
            leaf_file = csv.DictWriter(open(RPi.LEAF_FILE, 'wt'), fieldnames=self.leaf_header)
            leaf_file.writeheader()

        if not os.path.exists(RPi.SOIL_FILE):
            soil_file = csv.DictWriter(open(RPi.SOIL_FILE, 'wt'), fieldnames=self.soil_header)
            soil_file.writeheader()

    def set(self, code, value):
        category, name, replicate = self.config.SENSORS.get(code)

        if category == "ENV":
            if value < 0:
                print("ENV|{}: Invalid sensor value: {}".format(name, value))
                return

            if not self.env_data:
                self.env_data = self.env_format.copy()

            if name == "CO2":
                self.env_data[name] = int(value)
            else:
                self.env_data[name] = round(value, 2)

        elif category == "LEAF":
            if value > 100:
                print("LEAF|{}: Invalid sensor value: {}".format(name, value))
                return

            if not self.leaf_data.get(replicate):
                self.leaf_data[replicate] = self.leaf_format.copy()
                self.leaf_data[replicate]['rep'] = replicate

            self.leaf_data[replicate][name] = round(value, 2)

            if self.leaf_data[replicate]['LeafTemp'] != "NA" and self.leaf_data[replicate]['AmbientTemp'] != "NA":
                self.leaf_data[replicate]['Delta'] = round(self.leaf_data[replicate]['LeafTemp'] - self.leaf_data[replicate]['AmbientTemp'], 2)

        elif category == "SOIL":
            if value < 0:
                print("SOIL|{}: Invalid sensor value: {}".format(name, value))
                return

            if not self.soil_data.get(replicate):
                self.soil_data[replicate] = self.soil_format.copy()
                self.soil_data[replicate]['rep'] = replicate

            self.soil_data[replicate][name] = round(value, 2)

        else:
            print("ERROR! Sensor category not properly set: {}".format(category))

    def save(self):
        if self.env_data:
            env_file = csv.DictWriter(open(RPi.ENV_FILE, 'at'), fieldnames=self.env_header)
            env_file.writerow(self.env_data)

        if self.leaf_data:
            leaf_file = csv.DictWriter(open(RPi.LEAF_FILE, 'at'), fieldnames=self.leaf_header)
            leaf_file.writerows(self.leaf_data.values())

        if self.soil_data:
            soil_file = csv.DictWriter(open(RPi.SOIL_FILE, 'at'), fieldnames=self.soil_header)
            soil_file.writerows(self.soil_data.values())
