#include <Arduino.h>              //Analog communication (Soil moisture)
#include <SoftwareSerial.h>       //Serial communication (CO2)
#include <Wire.h>                 //SDA-SCL communication (Light & IR)
#include <OneWire.h>              //Digital bus communication (Soil temperatures)

#include <MHZ19.h>                //CO2 sensor (MH-Z19C)
#include <BH1750.h>               //Light sensor (BH1759)
#include <DHT.h>                  //Air temp & humidity (DHT22)
#include <Adafruit_MLX90614.h>    //Infrared sensor (Adafruit MLX90614)
#include <DallasTemperature.h>    //Soil temperature sensor (DS18B20)

#define CO2_BAUDRATE 9600
#define DHT_TYPE22 DHT22
#define DATACODE_LEN 4



struct DATAFRAME
{
  char code[DATACODE_LEN];
  float value {};
};

class Sensor {
  String sid;
  
  public:
    Sensor(String sid) { this->sid = sid; }
    String getSid() { return this->sid; }
      
    void transmit(float value) { this->transmit(this->sid, value); }
    void transmit(String code, float value) {
      DATAFRAME data;
      code.toCharArray(data.code, DATACODE_LEN + 1);
      data.value = value;
      
      unsigned long uBufSize = sizeof(data);
      char pBuffer[uBufSize];
      memcpy(pBuffer, &data, uBufSize);
    
      Serial.write('<');
      Serial.write((uint8_t*)&pBuffer, sizeof(pBuffer));
      Serial.write(">>");
    }

    void printErr(String err) { 
      Serial.println("ERROR (" + this->sid + "): " + err); 
    }

    virtual bool measure() {}
  
};

class CO2 : public Sensor {
  MHZ19 device;
  SoftwareSerial *serial;

  public:
    CO2 (String sid, int rx_pin, int tx_pin) : Sensor(sid) {
      pinMode(rx_pin, INPUT);
      pinMode(tx_pin, OUTPUT);
      
      serial = new SoftwareSerial(rx_pin, tx_pin);
      serial->begin(CO2_BAUDRATE);
      
      device.begin(*serial);
      device.autoCalibration();
      device.setFilter(true, false);
    }
  
    bool measure() {
      this->transmit(device.getCO2());
    }
  
};

class Light : public Sensor {
  BH1750 sensor;

  public:
    Light (String sid) : Sensor(sid) {
      Wire.begin();
      if (!sensor.begin())
        this->printErr("Unable to initialise BH1750 light sensor !!!");
        
    }

    bool measure() {
      this->transmit(sensor.readLightLevel(true));
    }
    
};

class Dht22 : public Sensor {
  DHT *sensor;

  public:
    Dht22 (String sid, int pin) : Sensor(sid) {
      sensor = new DHT(pin, DHT_TYPE22);
      sensor->begin();
    }

    bool measure() {
      this->transmit(this->getSid() + 't', sensor->readTemperature());
      this->transmit(this->getSid() + 'h', sensor->readHumidity());
    }
    
};

class IRThermo : public Sensor {
  Adafruit_MLX90614 sensor;

  public:
    IRThermo (String sid, byte address) : Sensor(sid) {
      sensor = Adafruit_MLX90614(address);
      sensor.begin();
    }

    bool measure() {
      this->transmit(this->getSid() + 'l', sensor.readObjectTempC());
      this->transmit(this->getSid() + 'a', sensor.readAmbientTempC());
    }

};

class SoilTemps : public Sensor {
  int N;
  String *sidX;
  size_t *addr;
  OneWire oneWire;
  DallasTemperature sensors;
  

  public:
    SoilTemps (String sid, int num, int bus) : Sensor(sid) {
      this->N = num;
      this->sidX = new String[N];
      this->addr = new size_t[N];
      this->oneWire = OneWire(bus);
      this->sensors = DallasTemperature(&oneWire);
      this->sensors.begin();

      if (sensors.getDeviceCount() != N)
        this->printErr("At least one SoilTemp sensor is missing !!!");

    }

    void setAddr(int id, String sidX, DeviceAddress addr) {
      this->sidX[id] = sidX;
      this->addr[id] = addr;
      
      if (!sensors.isConnected(addr))
        this->printErr("Sensor not connected: " + sidX);
    }

    bool measure() {
      sensors.requestTemperatures();

      for (int i=0; i<N; i++)
        this->transmit(this->getSid() + sidX[i], sensors.getTempC(addr[i]));
        
    }

};

class SoilMoisture : public Sensor {
  int pin;

  public:
    SoilMoisture (String sid, int pin) : Sensor(sid) {
      this->pin = pin;
    }
    
    bool measure() {
      float value = 0;
      
      for (int i = 0; i < 100; i++)
          value = value + analogRead(pin);
    
      value = value / 100.0;
      this->transmit(value);

    }
};
