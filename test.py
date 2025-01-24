import time

class PIDController:
    def __init__(self, Kp, Ti, Td, target_glucose):
        self.Kp = Kp
        self.Ti = Ti if Ti != 0 else 1e-6  # Vermeidung von Division durch Null
        self.Td = Td
        self.target_glucose = target_glucose
        self.integral = 0
        self.previous_error = 0
        self.previous_time = time.time()

    def calculate(self, current_glucose):
        error = self.target_glucose - current_glucose
        proportional = self.Kp * error

        current_time = time.time()
        delta_time = current_time - self.previous_time

        # Integrale Komponente
        self.integral += error * delta_time
        integral = (self.Kp / self.Ti) * self.integral

        # Derivative Komponente
        derivative = 0
        if delta_time > 0:
            derivative = (self.Kp * self.Td) * ((error - self.previous_error) / delta_time)

        # Speichern für die nächste Berechnung
        self.previous_error = error
        self.previous_time = current_time

        return proportional + integral + derivative


class GlucoseSensor:
    def __init__(self, delay):
        self.delay = delay
        self.history = []

    def measure(self, real_glucose):
        self.history.append((time.time(), real_glucose))
        return self._get_delayed_glucose()

    def _get_delayed_glucose(self):
        current_time = time.time()
        for timestamp, glucose in reversed(self.history):
            if current_time - timestamp >= self.delay:
                return glucose
        # Rückgabe des ältesten Wertes, falls Verzögerung nicht erreicht
        return self.history[0][1] if self.history else None


if __name__ == "__main__":
    Kp = 0.1
    Ti = 60
    Td = 30
    target_glucose = 120
    sensor_delay = 10

    pid_controller = PIDController(Kp, Ti, Td, target_glucose)
    glucose_sensor = GlucoseSensor(sensor_delay)

    real_glucose = 180
    simulation_duration = 60
    start_time = time.time()

    while time.time() - start_time < simulation_duration:
        sensor_glucose = glucose_sensor.measure(real_glucose)

        if sensor_glucose is not None:
            insulin_dose = pid_controller.calculate(sensor_glucose)
            print(f"Echte Glukose: {real_glucose:.2f} mg/dl, Sensorglukose: {sensor_glucose:.2f} mg/dl, Insulinabgabe: {insulin_dose:.2f}")
        else:
            print("Keine Sensordaten verfügbar")

        real_glucose = max(40, real_glucose - (insulin_dose * 0.01))
        time.sleep(1)
