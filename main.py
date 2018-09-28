import os, glob, pyexcel
import matplotlib.pyplot as plt
from numpy import *
from scipy.optimize import curve_fit
from scipy.signal import medfilt, butter, filtfilt, freqz, get_window
from scipy.stats

# Read Me
# Specify the main directory below under workDir
# Then place data in the specific order inside the main directory : Sensor/Scenario/execfile
# For example: Strain_Gauge/SG_dyn/dv1.xlsx

# Main working directory
workDir = 'c:/Users/Richard/Documents/GitHub/ME103_Lab_2'
# Change to main folder
os.chdir(workDir)
print('Now in {}'.format(workDir))

# Check all valid sensor folder
senList = [sen for sen in os.listdir(workDir) if os.path.isdir(sen)]
senList.remove('.git')
senList.remove('__pycache__')
report = dict.fromkeys(senList)
for sensor in senList:
    senDir = workDir + '/' + sensor
    os.chdir(senDir)
    print('Now in {}'.format(senDir))
    # Sensor Summary Template
    sensDict = {'static': [], 'quasi_static': [], 'transient': [],
                'dynamic': [], 'mean': []}
    # All Sensor summary list
    report[sensor] = sensDict

    # Check all valid scenrio folders
    dirList = [dir for dir in os.listdir(senDir) if os.path.isdir(dir)]
    for folderName in dirList:
        # Change to scenario folders
        os.chdir(senDir + '/' + folderName)
        print('Now in {}'.format(senDir+'/'+folderName))

        residual = zeros((3, 50000))
        xData = zeros((3, 50000))

        # Iterating through 3 trials
        for trial in range(1,4):
            print('Reading trial {}'.format(trial))
            trialList = ['td'+str(trial)+'.xlsx', 'ta'+str(trial)+'.xlsx']
            # Total Displacement from td.xlsx
            # Time from ta.xlsx
            # Voltage from ta.xlsx

            # Low pass filter
            def low_pass_filter(x, cutoff = 124, sampling = 10000, order = 6):
                nyquist = sampling/2
                normalize_freq = (cutoff/nyquist)
                b, a = butter(order, normalize_freq,
                            btype = 'lowpass', analog = False)
                # Filter Frequency Response
                # w, h = freqz(b, a)
                # plt.semilogx(w/(2*pi)*sampling, 20*log10(abs(h)))
                # plt.axvline(12.4, color='k')
                # plt.grid(which='both', axis='both')
                # plt.title('Butterworth filter frequency response')
                # plt.xlabel('Frequency [Hz]')
                # plt.ylabel('Amplitude [dB]')
                # plt.show()


                x[:] = filtfilt(b, a, x)
                return x
                # return filtfilt(b, a, x)

            # Iterating throught data files to collect data
            for fileName in trialList:
                # Getting sheets from file
                sheet = pyexcel.get_sheet(file_name = fileName)
                # Converting data in sheet to array (50001x2)
                rawData = sheet.to_array()
                # Removing axis entries
                del rawData[0]
                # Transform to 2x50000 matrix
                data = array(rawData, dtype = float).transpose()
                # Selectively collecting data from each file
                if fileName == trialList[0]:
                    time = data[0]
                    totalDisplacement = data[1]*0.0254

                elif fileName == trialList[1]:
                    voltage = low_pass_filter(medfilt(medfilt(data[1])))

                else:
                    raise ValueError('Error reading the wrong trial file')

            def rezeroing():
                global totalDisplacement
                shift = average(totalDisplacement)
                totalDisplacement-=shift

            def truncate():
                global time, voltage, totalDisplacement, sensor
                indexBool = time>=1
                voltage = voltage[indexBool]
                time = time[indexBool]
                if sensor == 'Accelerometer':
                    global acceleration
                    acceleration = acceleration[indexBool]
                elif sensor == 'Voice_Coil':
                    global velocity
                    velocity = velocity[indexBool]
                else:
                    totalDisplacement = totalDisplacement[indexBool]
                return indexBool

            def numerical_diff(x, dt, scale):
                global low_pass_filter
                xdot = zeros(len(x))
                xdotdot = zeros(len(x))
                for i in range(len(x)):
                    if i <= 1*scale:
                        # Forward differentiation
                        xdot[i] = (-3*x[i] + 4*x[i+1*scale]
                                    - x[i+2*scale])/(2*scale*dt)
                        xdotdot[i] = (x[i+2*scale] - 2*x[i+1*scale]
                                    + x[i])/((scale*dt)**2)
                    elif i >= (len(x) - 2*scale):
                        # Backward differentiation
                        xdot[i] = (3*x[i] - 4*x[i-1*scale]
                                    + x[i-2*scale])/(2*scale*dt)
                        xdotdot[i] = (-x[i-2*scale] + 2*x[i-1*scale]
                                    - x[i])/((scale*dt)**2)
                    else:
                        # Fourth order center differentiation
                        xdot[i] = (-x[i+2*scale]
                                    + 8*x[i+1*scale] - 8*x[i-1*scale]
                                    + x[i-2*scale])/(12*scale*dt)
                        xdotdot[i] = (-x[i+2*scale] + 16*x[i+1*scale] -30*x[i]
                                    + 16*x[i-1*scale]
                                    - x[i-2*scale])/(12*(scale*dt)**2)

                return low_pass_filter(medfilt(xdot)), xdotdot

            # Curve fitting based on scenario
            def checkScenario(f, x, y):
                global report, trial, sensor, residual, xData

                if 'static' in folderName:
                    y = low_pass_filter(y, cutoff = 10)
                    # y = average(x)*ones(len(x))
                    popt, _ = curve_fit(f, x, y)
                    report[sensor]['static'].append(popt)

                elif 'qui' in folderName:
                    y = low_pass_filter(y, cutoff = 20)
                    popt, _ = curve_fit(f, x, y)
                    report[sensor]['quasi_static'].append(popt)

                elif 'trans' in folderName:
                    index = truncate()
                    if len(residual[0]) == 50000:
                        residual = zeros((3, len(x[index])))
                        xData = zeros((3, len(x[index])))

                    popt, _ = curve_fit(f, x[index], y[index])
                    report[sensor]['transient'].append(popt)

                elif 'dyn' in folderName:
                    rezeroing()
                    popt, _ = curve_fit(f, x, y)
                    report[sensor]['dynamic'].append(popt)

                else:
                    raise ValueError('Scenario Not Found ')

                return popt

            # Plotting
            def plottingSubplot3(x1, y1, x2, y2, x3, y3, x4, y4,
                                xLabels, yLabels, labels):

                global sensor, folderName, trial
                plt.figure(figsize = (14, 7))

                plt.subplot(311)
                plt.gca().set_title('Sensor: {} Case : {} Trial {}'.format(
                sensor, folderName, str(trial)))
                plt.plot(x1, y1, label = labels[0])
                plt.xlabel(xLabels[0])
                plt.ylabel(yLabels[0])
                plt.grid(b = True, which = 'both')
                plt.legend()

                plt.subplot(312)
                plt.plot(x2, y2, label = labels[1])
                plt.xlabel(xLabels[1])
                plt.ylabel(yLabels[1])
                plt.grid(b = True, which = 'both')
                plt.legend()

                plt.subplot(313)
                plt.plot(x3, y3, label = labels[2])
                plt.plot(x4, y4, 'r', label = 'Fitted Curve')
                plt.xlabel(xLabels[2])
                plt.ylabel(yLabels[2])
                plt.grid(b = True, which = 'both')
                plt.legend()

                plt.savefig('Trial_{}.png'.format(trial))
                # plt.show()
                plt.close()

            # Data Anaysis
            if sensor == 'Strain_Gauge':
                # Governing equation
                # Vout = Vin + g*IR*delta_L + I*R
                sg = lambda x, a, b: 5 + 2*a*x + b
                coefficient = checkScenario(sg, totalDisplacement, voltage)

                fit = sg(totalDisplacement, coefficient[0], coefficient[1])

                yBar = sum(voltage)/len(voltage)
                ssreg = sum((fit - yBar)**2)
                sstot = sum((voltage - yBar)**2)

                print('R square is : {}'.format(ssreg/sstot))

                residual[trial-1] = voltage - fit

                xData[trial-1] = totalDisplacement

                plottingSubplot3(time, voltage, time, totalDisplacement,
                            totalDisplacement, voltage,
                            totalDisplacement, fit,
                            ['Time (s)', 'Time (s)', 'Distance (m)'],
                            ['Voltage (v)', 'Distance (m)', 'Voltage (v)'],
                            [sensor, 'Encoder', 'data'])

            elif sensor == 'IR':
                # Governing equation
                ir = lambda x, a, b: a*x + b
                coefficient = checkScenario(ir, totalDisplacement, voltage)

                fit = ir(totalDisplacement, coefficient[0], coefficient[1])

                yBar = sum(voltage)/len(voltage)
                ssreg = sum((fit - yBar)**2)
                sstot = sum((voltage - yBar)**2)

                print('R square is : {}'.format(ssreg/sstot))

                residual[trial-1] = voltage - fit

                xData[trial-1] = totalDisplacement

                plottingSubplot3(time, voltage, time, totalDisplacement,
                            totalDisplacement, voltage,
                            totalDisplacement, fit,
                            ['Time (s)', 'Time (s)', 'Distance (m)'],
                            ['Voltage (v)', 'Distance (m)', 'Voltage (v)'],
                            [sensor, 'Encoder', 'data'])

            elif sensor == 'Voice_Coil':
                velocity, _ = numerical_diff(totalDisplacement,
                                            dt = time[1]-time[0], scale = 40)
                # Governing Equation
                vc = lambda x, a, b: a*x + b
                coefficient = checkScenario(vc, velocity, voltage)

                fit = vc(velocity, coefficient[0], coefficient[1])

                yBar = sum(voltage)/len(voltage)
                ssreg = sum((fit - yBar)**2)
                sstot = sum((voltage - yBar)**2)

                print('R square is : {}'.format(ssreg/sstot))

                residual[trial-1] = voltage -fit

                xData[trial-1] = velocity

                plottingSubplot3(time, voltage, time, velocity,
                            velocity, voltage, velocity, fit,
                            ['Time (s)', 'Time (s)', 'Velocity (m/s)'],
                            ['Voltage (v)', 'Velocity (m/s)', 'Voltage (v)'],
                            [sensor, 'Encoder', 'data'])

            elif sensor == 'Accelerometer':
                velocity, _ = numerical_diff(totalDisplacement,
                                            dt = time[1]-time[0], scale = 40)
                acceleration, _ = numerical_diff(velocity,
                                    dt = time[1]-time[0], scale = 40)
                # Governing Equation
                acc = lambda x, a, b: a*x + b
                coefficient = checkScenario(acc, acceleration, voltage)

                fit = acc(acceleration, coefficient[0], coefficient[1])

                yBar = sum(voltage)/len(voltage)
                ssreg = sum((fit - yBar)**2)
                sstot = sum((voltage - yBar)**2)

                print('R square is : {}'.format(ssreg/sstot))

                residual[trial-1] = voltage - fit

                xData[trial-1] = acceleration

                plottingSubplot3(time, voltage, time, acceleration,
                        acceleration, voltage, acceleration, fit,
                        ['Time (s)', 'Time (s)', 'Acceleration (m s^-2)'],
                        ['Voltage (v)', 'Acceleration (m s^-2)', 'Voltage (v)'],
                        [sensor, 'Encoder', 'data'])

            else:
                raise ValueError('Sensor Not Regonized')

        plt.figure(figsize = (14, 7))
        plt.title('Residual : {}'.format(folderName))
        plt.plot(xData[0], residual[0], ':', label = 'Trial 1')
        plt.plot(xData[1], residual[1], ':', label = 'Trial 2')
        plt.plot(xData[2], residual[2], ':', label = 'Trial 3')
        plt.grid(b = True, which = 'both')
        plt.legend()
        plt.savefig('Residual.png')
        # plt.show()
        plt.close()

    print(report[sensor])
